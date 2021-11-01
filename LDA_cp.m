% Split the common-practice corpora (N=1000) into 50/50 train/test, 
% stratified by data set. The train data will then be divided into 10 folds
% stratified by data set. 

tic

%% IMPORT & CROSS-VALIDATION

% import
load datasets.mat

% train and test folds (50/50), stratified by dataset
numDocuments = numel(cp_excerpts.docs);
cvp = cvpartition(cp_excerpts.groups,'HoldOut',0.5);
train.docs = cp_excerpts.docs(cvp.training);
train.groups = cp_excerpts.groups(cvp.training);
test.docs = cp_excerpts.docs(cvp.test); 
test.modes = cp_excerpts.modes(cvp.test); 


%% TRAIN TOPIC MODEL -- N TOPICS

% use train data set to determine number of topics using 10-fold cross-validation
k=10;
cvp_train = cvpartition(train.groups,'KFold',k);
for i = 1:k
    
    % get partitions
    train_train = train.docs(cvp_train.training(i));
    train_test = train.docs(cvp_train.test(i));

    % bag of words
    bag = bagOfWords(train_train);
    bag = removeInfrequentWords(bag,5);
    bag = removeEmptyDocuments(bag);
    
    % calculate perplexity for range of topics
    numTopicsRange = [1:10];
    for j = 1:numel(numTopicsRange)
        numTopics = numTopicsRange(j);
        mdl = fitlda(bag,numTopics,'Solver','cgs','Verbose',0);
        [~,perplexity(i,j)] = logp(mdl,train_test);
        clear numTopics mdl
    end
    clear documentsTrain documentsValidation bag
end
figure; errorbar(mean(perplexity),std(perplexity)/sqrt(10));
ylabel("Perplexity (±1 SEM)")
xlabel("Number of Topics")

% use t-test with Bonferroni correction to determine last N that
% significantly differs from N-1.
for i = 2:numTopicsRange(end)
    [~,p(i-1,1),~,stats] = ttest2(perplexity(:,i),perplexity(:,i-1));
    pp(i-1,1) = p(i-1)*(i-1);
    t(i-1,1) = stats.tstat;
end
numTopics = find(p<.05,1,'last')+1; % Bonferroni correction for up to N topics.


%% TEST TOPIC MODEL -- EVALUATION

% bag of words
bag = bagOfWords(test.docs);
bag = removeInfrequentWords(bag,5);
bag = removeEmptyDocuments(bag);

% fit model on test data
mdl_test = fitlda(bag,numTopics,'Verbose',0);

% label predictions and observations.
pred = zeros(length(mdl_test.DocumentTopicProbabilities),1);
pred(mdl_test.DocumentTopicProbabilities(:,1)>mdl_test.DocumentTopicProbabilities(:,2)) = 1;
pred = logical(pred);
top = topkwords(mdl_test,10,1); % which topic probably represents the major mode?
if any(strcmp(top.Word,'I'))
    obs = strcmp(test.modes,'major');
else
    obs = strcmp(test.modes,'minor');
end

% calculate confusion matrix and stats
stats = confusionmatStats(obs,pred);

% Calculate balanced F measure weighted by number of major and minor mode movements.
WeightedFscore = (sum(strcmp(test.modes,'major')) * stats.Fscore(1) + sum(strcmp(test.modes,'minor')) * stats.Fscore(2)) / length(test.modes);


%% TRAIN SCALE FUNCTION -- WEIGHTING ROOT, BASS, OTHER SDs

% bag of words
bag = bagOfWords(train.docs);
bag = removeInfrequentWords(bag,2);
bag = removeEmptyDocuments(bag);

% fit model on test data
mdl_train = fitlda(bag,numTopics,'Verbose',0);

% create weights
allowed_values = .02 : .02 : 1;
num_val = length(allowed_values);
weights = [];
for P1_idx = 1 : num_val
   P1 = allowed_values(P1_idx);
   last_P2 = find(allowed_values <= 1 - P1, 1, 'last');
   for P2_idx = 1 : last_P2
     P2 = allowed_values(P2_idx);
     P3 = 1 - P1 - P2;
     weights = [weights; P1 P2 P3];
     clear P2 P3
   end
   clear P1 last_P2
end

% Run model for all possible weights.
numChords = 20;
maj = table;
min = table;
for j = 1:numTopics
    
    % find top words for topic i.
    top = topkwords(mdl_train,numChords,j);
        
    for i = 1:size(weights,1)
         
        % get scale stats.
        tmp = SDs_dist_v2(top.Word,top.Score,weights(i,:));
        root = weights(i,1);
        bass = weights(i,2);
        other = weights(i,3);
        if any(strcmp(top.Word,'I'))
            maj = [maj; addvars(tmp.corr.maj,root,bass,other)];
        else
            min = [min; addvars(tmp.corr.min,root,bass,other)];
        end
        clear tmp root bass other
        i
    end
    clear top
end
maj.Properties.VariableNames = {'r' 'p' 'root' 'bass' 'other'};
min.Properties.VariableNames = {'r' 'p' 'root' 'bass' 'other'};
maj = sortrows(maj,'r','descend');
min = sortrows(min,'r','descend');


%% TEST SCALE FUNCTION

fig_cloud = figure;
fig_scale = figure;
for i = 1:numTopics
    
    % word cloud
    set(0,'CurrentFigure',fig_cloud);
    if rem(numTopics,2)==0
        subplot(2,numTopics/2,i);
    else
        subplot(1,numTopics,i);
    end
    wordcloud(mdl_test,i);
    title("Topic " + i);
    
    % create table.
    name = strcat('topic_',num2str(i));
    tab.(name).top = topkwords(mdl_test,numChords,i);
    
    % create scale figure.
    set(0,'CurrentFigure',fig_scale);
    if rem(numTopics,2)==0
        subplot(2,numTopics/2,i);
    else
        subplot(1,numTopics,i);
    end
    if any(strcmp(tab.(name).top.Word,'I'))
        tab.(name).scale = RNs_2_SDs(tab.(name).top.Word,tab.(name).top.Score,table2array(maj(1,3:5)));
    else
        tab.(name).scale = RNs_2_SDs(tab.(name).top.Word,tab.(name).top.Score,table2array(min(1,3:5)));
    end
    title("Topic " + i);
    
    clear name
end


%% CLEAN & EXPORT
clearvars -except maj min mdl_train mdl_test weights stats WeightedFscore train test numTopics cvp t p pp perplexity tab fig_cloud fig_scale
save LDA_cp.mat
toc