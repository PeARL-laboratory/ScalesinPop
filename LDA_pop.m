% Split the popular corpora (N=938) into 50/50 train/test, stratified by 
% data set. The train data will then be divided into 10 folds stratified
% by data set.

tic

%% IMPORT & CROSS-VALIDATION

% import
load('datasets.mat','pop')

% train and test folds (50/50), stratified by dataset
cvp = cvpartition(pop.groups,'HoldOut',0.5);
train.groups = pop.groups(cvp.training);
train.docs = pop.docs(cvp.training);
test = pop.docs(cvp.test); 


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
numTopics = find(pp<.05,1,'last')+1; % Bonferroni correction for up to N topics.


%% TEST TOPIC MODEL -- EXPLORE!

% import weights
load('LDA_cp.mat','maj','min')

% bag of words
bag = bagOfWords(test);
bag = removeInfrequentWords(bag,5);
bag = removeEmptyDocuments(bag);

% fit model on test data
mdl_test = fitlda(bag,numTopics,'Verbose',0);

% figures for topics and scales.
fig_cloud = figure;
fig_scale = figure;
numChords = 20;
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
save LDA_pop.mat
toc