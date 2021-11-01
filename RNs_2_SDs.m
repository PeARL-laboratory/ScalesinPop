function [tab,fig] = RNs_2_SDs(chords,scores,varargin)
%
% [tab,fig] = RNs_2_SDs(chords,scores,varargin)
%
% Creates a table, bar plot, and correlation for the proportion of scale
% degrees that appear in the list of chords and their topic model scores.
% Users may specify weights for the root, bass, and other chord members in
% each chord, where the sum of the weights should equal 1 (default = 1/3
% for each).
%
% Input arguments:
%  chords         (required) -- An m-dimensional string array of Roman 
%                               numeral annotations (e.g., V43). Only
%                               extracts scale-degrees for triads and
%                               seventh chords (i.e., ignores other
%                               extensions, additions, etc.).
%
%  scores         (required) -- An m-dimensional numeric array of LDA
%                               scores. 
%
%  varargin       (optional) -- Weights for the root, bass, and other chord
%                               members reflected in each RN annotation.
%                               Default = 1/3 for each. 
%
% Output:
%  tab      -- Tables containing distributions for 35 (7 roots x 5
%              accidentals) and 12 possible scale degrees.
%  fig      -- Bar plots of SD distributions, and best-fitting KS major- or 
%              minor-key profile shown as a line plot.
%
% Example:
% [tab,fig] = RNs_2_SDs(chords,scores,[1/4 1/4 1/2]);
%
%          
%
% Change History :
% Date          Time	Prog            Note
% 11.01.2021  	12:26	David Sears	    Created under MATLAB 9.8 (R2020a, PC)

%% CHECK OPTIONAL ARGUMENT

if nargin > 2
    weights = varargin{1};
else
    weights = [1/3 1/3 1/3];
end


%% EXTRACT RN COMPONENTS

% Extract RN root.
idx = regexp(chords,'[vViI]');
sz = size(chords);
roots = strings(sz);
for i = 1:length(idx)
    tmp = char(chords(i));
    tmp_idx = cell2mat(idx(i));
    if tmp_idx(1) ~= 1
        tmp_idx = [1:tmp_idx(1)-1 tmp_idx];
    end
    roots(i,1) = string(tmp(tmp_idx));
end
roots = upper(roots);

% Extract quality.
quals(1:length(chords),1) = {'maj'};
idx = contains(chords,{'i' 'v'}); 
quals(idx) = {'min'};
idx = contains(chords,{'o'});
quals(idx) = {'dim'};
idx = contains(chords,{'h'});
quals(idx) = {'hdim'};
idx = contains(chords,{'+'});
quals(idx) = {'aug'};
idx = contains(chords,{'p'});
quals(idx) = {'pow'};

% Extract inversion.
idx = regexp(chords,'[^vViI]');
inversions = strings(sz);
for i = 1:length(idx)
    tmp = char(chords(i));
    tmp_idx = cell2mat(idx(i));
    if ~isempty(tmp_idx)
        inversions(i,1) = string(tmp(tmp_idx));
    end
end
inversions = erase(inversions,{'-' '#' 'sus4' 'sus2' 'Ger' 'Fr' 'It' 'Ct' 'p' 'd' 'M'});
idx = find(contains(chords,'('));
inversions(idx) = extractBefore(inversions(idx),'(');

% Extract extensions.
extensions = strings(sz);
idx = find(contains(chords,'('));
extensions(idx) = extractBetween(chords(idx),'(',')');
idx = find(contains(chords,'sus'));
extensions(idx) = strcat(extensions(idx),extractAfter(chords(idx),'sus'));


%% CONVERT CHORD MEMBERS TO SCALE DEGREES

RNs = {'--IV' '--I' '--V' '--II' '--VI' '--III' '--VII' ...
       '-IV' '-I' '-V' '-II' '-VI' '-III' '-VII' ...
       'IV' 'I' 'V' 'II' 'VI' 'III' 'VII' ...
       '#IV' '#I' '#V' '#II' '#VI' '#III' '#VII' ...
       '##IV' '##I' '##V' '##II' '##VI' '##III' '##VII'};
SDs = {'--4' '--1' '--5' '--2' '--6' '--3' '--7' ... 
      '-4' '-1' '-5' '-2' '-6' '-3' '-7' ...
      '4' '1' '5' '2' '6' '3' '7' ...
      '#4' '#1' '#5' '#2' '#6' '#3' '#7' ...
      '##4' '##1' '##5' '##2' '##6' '##3' '##7'};
Quals = {'aug' 'dim' 'hdim' 'maj' 'min' 'pow'};
members = cell(length(chords),6);

% Roots to SDs
[~,roots_idx] = ismember(roots,RNs);
SD_roots = SDs(roots_idx);
SD_roots = SD_roots';

% Qualities to SDs
qual_IDs = [4 8; -3 -6; -3 -6; 4 1; -3 1; 1 nan];
[~,qual_idx] = ismember(quals,Quals);
qual_num = qual_IDs(qual_idx,:);
quals_idx = roots_idx+qual_num;
SD_quals = cell(length(chords),2);
for i = 1:length(chords)
    tmp = quals_idx(i,:);
    if any(isnan(tmp))
        tmp(isnan(tmp)) = [];
        SD_quals(i,1) = SDs(tmp);
    else
    SD_quals(i,:) = SDs(tmp);
    end
    clear tmp
end

% Inversions to SDs
SDs_tot = [SD_roots SD_quals];
SDs_tot(:,4) = {'NaN'};
SDs_tot(find(cellfun(@isempty,SDs_tot))) = {'NaN'}; % for power chords, only two scale degrees.

% seventh? -- use quality of triad to determine quality of 7th, except for V
sev_IDs = [5 -9 -2 5 -2 -2]; % aug7, dim7, hdim7, maj7, min7, pow
sevenths = {'9' '7' '65' '43' '42'};
excl_idx = ismember(inversions,sevenths);
sev_num = [sev_IDs(qual_idx)]';
idx = strcmp(roots,'V'); % V is an exception; it requires a m7.
sev_num(idx) = -2;
idx = contains(chords,'d'); % d also requires a m7.
sev_num(idx) = -2;
idx = contains(chords,'M'); % d also requires a m7.
sev_num(idx) = 5;
sev_idx = roots_idx+sev_num;
SDs_tot(:,4) = SDs(sev_idx);
SDs_tot(excl_idx==0,4) = {'NaN'};

% identify bass
SD_bass = SD_roots;
idx = ismember(inversions,{'6' '65'}); % third
SD_bass(idx) = SDs_tot(idx,2);
idx = ismember(inversions,{'64' '43'}); % third
SD_bass(idx) = SDs_tot(idx,3);
idx = ismember(inversions,'42'); % third
SD_bass(idx) = SDs_tot(idx,3);

% Extensions to SDs -- Let's just deal with the cadential six-four
idx = strcmp(chords,'V(64)');
if any(idx==1)
    SDs_tot(idx,2:3) = {'1' '3'};
end


%% TABLE and PLOT(S) of SCALE-DEGREE DISTRIBUTIONS

% Sum weighted scores for scale-degrees. 
SDs_tot = SDs_tot';
SD = SDs_tot;
SD = SD(:);
scores = repelem(scores,4);
scores(strcmp(SD,'NaN')) = 0;
weight = repelem(weights(3),length(scores),1);
weight(1:4:end) = weights(1); % root
for i = 1:size(SDs_tot,2)
    [~,idx2(1,i)] = nanismember(SD_bass(i),SDs_tot(:,i));
end
weight_idx = sub2ind(size(SDs_tot),idx2,[1:i]);
weight(weight_idx) = weight(weight_idx) + weights(2); % bass
scores = scores.*weight;
idx = strcmp(SD,'NaN');
scores(idx) = [];
SD(idx) = [];

% Create table with 35 scale degrees (7 roots x 5 accidentals).
SDs_sorted = {'--1' '-1' '1' '#1' '##1' ...
              '--2' '-2' '2' '#2' '##2' ...
              '--3' '-3' '3' '#3' '##3' ...
              '--4' '-4' '4' '#4' '##4' ...
              '--5' '-5' '5' '#5' '##5' ...
              '--6' '-6' '6' '#6' '##6' ...
              '--7' '-7' '7' '#7' '##7'};
[~,idx] = ismember(SD,SDs_sorted);
probs = accumarray(idx,scores);
probs(max(idx)+1:length(SDs_sorted)) = 0;
probs = probs/sum(probs); % normalize
labels = SDs_sorted';
tab35 = table(labels,probs); 

% Create table with 12 scale degrees.
tet12_idx = [3 7 8 12 13 18 19 23 27 28 32 33];
tab12 = tab35(tet12_idx,:);
labs = [{'1' '#1/b2' '2' '#2/b3' '3' '4' '#4/b5' '5' '#5/b6' '6' '#6/b7' '7'}]';
tab12(:,1) = labs;

% Aggregate tables.
tab.tab35 = tab35;
tab.tab12 = tab12;

% Correlate probabilities to KS major- and minor-key profiles. (Install
% MIDIToolbox.)
kkmaj = refstat('kkmaj');
kkmin = refstat('kkmin');
kkmaj = kkmaj/sum(kkmaj);
kkmin = kkmin/sum(kkmin);
[r p] = corr(tab.tab12.probs,kkmaj');
tab.corr.maj = table(r,p);
[r p] = corr(tab.tab12.probs,kkmin');
tab.corr.min = table(r,p);

% Create plot.
fig = bar(tab12.probs);
hold on
if tab.corr.maj.r > tab.corr.min.r
    fig = plot(kkmaj,'k');
    legend(fig,'kkmaj');
else
    fig = plot(kkmin,'k');
    legend(fig,'kkmin');
end
set(gca,'xtick',1:length(labs),'xticklabel',labs);
xlabel('Scale','fontweight','bold');
ylabel('Proportion','fontweight','bold');


end

