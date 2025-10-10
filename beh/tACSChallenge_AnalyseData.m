function [all_ps,all_bs,all_hit_probs,all_bs_perm,all_pfs,group_level_p,group_level_comp,pf_p] = tACSChallenge_AnalyseData(data_path, subjs, conditions, perm)
%% script originally written by Benedikt Zoefel, CNRS Toulouse, in October 2021
%% modified in April 22 and June 25 (minor fixes)
%% added permutations and preferred phases in April 24

%% run this script with the following input argments:
%% - folder to your data, example: '../data/'
%% - subect initials to be analysed. The data for each subject must be in a separate folder that is labeled as such
%% (e.g., 'P01') and is located in the data folder. Example: {'P01, P02', P03'}
%% - condition labels to be analysed. Label must be part of the file name (including the stars in the example).
%% -"perm" is the number of permutations (typically 1000)
%% Example: {'*Montage A*','*Montage B*','*Montage C*'}

%% requires circular statistics toolbox

clc; close all;

addpath('./functions/');

no_phases = 8;

all_ps = zeros(length(conditions),length(subjs));
all_bs = zeros(length(conditions),length(subjs));
all_hit_probs = zeros(no_phases,length(conditions),length(subjs));
all_bs_perm = zeros(length(conditions),length(subjs),perm);
all_pfs = zeros(length(conditions),length(subjs)); % preferred phases

for s = 1:length(subjs)
    % load the data
    curr_data = tACSChallenge_SortData(data_path, subjs{s}, conditions);
    % and analyse it
    [all_ps(:,s), all_bs(:,s), all_hit_probs(:,:,s),all_bs_perm(:,s,:),all_pfs(:,s)] = tACSChallenge_EvalData(curr_data,perm,s); 
    
end

all_hit_probs(no_phases+1,:,:) = all_hit_probs(1,:,:); % duplicate first phase bin for visualisation

%% plots detection probability as a function of binned phases (only for visualisation of results)
figure
plot(-pi:pi/4:pi,squeeze(mean(all_hit_probs,3)))
xlabel('tACS phase'); ylabel('detection probability'); legend(conditions);

%% plot phasic modulation of target detection for each condition
figure
bar(mean(all_bs,2))
xlabel('conditions'); xticklabels(conditions); ylabel('modulation strength');

group_level_p = zeros(length(conditions),1);
group_level_comp = zeros(2,1); % A vs B, A vs C
pf_p = zeros(2,1);

%% group level statistics
%% test for phasic effect separately for each condition
for cond = 1:length(conditions)
    % observed mean
    currmean = mean(all_bs(cond,:),2);
    % average across subjects in surrogate
    currsurr = squeeze(mean(all_bs_perm(cond,:,:),2));
    % mean of surrogate
    currsurrmean = mean(currsurr);
    % std of surrogate
    currsurrstd = std(currsurr);
    % z-score
    [~,group_level_p(cond)] = ztest(currmean,currsurrmean,currsurrstd,'tail','right');
end
%% contrasting conditions, assuming that the tACS condition is the first one
[~,group_level_comp(1)] = ttest(all_bs(1,:),all_bs(2,:),'tail','right');
[~,group_level_comp(2)] = ttest(all_bs(1,:),all_bs(3,:),'tail','right');

%% compare difference in preferred phase against 0
pf_p(1) = circ_mtest(circ_dist(all_pfs(1,:),all_pfs(2,:)),0);
pf_p(2) = circ_mtest(circ_dist(all_pfs(1,:),all_pfs(3,:)),0);