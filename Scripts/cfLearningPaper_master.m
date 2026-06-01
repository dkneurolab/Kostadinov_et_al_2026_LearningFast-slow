%% Master file for figures of Kostadinov and Hausser, 2022
set(0,'DefaultFigureWindowStyle','docked');

% Paths to data
if ismac
    dataFold     = '/Users/dimitar/Desktop/CF_learning_paper'; % Where the raw data lives
    paperFold      = '/Users/dimitar/Dropbox/DK papers/2022_Kostadinov_FastSlowLearning'; % Where the figures and data are saved
elseif ispc
    dataFold     = 'E:\Analysis_current'; % Where the raw data lives
    paperFold      = 'C:\Users\Dimitar\Dropbox\DK papers\2022_Kostadinov_FastSlowLearning'; % Where the figures and data are saved
end

% Load in data summaries and set parameters
imDatasets = load(fullfile(dataFold,'imDatasets.mat'));
behavDatasets = load(fullfile(dataFold,'behavDatasets.mat'));
behavData_v1 = load(fullfile(dataFold,'sessionsv1_behav.mat')); behavData_v1 = behavData_v1.sessionsv1_Behav;
behavData_v4_1 = load(fullfile(dataFold,'sessionsv4_learning.mat')); behavData_v4_1 = behavData_v4_1.sessionsv4_learning;
behavData_v4_n = load(fullfile(dataFold,'sessionsv4_n.mat')); behavData_v4_n = behavData_v4_n.sessionsv4_n;
behavData_v5 = load(fullfile(dataFold,'sessionsv5.mat')); behavData_v5 = behavData_v5.sessionsv5;
matchColors = distinguishable_colors(100,[0 0 0; 1 1 1]); matchColors(1:6,:) = matchColors([6,5,2,3,1,4],:); matchColors = [.6 .6 .6; matchColors];

allDatasets = behavDatasets;
imNames = fieldnames(imDatasets);
for i = 1:numel(imNames); allDatasets.(imNames{i}) = imDatasets.(imNames{i}); end; clear i imNames

% Parse behavioral datasets
mouseStructs = fieldnames(allDatasets);
sessions_Behav = [];
for i = 1:numel(mouseStructs)
    seshLocal = allDatasets.(mouseStructs{i});
    if ~isfield(seshLocal,'fpblock')
        seshLocal(1).fpblock = [];
    end
    sessions_Behav = [sessions_Behav seshLocal]; %#ok<AGROW>
end

% Parse imaging datasets by session type
mouseStructs = fieldnames(imDatasets);
sessions_im = [];
for i = 1:numel(mouseStructs)
    sessions_im = [sessions_im imDatasets.(mouseStructs{i})]; %#ok<AGROW>
end 
sessionsv1_im = sessions_im(cat(1,sessions_im.version) == 1);
sessionsv4_im = sessions_im(cat(1,sessions_im.version) == 4);
sessionsv5_im = sessions_im(cat(1,sessions_im.version) == 5);


% Local settings for sessions:
sessionsv4_n    = sessionsv4_im([2,3,5,6,8,10,12,14,16,18,20,22]);
sessionsv4_1    = sessionsv4_im([1,4,7,9,11,13,15,17,19,21]);
pcaVar          = 80; % Percent of PCA variance to be used in for design matrix reduction

% Set session parameters:
savebool        = false; % Save data, etc.
closebool       = true; % Close figures after running scripts?

% GLM params:
% GLM parameters:
glmParams.dataType      = 'spk'; % 'fluo' or 'spk' for zF or spike events
glmParams.distr         = 'gaussian'; % 'gaussian' (for fluo data), 'poisson' (for spk counts), or 'binomial' (for binarized spk trains)
glmParams.link          = []; % Use default link function for given distribution
glmParams.alpha         = .5; % Others use 0.5 for gaussian GLMs and 0.8-0.9 for poisson or binomial GLMs
glmParams.nLambda       = 100; % Good number, don't change from 100

%% -------------------FIGURE 1: Setting up behaviour-------------------
%% Figure 1 - alignment and controls

% Old style figure showing endpoints, etc.
adaptTrialGrp = 30;
fig1_percMovCorrect(sessions_Behav,dataFold,paperFold,adaptTrialGrp,matchColors,closebool,savebool);

% Extract behavioural stability for v1 data
seshV1_movCorr = fig1_v1stability(behavData_v1,dataFold,paperFold,savebool);
% Plot v1 and v4 learning curves
nLearningSesh = 6;
seshV4_movCorr = fig1_v4learning(seshV1_movCorr,behavData_v4_1,behavData_v4_n,nLearningSesh,dataFold,paperFold,savebool);
% Plot v1-v4 example data and summary
fig1_v1to4learning(seshV1_movCorr, seshV4_movCorr,nLearningSesh+2,dataFold,paperFold,savebool);

% Plot v5 adaptation curves
adaptTrialGrp = 30;
fig1_adaptBehav(behavData_v5,adaptTrialGrp,dataFold,paperFold,savebool)

% Make example behav traces using DK169 naive, learning, expert, and adaptation sessions
fig1_exampleBehav(sessions_im(13:16),dataFold,paperFold,savebool)

% Manual data entry: days between v1v3end and v41v4n:
v1v3endDays = [25 21 31 31 16 11 20 20];
v41v4nDays = [6 12 8 9 7 13 22 15];

    
% DK070 138 103 105 169 171 194 199

%% -------------------FIGURE 2: Chronic imaging-------------------

% Calculate how many ROIs can be matched across days
matchROIs = fig2_matchedROIs_vs_days(dataFold,paperFold,matchColors,savebool);

% Calculate baseline FRs for matched ROIs and plot their onset, movement, and reward signals
[matchData, matchFOVs] = fig2_roiMatch(imDatasets,matchROIs,dataFold,paperFold,matchColors, savebool);
[matchSpont] = fig2_roiMatchSpont(matchData,dataFold,paperFold,matchColors, savebool);

% Combine v144 sessions into mega-response-matrix
matchPCs_v144 = fig2_match_v144(matchData,matchSpont,matchFOVs,dataFold,paperFold,matchColors,savebool);
% % Combine v145 sessions into mega-response-matrix
% matchPCs_v145 = fig1_match_v145(matchData,matchSpont,matchFOVs,dataFold,paperFold,matchColors,savebool);
% 
matchStats144 = fig2_matchStats(matchData, matchPCs_v144, dataFold, paperFold,'v144',savebool);
% matchStats145 = fig1_matchStats(matchData, matchPCs_v145, dataFold, paperFold,'v145',savebool);


%% -------------------FIGURE 3: Expert performance-------------------
% Figure 3 stuff starts here:

dataFold2 = fullfile(dataFold,'01_Behav+imaging');
tracebool = true; % Plot lots of traces?

% Make raw data example panels from DK169
fig3_exampleFig(sessionsv4_n(5:6),dataFold2,paperFold,savebool)

% Plot v4 performance
fig3_v4Performance(behavData_v4_n,dataFold,paperFold,savebool);

% Plot v4 heatmaps and traces
iV4n = fig3_v4Heatmaps(matchData, matchFOVs, dataFold, paperFold, savebool);

% Plot sessions split by UCO
fig3_v4UCO(sessionsv4_n,dataFold2,paperFold,savebool);

%% Figure 3 - Training on V4_n and testing on V4_n

% Run GLM pre-processing on trained v4 sessions first
fig3_GLMbatch(sessionsv4_n,dataFold,paperFold,comp,savebool);

% Make GLM data structures for R - start with v4_n
fig3_GLMtaskMod_v4n(sessionsv4_n, pcaVar, dataFold, paperFold, savebool);

% Run GLM fitting and testing on V4_n datasets
for iSesh = 1:numel(sessionsv4_n)
    tSesh = tic;
    glmFold = fullfile(dataFold,'01_Behav+imaging',sprintf('Version_%i',sessionsv4_n(iSesh).version),sessionsv4_n(iSesh).name,...
        sprintf('%s_%s',sessionsv4_n(iSesh).date,sessionsv4_n(iSesh).fov),'GLM','glmRfiles');
    fig3_processGLM(glmFold, glmParams, savebool)
    fprintf('Session took %i seconds.\n\n',round(toc(tSesh)));
end

%% Figure 3 - Make example figure of behaviour and responses for V4_n

% Plot raw traces and heatmaps for each FOV
fig3_GLMtrace_heat(sessionsv4_n,dataFold2,paperFold,glmParams,tracebool,savebool)

% Make summary glm bar graphs and save out significant IDs
fig3_sigSummary(sessionsv4_n,dataFold2,paperFold,glmParams,savebool)

% % Response similarity as a fxn of distance - single params and single FOVs
% % and combined within FOV and recording
% fig2_GLMdistDep(sessionsv4_n,dataFold2,paperFold,glmParams,savebool)

% Analyze movement parameters for each FOV and reward history
plotbool = true;
fig3_behavCorrelates(matchFOVs,glmParams,dataFold,paperFold,plotbool,savebool)

fig3_modDirection(matchData,glmParams,dataFold,paperFold,savebool)

%% -------------------FIGURE 4: slow learning ------------------
% Figure 4 stuff starts here:
% Set common model framework (should be same as in fig 2)

% Combine v144 sessions into mega-response-matrix
matchPCs_v144 = fig2_match_v144(matchData,matchSpont,matchFOVs,dataFold,paperFold,matchColors,savebool);

%% Figure 4 - make GLM design matrices, etc. for v1 and v4_1 sessions - DONE

% Making initial DM structures for v1 sessions
fig3_GLMbatch(sessionsv1_im,dataFold,paperFold,comp,savebool);

% Trimming structures for correct trials only
fig4_GLMtaskMod_v1(sessionsv1_im, pcaVar, dataFold, paperFold, savebool);

% Making initial DM structures
fig3_GLMbatch(sessionsv4_1,dataFold,paperFold,comp,savebool);

% Make GLM data structures for R - v4_1
fig3_GLMtaskMod_v4n(sessionsv4_1, pcaVar, dataFold, paperFold, savebool);

%% Figure 4 - Training on V4_1 and testing on V4_1 - DONE

% Run GLM fitting and testing on V4_1 datasets
for iSesh = 1:numel(sessionsv4_1)
    tSesh = tic;
    glmFold = fullfile(dataFold,'01_Behav+imaging',sprintf('Version_%i',sessionsv4_1(iSesh).version),sessionsv4_1(iSesh).name,...
        sprintf('%s_%s',sessionsv4_1(iSesh).date,sessionsv4_1(iSesh).fov),'GLM','glmRfiles');    
    fig3_processGLM(glmFold,glmParams,savebool)
    fprintf('Session took %i seconds.\n\n',round(toc(tSesh)));
end

dataFold2 = fullfile(dataFold,'01_Behav+imaging');

% Make summary glm bar graphs and save out significant IDs
fig4_sigSummary_v41(sessionsv4_1,dataFold2,paperFold,glmParams,savebool)

%% Figure 4 - Training on correct trials of V4_n and V4_1 to test on V1

% Process GLM correct trials only in V4_n sessions using cvglmnet, etc.
for iSesh = 1:numel(sessionsv4_n)
    tSesh = tic;
    glmFold = fullfile(dataFold,'01_Behav+imaging',sprintf('Version_%i',sessionsv4_n(iSesh).version),sessionsv4_n(iSesh).name,...
        sprintf('%s_%s',sessionsv4_n(iSesh).date,sessionsv4_n(iSesh).fov),'GLM','glmRfiles');
    fig4_processGLM_v4cor(glmFold,glmParams,savebool)
    fprintf('Session took %i seconds.\n\n',round(toc(tSesh)));
end

% Process GLM correct trials only in V4_1 sessions using cvglmnet, etc.
for iSesh = 1:numel(sessionsv4_1)
    tSesh = tic;
    glmFold = fullfile(dataFold,'01_Behav+imaging',sprintf('Version_%i',sessionsv4_1(iSesh).version),sessionsv4_1(iSesh).name,...
        sprintf('%s_%s',sessionsv4_1(iSesh).date,sessionsv4_1(iSesh).fov),'GLM','glmRfiles');    
    fig4_processGLM_v4cor(glmFold,glmParams,savebool)
    fprintf('Session took %i seconds.\n\n',round(toc(tSesh)));
end


% Fit v4n GLM on data from v1 (on correct trials)
[pSig_v1_v4n_lobv, pSig_v1_v4n_sim2] = fig4_matchGLMbars(sessionsv1_im,'v1',sessionsv4_n,'v4n','corrects',glmParams,dataFold,paperFold,savebool); % 'corrects' or 'all' trials to fit
% Fit v41 GLM on data from v1 (on correct trials)
[pSig_v1_v41_lobv, pSig_v1_v41_sim2] = fig4_matchGLMbars(sessionsv1_im,'v1',sessionsv4_1,'v41','corrects',glmParams,dataFold,paperFold,savebool); % 'corrects' or 'all' trials to fit
% Fit v4n GLM on data from v41 (on correct trials)
[pSig_v41_v4n_lobv, pSig_v41_v4n_sim2] = fig4_matchGLMbars(sessionsv4_1,'v41',sessionsv4_n,'v4n','corrects',glmParams,dataFold,paperFold,savebool); % 'corrects' or 'all' trials to fit

lobv_stats = cat(2,pSig_v1_v4n_lobv, pSig_v1_v41_lobv, pSig_v41_v4n_lobv);
sim2_stats = cat(2,pSig_v1_v4n_sim2, pSig_v1_v41_sim2, pSig_v41_v4n_sim2);

[p,tbl,stats] = friedman(lobv_stats);
c = multcompare(stats,'CType', 'bonferroni');

[p,tbl,stats] = friedman(sim2_stats);
c = multcompare(stats,'CType', 'bonferroni');


% Plot task-aligned activity for matched ROIs for example fig
sortFOV = 2; % Sort trial by v4_n session (which is second)
fig4_matchFOVs(sessionsv1_im,'v1',sessionsv4_n,'v4n',sortFOV,dataFold,paperFold,savebool)
fig4_matchFOVs(sessionsv1_im,'v1',sessionsv4_1,'v41',sortFOV,dataFold,paperFold,savebool)
fig4_matchFOVs(sessionsv4_1,'v41',sessionsv4_n,'v4n',sortFOV,dataFold,paperFold,savebool)

%% Figure 4 - make percent significant heatmaps

trainSesh = {'v41','v4n'}; % Use both v41 and v4n as prediction sessions
testSesh = {'v1','v41'}; % Use v1 and v41 as testing
fig4_matchPercSig(trainSesh,testSesh,'corrects',paperFold,savebool);


%% Figure 4 - fit movement and reward kernels using v4n data

fig4_fitKernels(sessionsv4_n,sessionsv1_im,sessionsv4_1,glmParams,dataFold,paperFold,savebool);

%% Figure 4 - kernel similarity as a fxn of distance

% Extract relative ROI positions
sessionsv4_n = fig4_mlPosition(sessionsv4_n, matchData);
sessionsv5_im = fig4_mlPosition(sessionsv5_im, matchData);
sessionsv4_1 = fig4_mlPosition(sessionsv4_1, matchData);
sessionsv1_im = fig4_mlPosition(sessionsv1_im, matchData);

% 
fig4_kernelRvsDist(sessionsv4_n,paperFold,savebool)
fig4_kernelRvsDist_v4only(sessionsv4_n,paperFold,savebool)


%% -------------------FIGURE 5: Fast learning-------------------
% Figure 5 stuff starts here:
% Set common model framework (should be same as in fig 3)

% % GLM parameters: this will determine which v4_n GLM to fish out as well as
% % what parameters to use to fit models in this figure
% glmParams.dataType      = 'spk'; % 'fluo' or 'spk' for zF or spike events
% glmParams.distr         = 'gaussian'; % 'gaussian' (for fluo data), 'poisson' (for spk counts), or 'binomial' (for binarized spk trains)
% glmParams.link          = []; % Use default link function for given distribution
% glmParams.alpha         = .5; % Others use 0.5 for gaussian GLMs and 0.8-0.9 for poisson or binomial GLMs
% glmParams.nLambda       = 100; % Good number, don't change from 100

% Run GLM pre-processing on trained v5 sessions first
fig3_GLMbatch(sessionsv5_im,dataFold,paperFold,comp,savebool);

% Make GLM data structures for R - start with v4_5
fig3_GLMtaskMod_v4n(sessionsv5_im, pcaVar, dataFold, paperFold, savebool);


%% Make example figures showing washin and washout using mov and reward cells from v4n and pool into summary figure
fig5_adaptAnalysis(matchData, glmParams, dataFold, paperFold, savebool)

% Calculate distance-dependence
fig5_adaptRatiovsDist(sessionsv4_n,paperFold,savebool)

%% -------------------FIGURE 6: Comparing fast and slow learning-------------------
% Figure 6 stuff starts here:

% Compare v1-v4 to changes to acute v4 changes (mov correlation and reward
% history)
fig6_part1

% Comparing reliabliity:
fig6_part2(matchData, glmParams, dataFold, paperFold, savebool)

% Comparing v1-4 chagnes to v5 adaptation:
fig6_analysis_sfn(matchData, paperFold, savebool)

