function fig2_processGLMcor(GLMroot, savebool)
%% GLM script for single dataset:
% set root dir
% GLMroot = '/Users/dimitar/Desktop/CF_learning_paper/01_Behav+imaging/Version_4/DK103/2018-02-22_lobv/GLM/glmRfiles';
% savebool = false;

%% Load in data
GLMdata = load(fullfile(GLMroot,'GLMcor_small.mat'));
v2struct(GLMdata.GLMcor_small);
DMpcadata = DMpca;
DMpca = DMpca.score(:,1:DMpca.nPCvar);
% iUCOsmall = readmatrix(fullfile(GLMroot,'FullData_small','iUCOSmall.csv'));
iFoldsmall = readmatrix(fullfile(GLMroot,'CorData_small','iFoldCsmall.csv'));

%% Set up inits
iEndTrial = DMdata.endTrialIndices;
nShuf = 20;
nPCvar = DMpcadata.nPCvar;
fPathSplit = strsplit(GLMroot,filesep);
cd(GLMroot);

%% Run full GLM
glmFile = fullfile(GLMroot,'GLMsigC.mat');
if ~exist(glmFile,'file')
    fprintf('Mouse name:       %s\n',fPathSplit{5});
    fprintf('Session date+FOV: %s\n',fPathSplit{6});
    fprintf('Session type:     Correct trial GLM\n');
    GLMoutput = fig2_makeLassoGLM(DMpca,Ytrials,iFoldsmall,iEndTrial,nShuf,0);
    if savebool
        save(glmFile,'GLMoutput');
    end
end
%% Remove predictors and test
% Visual predictors
glmFile = fullfile(GLMroot,'GLMsigC+Vis.mat');
if ~exist(glmFile,'file')
    fprintf('Mouse name:       %s\n',fPathSplit{5});
    fprintf('Session date+FOV: %s\n',fPathSplit{6});
    fprintf('Session type:     Visual predictors\n');
    permPreds = {'TrialON','BSwheelpos'};
    GLMoutput_visual = fig2_permutePreds(DM,DMdata,Ytrials,iFoldsmall,iEndTrial,nPCvar,permPreds,nShuf);
    for iCell = 1:numel(GLMoutput_visual)
        GLMoutput(iCell).yHat_shufVis_mu = GLMoutput_visual(iCell).yHat_shuf_mu;
        GLMoutput(iCell).yHat_shufVis_sd = GLMoutput_visual(iCell).yHat_shuf_sd;
        GLMoutput(iCell).devExp_shufVis = GLMoutput_visual(iCell).devExp_shuf;
        GLMoutput(iCell).sigBool_shufVis = GLMoutput_visual(iCell).sigBool;
    end
    if savebool
        save(glmFile,'GLMoutput');
    end
end

% Movement predictors
glmFile = fullfile(GLMroot,'GLMsigC+VisMov.mat');
if ~exist(glmFile,'file')
    fprintf('Mouse name:       %s\n',fPathSplit{5});
    fprintf('Session date+FOV: %s\n',fPathSplit{6});
    fprintf('Session type:     Movement predictors\n');
    permPreds = {'MovT','BSwheelvelL','BSwheelvelR'};
    GLMoutput_mov = fig2_permutePreds(DM,DMdata,Ytrials,iFoldsmall,iEndTrial,nPCvar,permPreds,nShuf);
    for iCell = 1:numel(GLMoutput_mov)
        GLMoutput(iCell).yHat_shufMov_mu = GLMoutput_mov(iCell).yHat_shuf_mu;
        GLMoutput(iCell).yHat_shufMov_sd = GLMoutput_mov(iCell).yHat_shuf_sd;
        GLMoutput(iCell).devExp_shufMov = GLMoutput_mov(iCell).devExp_shuf;
        GLMoutput(iCell).sigBool_shufMov = GLMoutput_mov(iCell).sigBool;
    end
    if savebool
        save(glmFile,'GLMoutput');
    end
end

% Reward predictors
glmFile = fullfile(GLMroot,'GLMsigC+VisMovRew.mat');
if ~exist(glmFile,'file')
    fprintf('Mouse name:       %s\n',fPathSplit{5});
    fprintf('Session date+FOV: %s\n',fPathSplit{6});
    fprintf('Session type:     Reward predictors\n');
    permPreds = {'RewTcorrect'};
    GLMoutput_rew = fig2_permutePreds(DM,DMdata,Ytrials,iFoldsmall,iEndTrial,nPCvar,permPreds,nShuf);
    for iCell = 1:numel(GLMoutput_rew)
        GLMoutput(iCell).yHat_shufRew_mu = GLMoutput_rew(iCell).yHat_shuf_mu;
        GLMoutput(iCell).yHat_shufRew_sd = GLMoutput_rew(iCell).yHat_shuf_sd;
        GLMoutput(iCell).devExp_shufRew = GLMoutput_rew(iCell).devExp_shuf;
        GLMoutput(iCell).sigBool_shufRew = GLMoutput_rew(iCell).sigBool;
    end
    if savebool
        save(glmFile,'GLMoutput');
    end
end

% Lick predictors
glmFile = fullfile(GLMroot,'GLMsigC+VisMovRewLick.mat');
if ~exist(glmFile,'file')
    fprintf('Mouse name:       %s\n',fPathSplit{5});
    fprintf('Session date+FOV: %s\n',fPathSplit{6});
    fprintf('Session type:     Licking predictors\n');
    permPreds = {'LickT'};
    GLMoutput_lick = fig2_permutePreds(DM,DMdata,Ytrials,iFoldsmall,iEndTrial,nPCvar,permPreds,nShuf);
    for iCell = 1:numel(GLMoutput_lick)
        GLMoutput(iCell).yHat_shufLick_mu = GLMoutput_lick(iCell).yHat_shuf_mu;
        GLMoutput(iCell).yHat_shufLick_sd = GLMoutput_lick(iCell).yHat_shuf_sd;
        GLMoutput(iCell).devExp_shufLick = GLMoutput_lick(iCell).devExp_shuf;
        GLMoutput(iCell).sigBool_shufLick = GLMoutput_lick(iCell).sigBool;
    end
    if savebool
        save(glmFile,'GLMoutput');
    end
end

