function fig3_processGLM_v4cor_rerun(GLMroot, savebool)
%% GLM script for single dataset:
% set root dir
% GLMroot = '/Users/dimitar/Desktop/CF_learning_paper/01_Behav+imaging/Version_4/DK103/2018-02-22_lobv/GLM/glmRfiles';
% savebool = false;

%% Load in data
GLMdata = load(fullfile(GLMroot,'GLMcor_small.mat'));
v2struct(GLMdata.GLMcor_small);
DMpcadata = DMpca; %#ok<*NODEF>
DMpca = DMpcadata.score;
% DMpca = DMpcadata.score(:,1:DMpca.nPCvar);
% iUCOsmall = readmatrix(fullfile(GLMroot,'CorData_small','iUCOSmall.csv'));
iFoldsmall = readmatrix(fullfile(GLMroot,'CorData_small','iFoldCSmall.csv'));

%% Set up inits
iEndTrial = DMdata.endTrialIndices;
nShuf = 20;
nPCvar = DMpcadata.nPCvar;
maxPred = find(cumsum(DMpcadata.exp) > 99,1);
fPathSplit = strsplit(GLMroot,filesep);
cd(GLMroot);

%% Run full GLM
glmFile = fullfile(GLMroot,'GLMsigCor.mat');
load(glmFile);
if ~isfield(GLMoutput,'GLMobj_shuf')
    fprintf('Mouse name:       %s\n',fPathSplit{5});
    fprintf('Session date+FOV: %s\n',fPathSplit{6});
    fprintf('Session type:     Full GLM\n');  
    permPreds = DMdata.fields;
    [GLMoutput.('devExp_shuf_orig')] = GLMoutput.('devExp_shuf');
    GLMoutput = fig2_permutePreds_rerun(GLMoutput,DM,DMdata,Ytrials,iFoldsmall,iEndTrial,nPCvar,permPreds,nShuf,0);    
    
    if savebool
        save(glmFile,'GLMoutput');
    end
end

%% Remove predictors and test
% Visual predictors
glmFile = fullfile(GLMroot,'GLMsigCor+Vis.mat');
glmOuttemp = load(glmFile);
if ~isfield(glmOuttemp.GLMoutput,'GLMobj_shufVis')
    fprintf('Mouse name:       %s\n',fPathSplit{5});
    fprintf('Session date+FOV: %s\n',fPathSplit{6});
    fprintf('Session type:     Visual predictors\n');
    permPreds = {'TrialON','BSwheelpos'};
    GLMoutput_visual = fig2_permutePreds_rerun(GLMoutput,DM,DMdata,Ytrials,iFoldsmall,iEndTrial,nPCvar,permPreds,nShuf,0);
    for iCell = 1:numel(GLMoutput_visual)
        GLMoutput(iCell).GLMobj_shufVis = GLMoutput_visual(iCell).GLMobj_shuf;
        GLMoutput(iCell).yHat_shufVis_mu = GLMoutput_visual(iCell).yHat_shuf_mu;
        GLMoutput(iCell).yHat_shufVis_sd = GLMoutput_visual(iCell).yHat_shuf_sd;
        GLMoutput(iCell).devExp_shufVis = GLMoutput_visual(iCell).devExp_shuf;
        GLMoutput(iCell).sigBool_shufVis = GLMoutput_visual(iCell).sigBool_shuf;
        GLMoutput(iCell).GLMobj_shufVis_proj = GLMoutput_visual(iCell).GLMobj_shuf_proj;
        GLMoutput(iCell).yHat_shufVis_proj_mu = GLMoutput_visual(iCell).yHat_shuf_proj_mu;
        GLMoutput(iCell).yHat_shufVis_proj_sd = GLMoutput_visual(iCell).yHat_shuf_proj_sd;
        GLMoutput(iCell).devExp_shufVis_proj = GLMoutput_visual(iCell).devExp_shuf_proj;
        GLMoutput(iCell).sigBool_shufVis_proj = GLMoutput_visual(iCell).sigBool_proj;
    end
    if savebool
        save(glmFile,'GLMoutput');
    end
else
    load(glmFile);
end

% Movement predictors
glmFile = fullfile(GLMroot,'GLMsigCor+VisMov.mat');
glmOuttemp = load(glmFile);
if ~isfield(glmOuttemp.GLMoutput,'GLMobj_shufMov')
    fprintf('Mouse name:       %s\n',fPathSplit{5});
    fprintf('Session date+FOV: %s\n',fPathSplit{6});
    fprintf('Session type:     Movement predictors\n');
    permPreds = {'MovT','BSwheelvelL','BSwheelvelR'};
    GLMoutput_mov = fig2_permutePreds_rerun(GLMoutput,DM,DMdata,Ytrials,iFoldsmall,iEndTrial,nPCvar,permPreds,nShuf,0);
    for iCell = 1:numel(GLMoutput_mov)
        GLMoutput(iCell).GLMobj_shufMov = GLMoutput_mov(iCell).GLMobj_shuf;
        GLMoutput(iCell).yHat_shufMov_mu = GLMoutput_mov(iCell).yHat_shuf_mu;
        GLMoutput(iCell).yHat_shufMov_sd = GLMoutput_mov(iCell).yHat_shuf_sd;
        GLMoutput(iCell).devExp_shufMov = GLMoutput_mov(iCell).devExp_shuf;
        GLMoutput(iCell).sigBool_shufMov = GLMoutput_mov(iCell).sigBool;
        GLMoutput(iCell).GLMobj_shufMov_proj = GLMoutput_mov(iCell).GLMobj_shuf_proj;
        GLMoutput(iCell).yHat_shufMov_proj_mu = GLMoutput_mov(iCell).yHat_shuf_proj_mu;
        GLMoutput(iCell).yHat_shufMov_proj_sd = GLMoutput_mov(iCell).yHat_shuf_proj_sd;
        GLMoutput(iCell).devExp_shufMov_proj = GLMoutput_mov(iCell).devExp_shuf_proj;
        GLMoutput(iCell).sigBool_shufMov_proj = GLMoutput_mov(iCell).sigBool_proj;
    end
    if savebool
        save(glmFile,'GLMoutput');
    end
else
    load(glmFile);
end

% Reward predictors
glmFile = fullfile(GLMroot,'GLMsigCor+VisMovRew.mat');
glmOuttemp = load(glmFile);
if ~isfield(glmOuttemp.GLMoutput,'GLMobj_shufRew')
    fprintf('Mouse name:       %s\n',fPathSplit{5});
    fprintf('Session date+FOV: %s\n',fPathSplit{6});
    fprintf('Session type:     Reward predictors\n');
    permPreds = {'RewTcorrect'};
    GLMoutput_rew = fig2_permutePreds(GLMoutput,DM,DMdata,Ytrials,iFoldsmall,iEndTrial,nPCvar,permPreds,nShuf,0);    
    for iCell = 1:numel(GLMoutput_rew)
        GLMoutput(iCell).GLMobj_shufRew = GLMoutput_rew(iCell).GLMobj_shuf;
        GLMoutput(iCell).yHat_shufRew_mu = GLMoutput_rew(iCell).yHat_shuf_mu;
        GLMoutput(iCell).yHat_shufRew_sd = GLMoutput_rew(iCell).yHat_shuf_sd;
        GLMoutput(iCell).devExp_shufRew = GLMoutput_rew(iCell).devExp_shuf;
        GLMoutput(iCell).sigBool_shufRew = GLMoutput_rew(iCell).sigBool;
        GLMoutput(iCell).GLMobj_shufRew_proj = GLMoutput_rew(iCell).GLMobj_shuf_proj;
        GLMoutput(iCell).yHat_shufRew_proj_mu = GLMoutput_rew(iCell).yHat_shuf_proj_mu;
        GLMoutput(iCell).yHat_shufRew_proj_sd = GLMoutput_rew(iCell).yHat_shuf_proj_sd;
        GLMoutput(iCell).devExp_shufRew_proj = GLMoutput_rew(iCell).devExp_shuf_proj;
        GLMoutput(iCell).sigBool_shufRew_proj = GLMoutput_rew(iCell).sigBool_proj;
    end
    if savebool
        save(glmFile,'GLMoutput');
    end
else
    load(glmFile);
end

% Lick predictors
glmFile = fullfile(GLMroot,'GLMsigCor+VisMovRewLick.mat');
glmOuttemp = load(glmFile);
if ~isfield(glmOuttemp.GLMoutput,'GLMobj_shufLick')
    fprintf('Mouse name:       %s\n',fPathSplit{5});
    fprintf('Session date+FOV: %s\n',fPathSplit{6});
    fprintf('Session type:     Licking predictors\n');
    permPreds = {'LickT'};
    GLMoutput_lick = fig2_permutePreds(GLMoutput,DM,DMdata,Ytrials,iFoldsmall,iEndTrial,nPCvar,permPreds,nShuf,0);    
    for iCell = 1:numel(GLMoutput_lick)
        GLMoutput(iCell).GLMobj_shufLick = GLMoutput_lick(iCell).GLMobj_shuf;
        GLMoutput(iCell).yHat_shufLick_mu = GLMoutput_lick(iCell).yHat_shuf_mu;
        GLMoutput(iCell).yHat_shufLick_sd = GLMoutput_lick(iCell).yHat_shuf_sd;
        GLMoutput(iCell).devExp_shufLick = GLMoutput_lick(iCell).devExp_shuf;
        GLMoutput(iCell).sigBool_shufLick = GLMoutput_lick(iCell).sigBool;
        GLMoutput(iCell).GLMobj_shufLick_proj = GLMoutput_lick(iCell).GLMobj_shuf_proj;
        GLMoutput(iCell).yHat_shufLick_proj_mu = GLMoutput_lick(iCell).yHat_shuf_proj_mu;
        GLMoutput(iCell).yHat_shufLick_proj_sd = GLMoutput_lick(iCell).yHat_shuf_proj_sd;
        GLMoutput(iCell).devExp_shufLick_proj = GLMoutput_lick(iCell).devExp_shuf_proj;
        GLMoutput(iCell).sigBool_shufLick_proj = GLMoutput_lick(iCell).sigBool_proj;
    end
    if savebool
        save(glmFile,'GLMoutput');
    end
% else
%     load(glmFile);
end

