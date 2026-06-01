function fig4_processGLM_v4cor(GLMroot, glmParams, savebool)
%% GLM script for single dataset:
% set root dir
% GLMroot = '/Users/dimitar/Desktop/CF_learning_paper/01_Behav+imaging/Version_4/DK103/2018-02-22_lobv/GLM/glmRfiles';

% % GLM parameters:
% glmParams.dataType      = 'spk'; % 'fluo' or 'spk' for zF or spike events
% glmParams.distr         = 'gaussian'; % 'gaussian' (for fluo data), 'poisson' (for spk counts), or 'binomial' (for binarized spk trains)
% glmParams.link          = []; % Use default link function for given distribution
% glmParams.alpha         = .5; % Others use 0.5 for gaussian GLMs and 0.8-0.9 for poisson or binomial GLMs
% glmParams.nLambda       = 100; % Good number, don't change from 100


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
% nPCvar = DMpcadata.nPCvar;
glmParams.maxPred = find(cumsum(DMpcadata.exp) > 99,1);
fPathSplit = strsplit(GLMroot,filesep);
cd(GLMroot);

switch glmParams.dataType
    case 'fluo'
        Ydata = double(zFtrials);
    case 'spk'
        Ydata = double(Ytrials);
        if strcmpi(glmParams.distr,'binomial')
            Ydata = double(logical(Ydata));
        end
end

glmFold = fullfile(GLMroot,sprintf('GLMfull_%s_%s_alpha%i',glmParams.dataType,glmParams.distr,glmParams.alpha*100));
if ~exist(glmFold,'dir'); mkdir(glmFold); end

%% Run full GLM
glmFile = fullfile(glmFold,'GLMsigCor.mat');
if ~exist(glmFile,'file')
    fprintf('Mouse name:       %s\n',fPathSplit{5});
    fprintf('Session date+FOV: %s\n',fPathSplit{6});
    fprintf('Session type:     Full GLM\n');    
    % Empirically determine PCA rank that gives best fit for each cell and run GLM + shuffles:        
    GLMoutput = fig2_makeLassoGLM(DMpcadata.score,DM,Ydata,glmParams,iFoldsmall,iEndTrial,nShuf);
    if savebool
        save(glmFile,'GLMoutput');
    end
else
    load(glmFile);
end

%% Remove predictors and test
% Visual predictors
glmFile = fullfile(glmFold,'GLMsigCor+Vis.mat');
if ~exist(glmFile,'file')
    fprintf('Mouse name:       %s\n',fPathSplit{5});
    fprintf('Session date+FOV: %s\n',fPathSplit{6});
    fprintf('Session type:     Visual predictors\n');
    permPreds = {'TrialON','BSwheelpos'};    
    GLMoutput_visual = fig2_permutePreds(GLMoutput,DM,DMdata,Ydata,glmParams,iFoldsmall,iEndTrial,permPreds,nShuf);
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
glmFile = fullfile(glmFold,'GLMsigCor+VisMov.mat');
if ~exist(glmFile,'file')
    fprintf('Mouse name:       %s\n',fPathSplit{5});
    fprintf('Session date+FOV: %s\n',fPathSplit{6});
    fprintf('Session type:     Movement predictors\n');
    permPreds = {'MovT','BSwheelvelL','BSwheelvelR'};
    GLMoutput_mov = fig2_permutePreds(GLMoutput,DM,DMdata,Ydata,glmParams,iFoldsmall,iEndTrial,permPreds,nShuf);    
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
glmFile = fullfile(glmFold,'GLMsigCor+VisMovRew.mat');
if ~exist(glmFile,'file')
    fprintf('Mouse name:       %s\n',fPathSplit{5});
    fprintf('Session date+FOV: %s\n',fPathSplit{6});
    fprintf('Session type:     Reward predictors\n');
    permPreds = {'RewTcorrect'};
    GLMoutput_rew = fig2_permutePreds(GLMoutput,DM,DMdata,Ydata,glmParams,iFoldsmall,iEndTrial,permPreds,nShuf);
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
glmFile = fullfile(glmFold,'GLMsigCor+VisMovRewLick.mat');
if ~exist(glmFile,'file')
    fprintf('Mouse name:       %s\n',fPathSplit{5});
    fprintf('Session date+FOV: %s\n',fPathSplit{6});
    fprintf('Session type:     Licking predictors\n');
    permPreds = {'LickT'};
    GLMoutput_lick = fig2_permutePreds(GLMoutput,DM,DMdata,Ydata,glmParams,iFoldsmall,iEndTrial,permPreds,nShuf);
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

