function GLMoutput = fig4_fitV4corr_kernels(GLMroot, glmParams)
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
iFoldsmall = readmatrix(fullfile(GLMroot,'CorData_small','iFoldCSmall.csv'));
glmFold = fullfile(GLMroot,sprintf('GLMfull_%s_%s_alpha%i',glmParams.dataType,glmParams.distr,glmParams.alpha*100));
if ~exist(glmFold,'dir'); mkdir(glmFold); end
GLMinput = load(fullfile(glmFold,'GLMsigCor+VisMovRewLick.mat')); GLMinput = GLMinput.GLMoutput;

%% Set up inits
glmFields = fieldnames(GLMinput);
keepFields = [1,3,12];
irmFields = true(numel(glmFields),1); irmFields(keepFields) = false;
GLMoutput = rmfield(GLMinput,glmFields(irmFields));

switch glmParams.dataType
    case 'fluo'
        Ydata = double(zFtrials);
    case 'spk'
        Ydata = double(Ytrials);
        if strcmpi(glmParams.distr,'binomial')
            Ydata = double(logical(Ydata));
        end
end

%% Get GLM fits on v4 for movement and reward
fPathSplit = strsplit(GLMroot,filesep);
fprintf('Mouse name:       %s\n',fPathSplit{5});
fprintf('Session date+FOV: %s\n',fPathSplit{6});
fprintf('Session type:     Expert movement kernel fitting\n');
permGrp = 'Mov'; % 'Vis','Mov','Rew',or 'Lick' (case-sensitive)
GLMmov = fig4_fitV4corr_fitData(GLMinput,DM,DMdata,Ydata,glmParams,permGrp);
fprintf('Session type:     Expert reward kernel fitting\n');
permGrp = 'Rew'; % 'Vis','Mov','Rew',or 'Lick' (case-sensitive)
GLMrew = fig4_fitV4corr_fitData(GLMinput,DM,DMdata,Ydata,glmParams,permGrp);
for iCell = 1:numel(GLMoutput)
    GLMoutput(iCell).yData = GLMmov(iCell).yData;
    GLMoutput(iCell).sigBoolMov = GLMmov(iCell).sigBoolLocal;
    GLMoutput(iCell).GLMobjMov = GLMmov(iCell).GLMobjLocal;
    GLMoutput(iCell).kernelMov = GLMmov(iCell).kernelLocal;
    GLMoutput(iCell).sigBoolRew = GLMrew(iCell).sigBoolLocal;
    GLMoutput(iCell).GLMobjRew = GLMrew(iCell).GLMobjLocal;
    GLMoutput(iCell).kernelRew = GLMrew(iCell).kernelLocal;
end

end
