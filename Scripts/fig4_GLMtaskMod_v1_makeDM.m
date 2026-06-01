function [GLMfull, GLMfull_small, GLMcor, GLMcor_small] = fig4_GLMtaskMod_v1_makeDM(trialBasis, inits,trialOutcomes,glmFold,savebool)
%% Extract things from inits that we care about
xvalFold = inits.xvalfold;

%% Take rewarded trials
nTrials0 = numel(trialBasis.expt);
trialsC = false(nTrials0,1);
trialsL = false(nTrials0,1);
for iTrial = 1:nTrials0
    if trialBasis.expt(iTrial).objpos(1) < 1
        trialsL(iTrial) = true;
    end
    if ~isempty(trialBasis.expt(iTrial).reward)
        trialsC(iTrial) = true;
    end
end

trialBasis.expt = trialBasis.expt(trialsC & trialsL);
trialBasis.bsbox = trialBasis.bsbox(trialsC & trialsL);
nTrials = numel(trialBasis.expt);
trialsLCall = 1:nTrials;

starti = 2;

trialsLC_holdout = trialsLCall(starti:3:end);

clear trialsALL0 trialsLU0 trialsLC0 trialsLO0 iTrials0 iRews0 trialout

%% Gather info to feed to glmnet - R
% Remove reward fields from design matrix
trialBasis.bsboxC = rmfield(trialBasis.bsbox,{'TrialOFFunder','TrialOFFover','RewTunder','RewTover','CueBeep','CueReward','RandReward'});
trialBasis = rmfield(trialBasis,'bsbox');

% Make correct trials only design matrix
[DM_C,DM_Cdata] = makeDesignMatrix_dk(trialBasis, trialsLCall, 'bsboxC');

% Trim off predictors that are during spont activity/witholds
[DM_CdataSmall,iPredCSmall, iTrialCSmall] = reducePredictorsBox(DM_Cdata,inits);
DM_CSmall = DM_C(iTrialCSmall,iPredCSmall);

% % Make matrix of spiking data
Ytrial_C = cat(1,trialBasis.expt(trialsLCall).spkData);

Ytrial_CSmall = Ytrial_C(iTrialCSmall,:);

%% use PCA to reduce design matrix
pcaVar = inits.pcaVar;

% % Correct trials only: long format - all predictors and delay periods
% DM_Cpca = struct;
% [DM_Cpca.coef, DM_Cpca.score, DM_Cpca.lat, ~, DM_Cpca.exp] = pca(DM_C);
% DM_Cpca.pcaVar = pcaVar;
% DM_Cpca.nPCvar = sum(cumsum(DM_Cpca.exp) < pcaVar)+1;
% DM = DM_C; DMdata = DM_Cdata; DMpca = DM_Cpca; Ytrials = Ytrial_C;
% GLMcor = v2struct(DM, DMdata, Ytrials, DMpca);
% clear DM DMdata DMpca Ytrials

% Correct trials only: short format - trimmed trials and predictors
DM_CSmallpca = struct;
[DM_CSmallpca.coef, DM_CSmallpca.score, DM_CSmallpca.lat, ~, DM_CSmallpca.exp] = pca(DM_CSmall);
DM_CSmallpca.pcaVar = pcaVar;
DM_CSmallpca.nPCvar = sum(cumsum(DM_CSmallpca.exp) < pcaVar)+1;
DM = DM_CSmall; DMdata = DM_CdataSmall; DMpca = DM_CSmallpca; Ytrials = Ytrial_CSmall;
GLMcor_small = v2struct(DM, DMdata, Ytrials, DMpca);
clear DM DMdata DMpca Ytrials

if savebool
    RglmFold = fullfile(glmFold,'glmRfiles');
    if ~exist(RglmFold,'dir'); mkdir(RglmFold); end
%     save(fullfile(RglmFold,'GLMfull.mat'),'GLMfull');
%     save(fullfile(RglmFold,'GLMfull_small.mat'),'GLMfull_small');
%     save(fullfile(RglmFold,'GLMcor.mat'),'GLMcor');
    save(fullfile(RglmFold,'GLMcor_small.mat'),'GLMcor_small');
end

%% Plot DMs to make sure it's all working
% tmax = 1000;
% nBases = DM_CdataSmall.nBases;
% nBasesEnd = cumsum(DM_CdataSmall.nBase);
% 
% figure;
% subplot(3,4,1); imagesc(DM_CSmall(1:tmax,:)',[-1 1]); colormap redblue
% hold on; for i = 1:20; xline(DM_CdataSmall.endTrialIndices(i),'--k','LineWidth',1.5); end
% hold on; for i = 1:numel(nBasesEnd); yline(nBasesEnd(i),'--k','LineWidth',1.5); end
% xlim([0 tmax]);
% subplot(3,4,1); imagesc(DM_CSmall(1:tmax,:)',[-1 1]); colormap redblue
% hold on; for i = 1:20; xline(DM_CdataSmall.endTrialIndices(i),'--k','LineWidth',1.5); end
% hold on; for i = 1:numel(nBasesEnd); yline(nBasesEnd(i),'--k','LineWidth',1.5); end
% xlim([0 tmax]);
% subplot(3,1,2);
% imagesc(DMboxCscore(1:tmax,:)',[-1 1]); colormap redblue
% hold on; for i = 1:20; xline(DM_CdataSmall.endTrialIndices(i),'--k','LineWidth',1.5); end
% subplot(3,1,3);
% imagesc(DMboxCscore(1:tmax,1:nPC80)',[-1 1]); colormap redblue
% hold on; for i = 1:20; xline(DM_CdataSmall.endTrialIndices(i),'--k','LineWidth',1.5); end


% %% Save maybe
% if savebool
%     % Correct trials only - trimmed trials and predictors
%     RglmFold_4 = fullfile(RglmFold,'CorData_small');
%     if ~exist(RglmFold_4,'dir'); mkdir(RglmFold_4); end
%     writematrix(iTrialsCsmall,fullfile(RglmFold_4,'iTrialsCsmall.csv'));
%     writematrix(iTrialsGlobalCsmall,fullfile(RglmFold_4,'iTrialsGlobalCsmall.csv'));
%     writematrix(DM_CSmall,fullfile(RglmFold_4,'DMcorSmall.csv'));
%     writematrix(DM_CSmallpca.score,fullfile(RglmFold_4,'DMcorSmallPCAscore.csv'));
%     writematrix(DM_CSmallpca.nPCvar,fullfile(RglmFold_4,'nPCvar.csv'));
%     writematrix(Ytrial_CSmall,fullfile(RglmFold_4,'YcorSmall.csv'));
%     writematrix(iFoldCsmall,fullfile(RglmFold_4,'iFoldCsmall.csv'));
%     writematrix(iFoldCsmall(iFoldCsmall <= xvalFold),fullfile(RglmFold_4,'iFoldCSmallTest.csv'));
%     writematrix(iFoldCsmall(iFoldCsmall > xvalFold),fullfile(RglmFold_4,'iFoldCSmallHoldout.csv'));
%     
% end

end

