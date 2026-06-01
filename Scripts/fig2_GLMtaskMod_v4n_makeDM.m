function [GLMfull, GLMfull_small, GLMcor, GLMcor_small] = fig2_GLMtaskMod_v4n_makeDM(trialBasis, inits,trialOutcomes,glmFold,savebool)
%% Extract things from inits that we care about
xvalFold = inits.xvalfold;

%% Split trials into uco
trialsLU0 = trialOutcomes.trialsLU0;
trialsLC0 = trialOutcomes.trialsLC0;
trialsLO0 = trialOutcomes.trialsLO0;

trialsALL0 = cat(1,trialBasis.expt.trialnum);
trialsGOOD = ismember(trialsALL0,[trialsLU0; trialsLC0; trialsLO0]);

trialBasis.expt = trialBasis.expt(trialsGOOD);
trialBasis.bsbox = trialBasis.bsbox(trialsGOOD);

trialsALL = cat(1,trialBasis.expt.trialnum);
trialsLUall = find(ismember(trialsALL,trialsLU0));
trialsLCall = find(ismember(trialsALL,trialsLC0));
trialsLOall = find(ismember(trialsALL,trialsLO0));
nTrials = numel(trialsALL);

starti = 2;

trialsLU_holdout = trialsLUall(starti:3:end);
trialsLC_holdout = trialsLCall(starti:3:end);
trialsLO_holdout = trialsLOall(starti:3:end);

clear trialsALL0 trialsLU0 trialsLC0 trialsLO0 iTrials0 iRews0 trialout

%% Gather info to feed to glmnet - R
% Remove reward fields from design matrix
trialBasis.bsboxUCO = rmfield(trialBasis.bsbox,{'CueBeep','CueReward','RandReward'});
trialBasis.bsboxC = rmfield(trialBasis.bsboxUCO,{'TrialOFFunder','TrialOFFover','RewTunder','RewTover'});
trialBasis = rmfield(trialBasis,'bsbox');

% Make full design matrix
[DM_UCO,DM_UCOdata] = makeDesignMatrix_dk(trialBasis, 1:nTrials, 'bsboxUCO');
% Make correct trials only design matrix
[DM_C,DM_Cdata] = makeDesignMatrix_dk(trialBasis, trialsLCall, 'bsboxC');

% Trim off predictors that are during spont activity/witholds
[DM_UCOdataSmall,iPredUCOSmall, iTrialUCOSmall] = reducePredictorsBox(DM_UCOdata,inits);
DM_UCOSmall = DM_UCO(iTrialUCOSmall,iPredUCOSmall);
[DM_CdataSmall,iPredCSmall, iTrialCSmall] = reducePredictorsBox(DM_Cdata,inits);
DM_CSmall = DM_C(iTrialCSmall,iPredCSmall);

% % Make matrix of spiking data
% Ytrial_UCO = zeros(size(DM_UCO,1),numel(spkNames));
% Ytrial_C = zeros(size(DM_C,1),numel(spkNames));
% for j = 1:numel(spkNames)
%     Ytrial_UCO(:,j) = getSpikes_dk(trialBasis,1:nTrials,spkNames{j});
%     Ytrial_C(:,j) = getSpikes_dk(trialBasis,trialsLCall,spkNames{j});
% end
Ytrial_UCO = cat(1,trialBasis.expt(:).spkData);
Ytrial_C = cat(1,trialBasis.expt(trialsLCall).spkData);
zFtrial_UCO = cat(1,trialBasis.expt(:).zFData);
zFtrial_C = cat(1,trialBasis.expt(trialsLCall).zFData);

Ytrial_UCOSmall = Ytrial_UCO(iTrialUCOSmall,:);
Ytrial_CSmall = Ytrial_C(iTrialCSmall,:);
zFtrial_UCOSmall = zFtrial_UCO(iTrialUCOSmall,:);
zFtrial_CSmall = zFtrial_C(iTrialCSmall,:);


% Make vectors that denote trial outcome, x-val fold, holdout status and
% local and global trial number
iUCO = zeros(size(DM_UCO,1),3);
iTrials = zeros(size(DM_UCO,1),1); iFold = iTrials; iTrialsGlobal = iTrials;
iU = 0; iC = 0; iO = 0;
for i = 1:nTrials
    iTrials(DM_UCOdata.endTrialIndices(i)+1:DM_UCOdata.endTrialIndices(i+1),1) = i;
    iTrialsGlobal(DM_UCOdata.endTrialIndices(i)+1:DM_UCOdata.endTrialIndices(i+1),1) = trialsALL(i);
    if ismember(i,trialsLUall)
        iUCO(DM_UCOdata.endTrialIndices(i)+1:DM_UCOdata.endTrialIndices(i+1),1) = 1;
        if ismember(i,trialsLU_holdout)
            iFold(DM_UCOdata.endTrialIndices(i)+1:DM_UCOdata.endTrialIndices(i+1),1) = xvalFold+1;
        else
            iU = iU+1;
            iFoldLocal = mod(iU, xvalFold); if iFoldLocal == 0; iFoldLocal = xvalFold; end
            iFold(DM_UCOdata.endTrialIndices(i)+1:DM_UCOdata.endTrialIndices(i+1),1) = iFoldLocal;
        end
    end
    if ismember(i,trialsLCall)
        iUCO(DM_UCOdata.endTrialIndices(i)+1:DM_UCOdata.endTrialIndices(i+1),2) = 1;
        if ismember(i,trialsLC_holdout)
            iFold(DM_UCOdata.endTrialIndices(i)+1:DM_UCOdata.endTrialIndices(i+1),1) = xvalFold+1;
        else
            iC = iC+1;
            iFoldLocal = mod(iC, xvalFold); if iFoldLocal == 0; iFoldLocal = xvalFold; end
            iFold(DM_UCOdata.endTrialIndices(i)+1:DM_UCOdata.endTrialIndices(i+1),1) = iFoldLocal;
        end
    end
    if ismember(i,trialsLOall)
        iUCO(DM_UCOdata.endTrialIndices(i)+1:DM_UCOdata.endTrialIndices(i+1),3) = 1;
        if ismember(i,trialsLO_holdout)
            iFold(DM_UCOdata.endTrialIndices(i)+1:DM_UCOdata.endTrialIndices(i+1),1) = xvalFold+1;
        else
            iO = iO+1;
            iFoldLocal = mod(iO, xvalFold); if iFoldLocal == 0; iFoldLocal = xvalFold; end
            iFold(DM_UCOdata.endTrialIndices(i)+1:DM_UCOdata.endTrialIndices(i+1),1) = iFoldLocal;
        end
    end
end

clear i j iC iFoldLocal iO iU

%% Use correct trial indices to pull out correct trials from big data matrices
% Get indices of correct trials form big matrix + Pull indices out of bigger matrices
iClong = logical(iUCO(:,2));
iTrialsC = iTrials(iClong);
iTrialsGlobalC = iTrialsGlobal(iClong);
iFoldC = iFold(iClong);

% Make small matrix and repeat steps above
iUCOsmall = iUCO(iTrialUCOSmall,:);
iCsmall = logical(iUCOsmall(:,2));
iTrialsCsmall = iTrials(iTrialUCOSmall(iCsmall));
iTrialsGlobalCsmall = iTrialsGlobal(iTrialUCOSmall(iCsmall));
iFoldCsmall = iFold(iTrialUCOSmall(iCsmall));

%% use PCA to reduce design matrix
pcaVar = inits.pcaVar;

% % Randomize DM to assess how many PCs to use bsBoxFields =
% fieldnames(trialBasis.bsboxUCO); nFields = zeros(size(bsBoxFields)); for
% iField = 1:numel(bsBoxFields)
%     nFields(iField) =
%     size(trialBasis.bsboxUCO(1).(bsBoxFields{iField}),2);
% end nFields = [0; cumsum(nFields)];
% 
% nFields = [0:290]';
% 
% nShuf = 500; DMexp_shuf = single(zeros(size(DM_UCO,2),nShuf)); rng(2349);
% for iShuf = 1:nShuf
%     DM_UCO_dummy = DM_UCO; shufSteps =
%     randi(size(DM_UCO,1),numel(nFields)-1,1); for jField =
%     1:numel(nFields)-1
%         DM_UCO_dummy(:,nFields(jField)+1:nFields(jField+1)) =
%         circshift(DM_UCO_dummy(:,nFields(jField)+1:nFields(jField+1)),shufSteps(jField),1);
%     end [~,~,~,~,expDummy] = pca(DM_UCO_dummy); DMexp_shuf(:,iShuf) =
%     expDummy;
% end DM_UCOpca = struct; [DM_UCOpca.coef, DM_UCOpca.score, DM_UCOpca.lat,
% ~, DM_UCOpca.exp] = pca(DM_UCO);
% 
% figure; plot(DM_UCOpca.exp,'-ok','LineWidth',2); hold on;
% plot(DMexp_shuf)
% 
% figure; plot(cumsum(DM_UCOpca.exp),'-ok','LineWidth',2); hold on;
% plot(cumsum(DMexp_shuf,1))

% Full dataset: long format - all predictors and delay periods
DM_UCOpca = struct;
[DM_UCOpca.coef, DM_UCOpca.score, DM_UCOpca.lat, ~, DM_UCOpca.exp] = pca(DM_UCO);
DM_UCOpca.pcaVar = pcaVar;
DM_UCOpca.nPCvar = sum(cumsum(DM_UCOpca.exp) < pcaVar)+1;
DM = DM_UCO; DMdata = DM_UCOdata; DMpca = DM_UCOpca; Ytrials = Ytrial_UCO; zFtrials = zFtrial_UCO;
GLMfull = v2struct(DM, DMdata, Ytrials, zFtrials, DMpca, iTrials, iTrialsGlobal, iUCO, iFold);
clear DM DMdata DMpca Ytrials

% Full dataset: short format - trimmed trials and predictors
DM_UCOSmallpca = struct;
[DM_UCOSmallpca.coef, DM_UCOSmallpca.score, DM_UCOSmallpca.lat, ~, DM_UCOSmallpca.exp] = pca(DM_UCOSmall);
DM_UCOSmallpca.pcaVar = pcaVar;
DM_UCOSmallpca.nPCvar = sum(cumsum(DM_UCOSmallpca.exp) < pcaVar)+1;
DM = DM_UCOSmall; DMdata = DM_UCOdataSmall; DMpca = DM_UCOSmallpca; Ytrials = Ytrial_UCOSmall; zFtrials = zFtrial_UCOSmall;
GLMfull_small = v2struct(DM, DMdata, Ytrials, zFtrials, DMpca, iTrials, iTrialsGlobal, iUCO, iFold);
clear DM DMdata DMpca Ytrials

% Correct trials only: long format - all predictors and delay periods
DM_Cpca = struct;
[DM_Cpca.coef, DM_Cpca.score, DM_Cpca.lat, ~, DM_Cpca.exp] = pca(DM_C);
DM_Cpca.pcaVar = pcaVar;
DM_Cpca.nPCvar = sum(cumsum(DM_Cpca.exp) < pcaVar)+1;
DM = DM_C; DMdata = DM_Cdata; DMpca = DM_Cpca; Ytrials = Ytrial_C; zFtrials = zFtrial_C;
GLMcor = v2struct(DM, DMdata, Ytrials, zFtrials, DMpca, iTrialsC, iTrialsGlobalC, iFoldC);
clear DM DMdata DMpca Ytrials

% Correct trials only: short format - trimmed trials and predictors
DM_CSmallpca = struct;
[DM_CSmallpca.coef, DM_CSmallpca.score, DM_CSmallpca.lat, ~, DM_CSmallpca.exp] = pca(DM_CSmall);
DM_CSmallpca.pcaVar = pcaVar;
DM_CSmallpca.nPCvar = sum(cumsum(DM_CSmallpca.exp) < pcaVar)+1;
DM = DM_CSmall; DMdata = DM_CdataSmall; DMpca = DM_CSmallpca; Ytrials = Ytrial_CSmall; zFtrials = zFtrial_CSmall;
GLMcor_small = v2struct(DM, DMdata, Ytrials, zFtrials, DMpca, iTrialsC, iTrialsGlobalC, iFoldC);
clear DM DMdata DMpca Ytrials

if savebool
    RglmFold = fullfile(glmFold,'glmRfiles');
    if ~exist(RglmFold,'dir'); mkdir(RglmFold); end
    save(fullfile(RglmFold,'GLMfull.mat'),'GLMfull');
    save(fullfile(RglmFold,'GLMfull_small.mat'),'GLMfull_small');
    save(fullfile(RglmFold,'GLMcor.mat'),'GLMcor');
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


%% Save maybe
if savebool
    % Full data (all trials)
    RglmFold_1 = fullfile(RglmFold,'FullData');
    if ~exist(RglmFold_1,'dir'); mkdir(RglmFold_1); end
    writematrix(iTrials,fullfile(RglmFold_1,'iTrials.csv'));
    writematrix(iTrialsGlobal,fullfile(RglmFold_1,'iTrialsGlobal.csv'));
    writematrix(iUCO,fullfile(RglmFold_1,'iUCO.csv'));
    writematrix(DM_UCO,fullfile(RglmFold_1,'DMfull.csv'));
    writematrix(DM_UCOpca.score,fullfile(RglmFold_1,'DMfullPCAscore.csv'));
    writematrix(DM_UCOpca.nPCvar,fullfile(RglmFold_1,'nPCvar.csv'));
    writematrix(Ytrial_UCO,fullfile(RglmFold_1,'Yfull.csv'));
    writematrix(zFtrial_UCO,fullfile(RglmFold_1,'zFfull.csv'));
    writematrix(iFold,fullfile(RglmFold_1,'iFold.csv'));
    writematrix(iFold(iFold <= xvalFold),fullfile(RglmFold_1,'iFoldTest.csv'));
    writematrix(iFold(iFold > xvalFold),fullfile(RglmFold_1,'iFoldHoldout.csv'));
    
    % Full data (all trials) - trimmed trials and predictors
    RglmFold_2 = fullfile(RglmFold,'FullData_small');
    if ~exist(RglmFold_2,'dir'); mkdir(RglmFold_2); end
    writematrix(iTrials(iTrialUCOSmall),fullfile(RglmFold_2,'iTrialsSmall.csv'));
    writematrix(iTrialsGlobal(iTrialUCOSmall),fullfile(RglmFold_2,'iTrialsGlobalSmall.csv'));
    writematrix(iUCO(iTrialUCOSmall),fullfile(RglmFold_2,'iUCOSmall.csv'));
    writematrix(DM_UCOSmall,fullfile(RglmFold_2,'DMfullSmall.csv'));
    writematrix(DM_UCOSmallpca.score,fullfile(RglmFold_2,'DMfullSmallPCAscore.csv'));
    writematrix(Ytrial_UCOSmall,fullfile(RglmFold_2,'YfullSmall.csv'));
    writematrix(zFtrial_UCOSmall,fullfile(RglmFold_1,'zFfullSmall.csv'));
    iFoldSmall = iFold(iTrialUCOSmall);
    writematrix(iFoldSmall,fullfile(RglmFold_2,'iFoldSmall.csv'));
    writematrix(iFoldSmall(iFoldSmall <= xvalFold),fullfile(RglmFold_2,'iFoldSmallTest.csv'));
    writematrix(iFoldSmall(iFoldSmall > xvalFold),fullfile(RglmFold_2,'iFoldSmallHoldout.csv'));
    
    
    % Correct trials only
    RglmFold_3 = fullfile(RglmFold,'CorData');
    if ~exist(RglmFold_3,'dir'); mkdir(RglmFold_3); end
    writematrix(iTrialsC,fullfile(RglmFold_3,'iTrialsC.csv'));
    writematrix(iTrialsGlobalC,fullfile(RglmFold_3,'iTrialsGlobalC.csv'));
    writematrix(DM_C,fullfile(RglmFold_3,'DMcor.csv'));
    writematrix(DM_Cpca.score,fullfile(RglmFold_3,'DMcorPCAscore.csv'));
    writematrix(DM_Cpca.nPCvar,fullfile(RglmFold_3,'nPCvar.csv'));
    writematrix(Ytrial_C,fullfile(RglmFold_3,'Ycor.csv'));
    writematrix(zFtrial_C,fullfile(RglmFold_3,'zFcor.csv'));
    writematrix(iFoldC,fullfile(RglmFold_3,'iFoldC.csv'));
    writematrix(iFoldC(iFoldC <= xvalFold),fullfile(RglmFold_3,'iFoldCTest.csv'));
    writematrix(iFoldC(iFoldC > xvalFold),fullfile(RglmFold_3,'iFoldCHoldout.csv'));
    
    % Correct trials only - trimmed trials and predictors
    RglmFold_4 = fullfile(RglmFold,'CorData_small');
    if ~exist(RglmFold_4,'dir'); mkdir(RglmFold_4); end
    writematrix(iTrialsCsmall,fullfile(RglmFold_4,'iTrialsCsmall.csv'));
    writematrix(iTrialsGlobalCsmall,fullfile(RglmFold_4,'iTrialsGlobalCsmall.csv'));
    writematrix(DM_CSmall,fullfile(RglmFold_4,'DMcorSmall.csv'));
    writematrix(DM_CSmallpca.score,fullfile(RglmFold_4,'DMcorSmallPCAscore.csv'));
    writematrix(DM_CSmallpca.nPCvar,fullfile(RglmFold_4,'nPCvar.csv'));
    writematrix(Ytrial_CSmall,fullfile(RglmFold_4,'YcorSmall.csv'));
    writematrix(zFtrial_CSmall,fullfile(RglmFold_4,'zFcorSmall.csv'));
    writematrix(iFoldCsmall,fullfile(RglmFold_4,'iFoldCsmall.csv'));
    writematrix(iFoldCsmall(iFoldCsmall <= xvalFold),fullfile(RglmFold_4,'iFoldCSmallTest.csv'));
    writematrix(iFoldCsmall(iFoldCsmall > xvalFold),fullfile(RglmFold_4,'iFoldCSmallHoldout.csv'));
    
end

end

