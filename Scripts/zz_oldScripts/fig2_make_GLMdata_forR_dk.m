function [GLMdata, GLMholdout,DMdataBoxSmall] = fig2_make_GLMdata_forR_dk(DMboxdata0, trialbasis, spkNames, xvalfold,trialout,savebool)
% function [GLMdata, GLMholdout,DMdataBoxSmall] = fig2_make_GLMdata_forR_dk(DMdata, DMboxdata, trialbasis, rewbasis, mzones, xvalfold,trialout)

trialsLU0 = trialout.trialsLU0;
trialsLC0 = trialout.trialsLC0;
trialsLO0 = trialout.trialsLO0;

%% Get indices of spk variables 
% Assume that the spk vars are all at the end!!!!
% Note that we probably won't use these
% DMspki = [0; cumsum(DMdata.nBases)];
% DMspki = DMspki([~DMdata.behavfields; true]);
DMboxspki = [0; cumsum(DMboxdata0.nBases)];
DMboxspki = DMboxspki([~DMboxdata0.behavfields; true]);

trialsALL0 = cat(1,trialbasis.expt.trialnum);
trialsGOOD = ismember(trialsALL0,[trialsLU0; trialsLC0; trialsLO0]);

trialbasis.expt = trialbasis.expt(trialsGOOD);
trialbasis.bs = trialbasis.bs(trialsGOOD);
trialbasis.bsbox = trialbasis.bsbox(trialsGOOD);

trialsALL = cat(1,trialbasis.expt.trialnum);
trialsLUall = find(ismember(trialsALL,trialsLU0));
trialsLCall = find(ismember(trialsALL,trialsLC0));
trialsLOall = find(ismember(trialsALL,trialsLO0));
nTrials = numel(trialsALL);

starti = 2;

trialsLU_holdout = trialsLUall(starti:3:end);
trialsLC_holdout = trialsLCall(starti:3:end);
trialsLO_holdout = trialsLOall(starti:3:end);

trialsLU = trialsLUall(~ismember(trialsLUall,trialsLU_holdout));
trialsLC = trialsLCall(~ismember(trialsLCall,trialsLC_holdout));
trialsLO = trialsLOall(~ismember(trialsLOall,trialsLO_holdout));

% iTrials0 = (1:numel(trialbasis.expt))';
% iHoldTrials = sort([trialsLU_holdout; trialsLC_holdout; trialsLO_holdout]);
% iTrials = iTrials0(~ismember(iTrials0,iHoldTrials));

% iRews0 = (1:numel(rewbasis.expt))';
% iHoldRews = iRews0(starti:3:end);
% iRews = iRews0(~ismember(iRews0,iHoldRews));

clear trialsALL0 trialsLU0 trialsLC0 trialsLO0 iTrials0 iRews0 trialout

%% Gather info to feed to glmnet - R
% Remove reward fields from design matrix
bsbox2 = rmfield(trialbasis.bsbox,{'CueBeep','CueReward','RandReward'});
trialbasis.bsbox2 = bsbox2; clear bsbox2;

% Make full design matrix
[DMboxAll,DMboxAllData] = makeDesignMatrix_dk(trialbasis, find(trialsGOOD), 'bsbox2');
% Trim off predictors that are during spont activity/witholds
[DMdataBoxSmall,iBoxSmall] = reducePredictorsBox(DMboxAllData);
DMboxSmall = DMboxAll(:,iBoxSmall);

% Make matrix of spiking data
YtrialAll = zeros(size(DMboxAll,1),numel(spkNames));
for j = 1:numel(spkNames)        
    YtrialAll(:,j) = getSpikes_dk(trialbasis,find(trialsGOOD),spkNames{j});
end

% Make vectors that denote trial outcome, x-val fold, holdout status and
% local and global trial number
iUCO = zeros(size(DMboxAll,1),3);
iTrials = zeros(size(DMboxAll,1),1); iFold = iTrials; iTrialsGlobal = iTrials;
iU = 0; iC = 0; iO = 0;
for i = 1:nTrials
    iTrials(DMdataBoxSmall.endTrialIndices(i)+1:DMdataBoxSmall.endTrialIndices(i+1),1) = i;
    iTrialsGlobal(DMdataBoxSmall.endTrialIndices(i)+1:DMdataBoxSmall.endTrialIndices(i+1),1) = trialsALL(i);
    if ismember(i,trialsLUall)
        iUCO(DMdataBoxSmall.endTrialIndices(i)+1:DMdataBoxSmall.endTrialIndices(i+1),1) = 1;
        if ismember(i,trialsLU_holdout)
            iFold(DMdataBoxSmall.endTrialIndices(i)+1:DMdataBoxSmall.endTrialIndices(i+1),1) = xvalfold+1;
        else
            iU = iU+1;
            iFoldLocal = mod(iU, xvalfold); if iFoldLocal == 0; iFoldLocal = xvalfold; end
            iFold(DMdataBoxSmall.endTrialIndices(i)+1:DMdataBoxSmall.endTrialIndices(i+1),1) = iFoldLocal;
        end
    end
    if ismember(i,trialsLCall)
        iUCO(DMdataBoxSmall.endTrialIndices(i)+1:DMdataBoxSmall.endTrialIndices(i+1),2) = 1;
        if ismember(i,trialsLC_holdout)
            iFold(DMdataBoxSmall.endTrialIndices(i)+1:DMdataBoxSmall.endTrialIndices(i+1),1) = xvalfold+1;
        else
            iC = iC+1;
            iFoldLocal = mod(iC, xvalfold); if iFoldLocal == 0; iFoldLocal = xvalfold; end
            iFold(DMdataBoxSmall.endTrialIndices(i)+1:DMdataBoxSmall.endTrialIndices(i+1),1) = iFoldLocal;
        end
    end
    if ismember(i,trialsLOall)
        iUCO(DMdataBoxSmall.endTrialIndices(i)+1:DMdataBoxSmall.endTrialIndices(i+1),3) = 1;
        if ismember(i,trialsLO_holdout)
            iFold(DMdataBoxSmall.endTrialIndices(i)+1:DMdataBoxSmall.endTrialIndices(i+1),1) = xvalfold+1;
        else
            iO = iO+1;
            iFoldLocal = mod(iO, xvalfold); if iFoldLocal == 0; iFoldLocal = xvalfold; end
            iFold(DMdataBoxSmall.endTrialIndices(i)+1:DMdataBoxSmall.endTrialIndices(i+1),1) = iFoldLocal;
        end
    end
end




clear i j iC iFoldLocal iO iU

if savebool
    writematrix(iTrials,'iTrials.csv');
    writematrix(iTrialsGlobal,'iTrialsGlobal.csv');
    writematrix(iUCO,'iUCO.csv');
    writematrix(iFold,'iFold.csv');
    writematrix(DMboxAll,'DMboxAll.csv');
    writematrix(DMboxSmall,'DMboxSmall.csv');
    writematrix(YtrialAll,'YtrialAll.csv');
    
    % Split into test and holdout
    writematrix(iFold(iFold <= xvalfold),'iFoldTest.csv');
    writematrix(DMboxAll(iFold <= xvalfold,:),'DMboxAllTest.csv');
    writematrix(DMboxSmall(iFold <= xvalfold,:),'DMboxSmallTest.csv');
    writematrix(YtrialAll(iFold <= xvalfold,:),'YtrialAllTest.csv');
    
    writematrix(iFold(iFold > xvalfold),'iFoldHoldout.csv');
    writematrix(DMboxAll(iFold > xvalfold,:),'DMboxAllHoldout.csv');
    writematrix(DMboxSmall(iFold > xvalfold,:),'DMboxSmallHoldout.csv');
    writematrix(YtrialAll(iFold > xvalfold,:),'YtrialAllHoldout.csv');
end

