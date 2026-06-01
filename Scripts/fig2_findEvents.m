function eventsOut = fig2_findEvents(GLMsmall,GLMsig,expt,inits)
%% Figure out the holdout trials

iHoldout = unique(GLMsmall.iTrials(GLMsmall.iFold == max(GLMsmall.iFold)));
exptHoldout = expt(iHoldout);

% Note that this is the same as the settings in 'reducePredictorsBox.m' -
% check this fxn for correspondence
preClip = round(inits.imfs*(inits.prebuffsec-0.2));
if preClip < 0; preClip = 0; end
postClip = round(inits.imfs*(inits.postbuffsec-1));
if postClip < 0; postClip = 0; end

%% Iterate through
% Make empty arrays, will populate trial by trial
% trialON = []; movON = trialON; rewT = trialON; lickT = trialON;
trialON = false(size(GLMsig(1).yData));
movON = trialON;
rewT = trialON;
lickT = trialON;
wheelV = nan(size(GLMsig(1).yData));

PC1 = nan(size(GLMsig(1).yData));

trialEnd = 0;
for iTrial = 1:numel(exptHoldout)
    % Trial onset detection
    traceLocal = false(exptHoldout(iTrial).duration-preClip-postClip,1);
    if ~isempty(exptHoldout(iTrial).trialon)
        iLocal = exptHoldout(iTrial).trialon-exptHoldout(iTrial).start-preClip+1;
        iLocal(iLocal > numel(traceLocal)) = [];
        traceLocal(iLocal) = true;
    end
    trialON(trialEnd+1:trialEnd+numel(traceLocal),1) = traceLocal;
    PC1(trialEnd+1:trialEnd+numel(traceLocal),1) = exptHoldout(iTrial).spkData(preClip+1:end-postClip,1);
    
    % Mov onset detection
    traceLocal(:) = false;
    if ~isempty(exptHoldout(iTrial).movonset)
        iLocal = exptHoldout(iTrial).movonset-exptHoldout(iTrial).start-preClip+1;
        iLocal(iLocal > numel(traceLocal)) = [];
        traceLocal(iLocal) = true;
    end
    movON(trialEnd+1:trialEnd+numel(traceLocal),1) = traceLocal;
    
    % Rew onset detection
    traceLocal(:) = false;
    if ~isempty(exptHoldout(iTrial).reward)
        iLocal = exptHoldout(iTrial).reward-exptHoldout(iTrial).start-preClip+1;
        iLocal(iLocal > numel(traceLocal)) = [];
        traceLocal(iLocal) = true;
    end
    rewT(trialEnd+1:trialEnd+numel(traceLocal),1) = traceLocal;
    
    % Licking onset detection
    traceLocal(:) = false;
    if ~isempty(exptHoldout(iTrial).lick)
        iLocal = exptHoldout(iTrial).lick-exptHoldout(iTrial).start-preClip+1;
        iLocal(iLocal > numel(traceLocal)) = [];
        traceLocal(iLocal) = true;
    end
    lickT(trialEnd+1:trialEnd+numel(traceLocal),1) = traceLocal;
    
    wheelV(trialEnd+1:trialEnd+numel(traceLocal),1) = -exptHoldout(iTrial).wheelvel(preClip+1:end-postClip);
    
    trialEnd = trialEnd+numel(traceLocal);
end

%% Get raw data out and make traces and matrices!
yData = cat(2,GLMsig.yData)*inits.imfs;
yHat = cat(2,GLMsig.yHat)*inits.imfs;
yShufVis = cat(2,GLMsig.yHat_shufVis_proj_mu)*inits.imfs;
yShufMov = cat(2,GLMsig.yHat_shufMov_proj_mu)*inits.imfs;
yShufRew = cat(2,GLMsig.yHat_shufRew_proj_mu)*inits.imfs;
yShufLick = cat(2,GLMsig.yHat_shufLick_proj_mu)*inits.imfs;

% Trial ON - 200 ms before, 800 ms after
iON = find(trialON);
yData_mat = nan(inits.imfs,size(yData,2),numel(iON)); yHat_mat = yData_mat;
yShufVis_mat = yData_mat; yShufMov_mat = yData_mat;
yShufRew_mat = yData_mat; yShufLick_mat = yData_mat;
movMat = nan(inits.imfs,numel(iON)); rewMat = movMat; lickMat = movMat; velMat = movMat;
for iTrial = 1:numel(iON)
    iData = iON(iTrial)-inits.imfs/5:iON(iTrial)+inits.imfs/5*4-1;
    yData_mat(:,:,iTrial) = yData(iData,:);
    yHat_mat(:,:,iTrial) = yHat(iData,:);
    yShufVis_mat(:,:,iTrial) = yShufVis(iData,:);
    yShufMov_mat(:,:,iTrial) = yShufMov(iData,:);
    yShufRew_mat(:,:,iTrial) = yShufRew(iData,:);
    yShufLick_mat(:,:,iTrial) = yShufLick(iData,:);
    movMat(:,iTrial) = movON(iData);
    rewMat(:,iTrial) = rewT(iData);
    lickMat(:,iTrial) = lickT(iData);
    velMat(:,iTrial) = wheelV(iData);
end
iCells = zeros(size(GLMsig))';
for j = 1:numel(GLMsig); iCells(j) = ~isempty(GLMsig(j).GLMrankPCA); end
eventsOut(1).event = 'Vis';
eventsOut(1).tEvent = (-inits.imfs/5:inits.imfs/5*4-1)'/inits.imfs+1/inits.imfs/2;
eventsOut(1).iCells = find(iCells);
eventsOut(1).yDataMu = mean(yData_mat,3,'omitnan');
eventsOut(1).yDataSD = std(yData_mat,[],3,'omitnan');
eventsOut(1).yHatMu = mean(yHat_mat,3,'omitnan');
eventsOut(1).yHatSD = std(yHat_mat,[],3,'omitnan');
eventsOut(1).yShufVisMu = mean(yShufVis_mat,3,'omitnan');
eventsOut(1).yShufVisSD = std(yShufVis_mat,[],3,'omitnan');
eventsOut(1).yShufMovMu = mean(yShufMov_mat,3,'omitnan');
eventsOut(1).yShufMovSD = std(yShufMov_mat,[],3,'omitnan');
eventsOut(1).yShufRewMu = mean(yShufRew_mat,3,'omitnan');
eventsOut(1).yShufRewSD = std(yShufRew_mat,[],3,'omitnan');
eventsOut(1).yShufLickMu = mean(yShufLick_mat,3,'omitnan');
eventsOut(1).yShufLickSD = std(yShufLick_mat,[],3,'omitnan');
eventsOut(1).movMu = mean(movMat,2,'omitnan');
eventsOut(1).movSD = std(movMat,[],2,'omitnan');
eventsOut(1).rewMu = mean(rewMat,2,'omitnan');
eventsOut(1).rewSD = std(rewMat,[],2,'omitnan');
eventsOut(1).lickMu = mean(lickMat,2,'omitnan');
eventsOut(1).lickSD = std(lickMat,[],2,'omitnan');
eventsOut(1).velMu = mean(velMat,2,'omitnan');
eventsOut(1).velSD = std(velMat,[],2,'omitnan');
eventsOut(1).nEvents = numel(iON);

% Mov onset - 500 ms on each side
iMov = find(movON);
yData_mat = nan(inits.imfs,size(yData,2),numel(iMov)); yHat_mat = yData_mat;
yShufVis_mat = yData_mat; yShufMov_mat = yData_mat;
yShufRew_mat = yData_mat; yShufLick_mat = yData_mat;
movMat = nan(inits.imfs,numel(iMov)); rewMat = movMat; lickMat = movMat; velMat = movMat;
for iTrial = 1:numel(iMov)
    if iMov(iTrial)>inits.imfs/2
        iData = iMov(iTrial)-inits.imfs/2:iMov(iTrial)+inits.imfs/2-1;
        iMat = 1:inits.imfs;
    else        
        iData = 1:iMov(iTrial)+inits.imfs/2-1;
        iMat = inits.imfs/2-iMov(iTrial)+2:inits.imfs; % +2 to account for index offsets
    end
    yData_mat(iMat,:,iTrial) = yData(iData,:);
    yHat_mat(iMat,:,iTrial) = yHat(iData,:);
    yShufVis_mat(iMat,:,iTrial) = yShufVis(iData,:);
    yShufMov_mat(iMat,:,iTrial) = yShufMov(iData,:);
    yShufRew_mat(iMat,:,iTrial) = yShufRew(iData,:);
    yShufLick_mat(iMat,:,iTrial) = yShufLick(iData,:);
    movMat(iMat,iTrial) = movON(iData);
    rewMat(iMat,iTrial) = rewT(iData);
    lickMat(iMat,iTrial) = lickT(iData);
    velMat(iMat,iTrial) = wheelV(iData);
end
eventsOut(2).event = 'Mov';
eventsOut(2).tEvent = (-inits.imfs/2:inits.imfs/2-1)'/inits.imfs+1/inits.imfs/2;
eventsOut(2).yDataMu = mean(yData_mat,3,'omitnan');
eventsOut(2).yDataSD = std(yData_mat,[],3,'omitnan');
eventsOut(2).yHatMu = mean(yHat_mat,3,'omitnan');
eventsOut(2).yHatSD = std(yHat_mat,[],3,'omitnan');
eventsOut(2).yShufVisMu = mean(yShufVis_mat,3,'omitnan');
eventsOut(2).yShufVisSD = std(yShufVis_mat,[],3,'omitnan');
eventsOut(2).yShufMovMu = mean(yShufMov_mat,3,'omitnan');
eventsOut(2).yShufMovSD = std(yShufMov_mat,[],3,'omitnan');
eventsOut(2).yShufRewMu = mean(yShufRew_mat,3,'omitnan');
eventsOut(2).yShufRewSD = std(yShufRew_mat,[],3,'omitnan');
eventsOut(2).yShufLickMu = mean(yShufLick_mat,3,'omitnan');
eventsOut(2).yShufLickSD = std(yShufLick_mat,[],3,'omitnan');
eventsOut(2).movMu = mean(movMat,2,'omitnan');
eventsOut(2).movSD = std(movMat,[],2,'omitnan');
eventsOut(2).rewMu = mean(rewMat,2,'omitnan');
eventsOut(2).rewSD = std(rewMat,[],2,'omitnan');
eventsOut(2).lickMu = mean(lickMat,2,'omitnan');
eventsOut(2).lickSD = std(lickMat,[],2,'omitnan');
eventsOut(2).velMu = mean(velMat,2,'omitnan');
eventsOut(2).velSD = std(velMat,[],2,'omitnan');
eventsOut(2).nEvents = numel(iMov);

% RewT - 500 ms on each side
iRew = find(rewT);
yData_mat = nan(inits.imfs,size(yData,2),numel(iRew)); yHat_mat = yData_mat;
yShufVis_mat = yData_mat; yShufMov_mat = yData_mat;
yShufRew_mat = yData_mat; yShufLick_mat = yData_mat;
movMat = nan(inits.imfs,numel(iRew)); rewMat = movMat; lickMat = movMat; velMat = movMat;
for iTrial = 1:numel(iRew)
    iData = iRew(iTrial)-inits.imfs/2:iRew(iTrial)+inits.imfs/2-1;
    yData_mat(:,:,iTrial) = yData(iData,:);
    yHat_mat(:,:,iTrial) = yHat(iData,:);
    yShufVis_mat(:,:,iTrial) = yShufVis(iData,:);
    yShufMov_mat(:,:,iTrial) = yShufMov(iData,:);
    yShufRew_mat(:,:,iTrial) = yShufRew(iData,:);
    yShufLick_mat(:,:,iTrial) = yShufLick(iData,:);
    movMat(:,iTrial) = movON(iData);
    rewMat(:,iTrial) = rewT(iData);
    lickMat(:,iTrial) = lickT(iData);
    velMat(:,iTrial) = wheelV(iData);
end
eventsOut(3).event = 'Rew';
eventsOut(3).tEvent = (-inits.imfs/2:inits.imfs/2-1)'/inits.imfs+1/inits.imfs/2;
eventsOut(3).yDataMu = mean(yData_mat,3,'omitnan');
eventsOut(3).yDataSD = std(yData_mat,[],3,'omitnan');
eventsOut(3).yHatMu = mean(yHat_mat,3,'omitnan');
eventsOut(3).yHatSD = std(yHat_mat,[],3,'omitnan');
eventsOut(3).yShufVisMu = mean(yShufVis_mat,3,'omitnan');
eventsOut(3).yShufVisSD = std(yShufVis_mat,[],3,'omitnan');
eventsOut(3).yShufMovMu = mean(yShufMov_mat,3,'omitnan');
eventsOut(3).yShufMovSD = std(yShufMov_mat,[],3,'omitnan');
eventsOut(3).yShufRewMu = mean(yShufRew_mat,3,'omitnan');
eventsOut(3).yShufRewSD = std(yShufRew_mat,[],3,'omitnan');
eventsOut(3).yShufLickMu = mean(yShufLick_mat,3,'omitnan');
eventsOut(3).yShufLickSD = std(yShufLick_mat,[],3,'omitnan');
eventsOut(3).movMu = mean(movMat,2,'omitnan');
eventsOut(3).movSD = std(movMat,[],2,'omitnan');
eventsOut(3).rewMu = mean(rewMat,2,'omitnan');
eventsOut(3).rewSD = std(rewMat,[],2,'omitnan');
eventsOut(3).lickMu = mean(lickMat,2,'omitnan');
eventsOut(3).lickSD = std(lickMat,[],2,'omitnan');
eventsOut(3).velMu = mean(velMat,2,'omitnan');
eventsOut(3).velSD = std(velMat,[],2,'omitnan');
eventsOut(3).nEvents = numel(iRew);

% LickT - 500 ms on each side
iLick = find(lickT); iLick(iLick > size(yData,1)-inits.imfs/2) = [];
yData_mat = nan(inits.imfs,size(yData,2),numel(iLick)); yHat_mat = yData_mat;
yShufVis_mat = yData_mat; yShufMov_mat = yData_mat;
yShufRew_mat = yData_mat; yShufLick_mat = yData_mat;
movMat = nan(inits.imfs,numel(iLick)); rewMat = movMat; lickMat = movMat; velMat = movMat;
for iTrial = 1:numel(iLick)
    iData = iLick(iTrial)-inits.imfs/2:iLick(iTrial)+inits.imfs/2-1;
    yData_mat(:,:,iTrial) = yData(iData,:);
    yHat_mat(:,:,iTrial) = yHat(iData,:);
    yShufVis_mat(:,:,iTrial) = yShufVis(iData,:);
    yShufMov_mat(:,:,iTrial) = yShufMov(iData,:);
    yShufRew_mat(:,:,iTrial) = yShufRew(iData,:);
    yShufLick_mat(:,:,iTrial) = yShufLick(iData,:);
    movMat(:,iTrial) = movON(iData);
    rewMat(:,iTrial) = rewT(iData);
    lickMat(:,iTrial) = lickT(iData);
    velMat(:,iTrial) = wheelV(iData);
end
eventsOut(4).event = 'Lick';
eventsOut(4).tEvent = (-inits.imfs/2:inits.imfs/2-1)'/inits.imfs+1/inits.imfs/2;
eventsOut(4).yDataMu = mean(yData_mat,3,'omitnan');
eventsOut(4).yDataSD = std(yData_mat,[],3,'omitnan');
eventsOut(4).yHatMu = mean(yHat_mat,3,'omitnan');
eventsOut(4).yHatSD = std(yHat_mat,[],3,'omitnan');
eventsOut(4).yShufVisMu = mean(yShufVis_mat,3,'omitnan');
eventsOut(4).yShufVisSD = std(yShufVis_mat,[],3,'omitnan');
eventsOut(4).yShufMovMu = mean(yShufMov_mat,3,'omitnan');
eventsOut(4).yShufMovSD = std(yShufMov_mat,[],3,'omitnan');
eventsOut(4).yShufRewMu = mean(yShufRew_mat,3,'omitnan');
eventsOut(4).yShufRewSD = std(yShufRew_mat,[],3,'omitnan');
eventsOut(4).yShufLickMu = mean(yShufLick_mat,3,'omitnan');
eventsOut(4).yShufLickSD = std(yShufLick_mat,[],3,'omitnan');
eventsOut(4).movMu = mean(movMat,2,'omitnan');
eventsOut(4).movSD = std(movMat,[],2,'omitnan');
eventsOut(4).rewMu = mean(rewMat,2,'omitnan');
eventsOut(4).rewSD = std(rewMat,[],2,'omitnan');
eventsOut(4).lickMu = mean(lickMat,2,'omitnan');
eventsOut(4).lickSD = std(lickMat,[],2,'omitnan');
eventsOut(4).velMu = mean(velMat,2,'omitnan');
eventsOut(4).velSD = std(velMat,[],2,'omitnan');
eventsOut(4).nEvents = numel(iLick);

end