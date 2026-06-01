function seshOut = fig3_movAnalysis(sesh, imfs, vrfs, figFold,plotbool, savebool)
%% This parametrizes movement in each session by:
%     1. UCO outcome
%     2. Displacement
%     3. Max velocity
%     4. Max accelaration
%     5. correctness
%     6. Activity amount

%% Set out variables that will be used to organize data
nTrials = numel(sesh.iTrial);
spkMov = smoothdata(sesh.spkMov(:,sesh.iPCmov,:),1,'gaussian',3);
spkMovBin = smoothdata(single(logical(sesh.spkMov(:,sesh.iPCmov,:))),1,'gaussian',3);
seshNames = fieldnames(sesh);
seshInfo = rmfield(sesh,seshNames(4:end));

iUCO = sesh.iUCO;               % 1. Outcome
maxDisp = zeros(nTrials,1);     % 2. Max displacement
maxVel = maxDisp;               % 3. Max velocity
maxAccel = maxDisp;             % 4. Max accelaration
RmovC = sesh.RmovC;             % 5. Correctness
maxSpk = maxDisp;               % 6. Max activity
maxSpkBin = maxDisp;            % 7. Max sychrony

coeff = pca(mean(spkMov(end/2-imfs/2+1:end/2+imfs/2,:,:),3));
[~, iSort] = sort(coeff(:,1),'descend');
spkMov = spkMov(:,iSort,:);
spkMovBin = spkMovBin(:,iSort,:);

seshOut = seshInfo;
seshOut.spkMovTrial = squeeze(imfs*max(mean(spkMov(end/2-imfs/2+1:end/2,:,:),2)));
seshOut.spkMovTrace = imfs*mean(mean(spkMov,3),2);
%% Tabulate parameters that require new calculations

for iTrial = 1:nTrials
    maxDisp(iTrial) = max(sesh.xMov(end/2+1:end/2+vrfs,iTrial));
    maxVel(iTrial) = max(sesh.vMov(end/2+1:end/2+vrfs,iTrial));
    maxAccel(iTrial) = max(diff(sesh.vMov(end/2+1:end/2+vrfs,iTrial)))*100;
    maxSpk(iTrial) = imfs*max(mean(spkMov(end/2-imfs/2+1:end/2,:,iTrial),2));
    maxSpkBin(iTrial) = imfs*max(mean(spkMovBin(end/2-imfs/2+1:end/2,:,iTrial),2));

end

%% 1. Outcome
spkMat = zeros(size(spkMov,1),size(spkMov,2),3);
xMat = cell(3,1);
vMat = xMat;
colMat = [158 31 99; 0 0 0; 194 117 20]/255;
for iOutcome = 1:3
    spkMat(:,:,iOutcome) = mean(spkMov(:,:,iUCO(:,iOutcome)),3);
    xMat{iOutcome,1} = sesh.xMov(:,iUCO(:,iOutcome));
    vMat{iOutcome,1} = sesh.vMov(:,iUCO(:,iOutcome));
%     seshOut.Data.Outcome.spkMax(iOutcome,1) = imfs*mean(mean(spkMat(end/2-imfs/2+1:end/2,:,iOutcome),2));
%     seshOut.Data.Outcome.spkMaxCell(:,iOutcome) = imfs*mean(spkMat(end/2-imfs/2+1:end/2,:,iOutcome));
    seshOut.Data.Outcome.spkMax(iOutcome,1) = imfs*max(mean(spkMat(end/2-imfs/2+1:end/2,:,iOutcome),2));
    seshOut.Data.Outcome.spkMaxCell(:,iOutcome) = imfs*max(spkMat(end/2-imfs/2+1:end/2,:,iOutcome));
    seshOut.Data.Outcome.xMax(iOutcome,1) = max(mean(sesh.xMov(end/2+1:end/2+vrfs,iUCO(:,iOutcome)),2));
    seshOut.Data.Outcome.vMax(iOutcome,1) = max(mean(sesh.vMov(end/2+1:end/2+vrfs,iUCO(:,iOutcome)),2));
end
seshOut.outcomeSpk = seshOut.Data.Outcome.spkMax;
seshOut.outcomeSpkCell = seshOut.Data.Outcome.spkMaxCell;
seshInfo.splitType = 'Outcome';
seshInfo.range = {'Undershoot';'Correct';'Overshoot'};
if plotbool; figOutcome = fig3_movAnalysis_plot(imfs*spkMat,xMat,vMat,imfs,vrfs,seshInfo,colMat); end

%% 2. Max displacement
nDisps = 5;
spkMat = zeros(size(spkMov,1),size(spkMov,2),nDisps);
xMat = cell(nDisps,1);
vMat = xMat;
colMat = repmat(linspace(.8,0,5),3,1)';
[~,iMaxDisp] = sort(maxDisp,'ascend');
edges = round(linspace(1,nTrials,nDisps+1));
seshInfo.range = [];
for iDisp = 1:nDisps
    spkMat(:,:,iDisp) = mean(spkMov(:,:,iMaxDisp(edges(iDisp):edges(iDisp+1))),3);
    xMat{iDisp,1} = sesh.xMov(:,iMaxDisp(edges(iDisp):edges(iDisp+1)));
    vMat{iDisp,1} = sesh.vMov(:,iMaxDisp(edges(iDisp):edges(iDisp+1)));
%     seshOut.Data.xMax.spkMax(iDisp,1) = imfs*mean(mean(spkMat(end/2-imfs/2+1:end/2,:,iDisp),2));
%     seshOut.Data.xMax.spkMaxCell(:,iDisp) = imfs*mean(spkMat(end/2-imfs/2+1:end/2,:,iDisp));
    seshOut.Data.xMax.spkMax(iDisp,1) = imfs*max(mean(spkMat(end/2-imfs/2+1:end/2,:,iDisp),2));
    seshOut.Data.xMax.spkMaxCell(:,iDisp) = imfs*max(spkMat(end/2-imfs/2+1:end/2,:,iDisp));
    seshOut.Data.xMax.xMax(iDisp,1) = max(mean(sesh.xMov(end/2+1:end/2+vrfs,iMaxDisp(edges(iDisp):edges(iDisp+1))),2));
    seshOut.Data.xMax.vMax(iDisp,1) = max(mean(sesh.vMov(end/2+1:end/2+vrfs,iMaxDisp(edges(iDisp):edges(iDisp+1))),2));
    seshInfo.range{iDisp,1} = median(maxDisp(iMaxDisp(edges(iDisp):edges(iDisp+1)))); %maxDisp(iMaxDisp(edges(iDisp)));
end
seshOut.xMaxSpk = seshOut.Data.xMax.spkMax;
seshOut.xMaxSpkCell = seshOut.Data.xMax.spkMaxCell;
seshInfo.splitType = 'maxDisplacement';
if plotbool; figDisplacement = fig3_movAnalysis_plot(imfs*spkMat,xMat,vMat,imfs,vrfs,seshInfo,colMat); end

%% 3. Max velocity
nDisps = 5;
spkMat = zeros(size(spkMov,1),size(spkMov,2),nDisps);
xMat = cell(nDisps,1);
vMat = xMat;
colMat = repmat(linspace(.8,0,5),3,1)';
[~,iMaxVel] = sort(maxVel,'ascend');
edges = round(linspace(1,nTrials,nDisps+1));
seshInfo.range = [];
for iDisp = 1:nDisps
    spkMat(:,:,iDisp) = mean(spkMov(:,:,iMaxVel(edges(iDisp):edges(iDisp+1))),3);
    xMat{iDisp,1} = sesh.xMov(:,iMaxVel(edges(iDisp):edges(iDisp+1)));
    vMat{iDisp,1} = sesh.vMov(:,iMaxVel(edges(iDisp):edges(iDisp+1)));
%     seshOut.Data.vMax.spkMax(iDisp,1) = imfs*mean(mean(spkMat(end/2-imfs/2+1:end/2,:,iDisp),2));
%     seshOut.Data.vMax.spkMaxCell(:,iDisp) = imfs*mean(spkMat(end/2-imfs/2+1:end/2,:,iDisp));
    seshOut.Data.vMax.spkMax(iDisp,1) = imfs*max(mean(spkMat(end/2-imfs/2+1:end/2,:,iDisp),2));
    seshOut.Data.vMax.spkMaxCell(:,iDisp) = imfs*max(spkMat(end/2-imfs/2+1:end/2,:,iDisp));
    seshOut.Data.vMax.xMax(iDisp,1) = max(mean(sesh.xMov(end/2+1:end/2+vrfs,iMaxVel(edges(iDisp):edges(iDisp+1))),2));
    seshOut.Data.vMax.vMax(iDisp,1) = max(mean(sesh.vMov(end/2+1:end/2+vrfs,iMaxVel(edges(iDisp):edges(iDisp+1))),2));
    seshInfo.range{iDisp,1} = median(maxVel(iMaxVel(edges(iDisp):edges(iDisp+1)))); % maxVel(iMaxVel(edges(iDisp)));
end
seshOut.vMaxSpk = seshOut.Data.vMax.spkMax;
seshOut.vMaxSpkCell = seshOut.Data.vMax.spkMaxCell;
seshInfo.splitType = 'maxVelocity';
if plotbool; figVelocity = fig3_movAnalysis_plot(imfs*spkMat,xMat,vMat,imfs,vrfs,seshInfo,colMat); end

%% 4. Max accelaration
nDisps = 5;
spkMat = zeros(size(spkMov,1),size(spkMov,2),nDisps);
xMat = cell(nDisps,1);
vMat = xMat;
colMat = repmat(linspace(.8,0,5),3,1)';
[~,iMaxAccel] = sort(maxAccel,'ascend');
edges = round(linspace(1,nTrials,nDisps+1));
seshInfo.range = [];
for iDisp = 1:nDisps
    spkMat(:,:,iDisp) = mean(spkMov(:,:,iMaxAccel(edges(iDisp):edges(iDisp+1))),3);
    xMat{iDisp,1} = sesh.xMov(:,iMaxAccel(edges(iDisp):edges(iDisp+1)));
    vMat{iDisp,1} = sesh.vMov(:,iMaxAccel(edges(iDisp):edges(iDisp+1)));
%     seshOut.Data.aMax.spkMax(iDisp,1) = imfs*mean(mean(spkMat(end/2-imfs/2+1:end/2,:,iDisp),2));
%     seshOut.Data.aMax.spkMaxCell(:,iDisp) = imfs*mean(spkMat(end/2-imfs/2+1:end/2,:,iDisp));
    seshOut.Data.aMax.spkMax(iDisp,1) = imfs*max(mean(spkMat(end/2-imfs/2+1:end/2,:,iDisp),2));
    seshOut.Data.aMax.spkMaxCell(:,iDisp) = imfs*max(spkMat(end/2-imfs/2+1:end/2,:,iDisp));
    seshOut.Data.aMax.xMax(iDisp,1) = max(mean(sesh.xMov(end/2+1:end/2+vrfs,iMaxAccel(edges(iDisp):edges(iDisp+1))),2));
    seshOut.Data.aMax.vMax(iDisp,1) = max(mean(sesh.vMov(end/2+1:end/2+vrfs,iMaxAccel(edges(iDisp):edges(iDisp+1))),2));
    seshInfo.range{iDisp,1} = median(maxAccel(iMaxAccel(edges(iDisp):edges(iDisp+1))));
end
seshOut.aMaxSpk = seshOut.Data.aMax.spkMax;
seshOut.aMaxSpkCell = seshOut.Data.aMax.spkMaxCell;
seshInfo.splitType = 'maxAccelaration';
if plotbool; figAccel = fig3_movAnalysis_plot(imfs*spkMat,xMat,vMat,imfs,vrfs,seshInfo,colMat); end

%% 5. Max correctness
nDisps = 5;
spkMat = zeros(size(spkMov,1),size(spkMov,2),nDisps);
xMat = cell(nDisps,1);
vMat = xMat;
colMat = repmat(linspace(.8,0,5),3,1)';
[~,iCorr] = sort(RmovC,'ascend');
edges = round(linspace(1,nTrials,nDisps+1));
seshInfo.range = [];
for iDisp = 1:nDisps
    spkMat(:,:,iDisp) = mean(spkMov(:,:,iCorr(edges(iDisp):edges(iDisp+1))),3);
    xMat{iDisp,1} = sesh.xMov(:,iCorr(edges(iDisp):edges(iDisp+1)));
    vMat{iDisp,1} = sesh.vMov(:,iCorr(edges(iDisp):edges(iDisp+1)));
%     seshOut.Data.Rcorr.spkMax(iDisp,1) = imfs*mean(mean(spkMat(end/2-imfs/2+1:end/2,:,iDisp),2));
%     seshOut.Data.Rcorr.spkMaxCell(:,iDisp) = imfs*mean(spkMat(end/2-imfs/2+1:end/2,:,iDisp));
    seshOut.Data.Rcorr.spkMax(iDisp,1) = imfs*max(mean(spkMat(end/2-imfs/2+1:end/2,:,iDisp),2));
    seshOut.Data.Rcorr.spkMaxCell(:,iDisp) = imfs*max(spkMat(end/2-imfs/2+1:end/2,:,iDisp));
    seshOut.Data.Rcorr.xMax(iDisp,1) = max(mean(sesh.xMov(end/2+1:end/2+vrfs,iCorr(edges(iDisp):edges(iDisp+1))),2));
    seshOut.Data.Rcorr.vMax(iDisp,1) = max(mean(sesh.vMov(end/2+1:end/2+vrfs,iCorr(edges(iDisp):edges(iDisp+1))),2));
    seshInfo.range{iDisp,1} = median(RmovC(iCorr(edges(iDisp):edges(iDisp+1)))); % RmovC(iCorr(edges(iDisp)));
end
seshOut.RcorrSpk = seshOut.Data.Rcorr.spkMax;
seshOut.RcorrSpkCell = seshOut.Data.Rcorr.spkMaxCell;
seshInfo.splitType = 'maxCorrectness';
if plotbool; figRcorr = fig3_movAnalysis_plot(imfs*spkMat,xMat,vMat,imfs,vrfs,seshInfo,colMat); end

%% 5b. Max correctness - by actual correlation
nDisps = 3;
spkMat = zeros(size(spkMov,1),size(spkMov,2),nDisps);
xMat = cell(nDisps,1);
vMat = xMat;
colMat = repmat(linspace(.8,0,5),3,1)';
[~,iCorr] = sort(RmovC,'ascend');

dispBins = [0 .5, .9];
% dispBins = [0,linspace(0.5,0.9,nDisps-1)];
edges = zeros(nDisps+1,1); edges(end) = nTrials;
for iDispMid = 1:numel(dispBins)
    edges(iDispMid) = find(RmovC(iCorr) > dispBins(iDispMid),1);
end
seshInfo.range = [];
for iDisp = 1:nDisps
    spkMat(:,:,iDisp) = mean(spkMov(:,:,iCorr(edges(iDisp):edges(iDisp+1))),3);
    xMat{iDisp,1} = sesh.xMov(:,iCorr(edges(iDisp):edges(iDisp+1)));
    vMat{iDisp,1} = sesh.vMov(:,iCorr(edges(iDisp):edges(iDisp+1)));
%     seshOut.Data.RcorrB.spkMax(iDisp,1) = imfs*mean(mean(spkMat(end/2-imfs*.3+1:end/2,:,iDisp),2));
%     seshOut.Data.RcorrB.spkMaxCell(:,iDisp) = imfs*mean(spkMat(end/2-imfs*.3+1:end/2,:,iDisp));
    seshOut.Data.RcorrB.spkMax(iDisp,1) = imfs*max(mean(spkMat(end/2-imfs/2+1:end/2,:,iDisp),2));
    seshOut.Data.RcorrB.spkMaxCell(:,iDisp) = imfs*max(spkMat(end/2-imfs/2+1:end/2,:,iDisp));
    seshOut.Data.RcorrB.xMax(iDisp,1) = max(mean(sesh.xMov(end/2+1:end/2+vrfs,iCorr(edges(iDisp):edges(iDisp+1))),2));
    seshOut.Data.RcorrB.vMax(iDisp,1) = max(mean(sesh.vMov(end/2+1:end/2+vrfs,iCorr(edges(iDisp):edges(iDisp+1))),2));
    seshInfo.range{iDisp,1} = median(RmovC(iCorr(edges(iDisp):edges(iDisp+1)))); % RmovC(iCorr(edges(iDisp)));
end
if numel(seshOut.Data.RcorrB.spkMax) < 5
    seshOut.Data.RcorrB.spkMax(end+1:5) = nan;
end
seshOut.RcorrSpkB = seshOut.Data.RcorrB.spkMax;
seshOut.RcorrSpkBCell = seshOut.Data.RcorrB.spkMaxCell;
seshInfo.splitType = 'maxCorrectness version B';
if plotbool; figRcorrB = fig3_movAnalysis_plot(imfs*spkMat,xMat,vMat,imfs,vrfs,seshInfo,colMat); end
% figRcorrB = fig3_movAnalysis_plot(imfs*spkMat,xMat,vMat,imfs,vrfs,seshInfo,colMat);
%% 6. Max activity
nDisps = 5;
spkMat = zeros(size(spkMov,1),size(spkMov,2),nDisps);
xMat = cell(nDisps,1);
vMat = xMat;
colMat = repmat(linspace(.8,0,5),3,1)';
[~,iSpk] = sort(maxSpk,'ascend');
edges = round(linspace(1,nTrials,nDisps+1));
seshInfo.range = [];
for iDisp = 1:nDisps
    spkMat(:,:,iDisp) = mean(spkMov(:,:,iSpk(edges(iDisp):edges(iDisp+1))),3);
    xMat{iDisp,1} = sesh.xMov(:,iSpk(edges(iDisp):edges(iDisp+1)));
    vMat{iDisp,1} = sesh.vMov(:,iSpk(edges(iDisp):edges(iDisp+1)));
    seshOut.Data.spkMax.spkMax(iDisp,1) = imfs*mean(mean(spkMat(end/2-imfs/2+1:end/2,:,iDisp),2));
%     seshOut.Data.spkMax.spkMax(iDisp,1) = imfs*mean(mean(spkMat(end/2-imfs/2+1:end/2,:,iDisp),2));
    seshOut.Data.spkMax.xMax(iDisp,1) = max(mean(sesh.xMov(end/2+1:end/2+vrfs,iSpk(edges(iDisp):edges(iDisp+1))),2));
    seshOut.Data.spkMax.vMax(iDisp,1) = max(mean(sesh.vMov(end/2+1:end/2+vrfs,iSpk(edges(iDisp):edges(iDisp+1))),2));
    seshInfo.range{iDisp,1} = median(maxSpk(iSpk(edges(iDisp):edges(iDisp+1)))); %maxSpk(iSpk(edges(iDisp)));
end
seshOut.spkMaxvMax = seshOut.Data.spkMax.vMax;
seshInfo.splitType = 'maxActivity';
if plotbool; figActivity = fig3_movAnalysis_plot(imfs*spkMat,xMat,vMat,imfs,vrfs,seshInfo,colMat); end

%% 7. Max activity (binary)
nDisps = 5;
spkMat = zeros(size(spkMov,1),size(spkMov,2),nDisps);
xMat = cell(nDisps,1);
vMat = xMat;
colMat = repmat(linspace(.8,0,5),3,1)';
[~,iSpk] = sort(maxSpkBin,'ascend');
edges = round(linspace(1,nTrials,nDisps+1));
seshInfo.range = [];
for iDisp = 1:nDisps
    spkMat(:,:,iDisp) = mean(spkMov(:,:,iSpk(edges(iDisp):edges(iDisp+1))),3);
    xMat{iDisp,1} = sesh.xMov(:,iSpk(edges(iDisp):edges(iDisp+1)));
    vMat{iDisp,1} = sesh.vMov(:,iSpk(edges(iDisp):edges(iDisp+1)));
    seshOut.Data.syncMax.spkMax(iDisp,1) = imfs*mean(mean(spkMat(end/2-imfs/2+1:end/2,:,iDisp),2));
%     seshOut.Data.syncMax.spkMax(iDisp,1) = imfs*max(mean(spkMat(end/2-imfs/2+1:end/2,:,iDisp),2));
    seshOut.Data.syncMax.xMax(iDisp,1) = max(mean(sesh.xMov(end/2+1:end/2+vrfs,iSpk(edges(iDisp):edges(iDisp+1))),2));
    seshOut.Data.syncMax.vMax(iDisp,1) = max(mean(sesh.vMov(end/2+1:end/2+vrfs,iSpk(edges(iDisp):edges(iDisp+1))),2));
    seshInfo.range{iDisp,1} = median(maxSpkBin(iSpk(edges(iDisp):edges(iDisp+1)))); % maxSpkBin(iSpk(edges(iDisp)));
end
seshOut.syncMaxvMax = seshOut.Data.syncMax.vMax;
seshInfo.splitType = 'maxSynchrony';
if plotbool; figSynchrony = fig3_movAnalysis_plot(imfs*spkMat,xMat,vMat,imfs,vrfs,seshInfo,colMat); end


%% Save maybe
if savebool && plotbool
    figname = sprintf('%s_01%s_%s',seshInfo.fov,'UCO',seshInfo.name);
    savefig(figOutcome, fullfile(figFold,[figname,'.fig']));
    saveas(figOutcome,fullfile(figFold,[figname,'.png']));
    print(figOutcome,fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    
    figname = sprintf('%s_02%s_%s',seshInfo.fov,'xMax',seshInfo.name);
    savefig(figDisplacement, fullfile(figFold,[figname,'.fig']));
    saveas(figDisplacement,fullfile(figFold,[figname,'.png']));
    print(figDisplacement,fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    
    figname = sprintf('%s_03%s_%s',seshInfo.fov,'vMax',seshInfo.name);
    savefig(figVelocity, fullfile(figFold,[figname,'.fig']));
    saveas(figVelocity,fullfile(figFold,[figname,'.png']));
    print(figVelocity,fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    
    figname = sprintf('%s_04%s_%s',seshInfo.fov,'aMax',seshInfo.name);
    savefig(figAccel, fullfile(figFold,[figname,'.fig']));
    saveas(figAccel,fullfile(figFold,[figname,'.png']));
    print(figAccel,fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    
    figname = sprintf('%s_05%s_%s',seshInfo.fov,'Rcorr',seshInfo.name);
    savefig(figRcorr, fullfile(figFold,[figname,'.fig']));
    saveas(figRcorr,fullfile(figFold,[figname,'.png']));
    print(figRcorr,fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    
    figname = sprintf('%s_05b%s_%s',seshInfo.fov,'Rcorrb',seshInfo.name);
    savefig(figRcorrB, fullfile(figFold,[figname,'.fig']));
    saveas(figRcorrB,fullfile(figFold,[figname,'.png']));
    print(figRcorrB,fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    
    figname = sprintf('%s_06%s_%s',seshInfo.fov,'spkMax',seshInfo.name);
    savefig(figActivity, fullfile(figFold,[figname,'.fig']));
    saveas(figActivity,fullfile(figFold,[figname,'.png']));
    print(figActivity,fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    
    figname = sprintf('%s_07%s_%s',seshInfo.fov,'syncMax',seshInfo.name);
    savefig(figSynchrony, fullfile(figFold,[figname,'.fig']));
    saveas(figSynchrony,fullfile(figFold,[figname,'.png']));
    print(figSynchrony,fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');    

end