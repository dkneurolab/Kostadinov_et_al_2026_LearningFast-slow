function seshOut = fig3_rewAnalysis(sesh, seshIn, imfs, vrfs, figFold, plotbool, savebool)
%% This parametrizes movement in each session by:
%     1. Last trial outcome

%% Set out variables that will be used to organize data
nTrials = numel(sesh.iTrial);
spkRew = smoothdata(sesh.spkRew(:,sesh.iPCrew,:),1,'gaussian',3);
spkRewBin = smoothdata(single(logical(sesh.spkRew(:,sesh.iPCrew,:))),1,'gaussian',3);
seshNames = fieldnames(sesh);
seshInfo = rmfield(sesh,seshNames(4:end));

% 1. Correct on this trial split by last trial's outcome
iUCOlast = [false(1,3); sesh.iUCO(1:end-1,:)];
iUCOlast(~sesh.iUCO(:,2),:) = false;

coeff = pca(mean(spkRew(end/2-imfs/2+1:end/2+imfs/2,:,:),3));
[~, iSort] = sort(coeff(:,1),'descend');
spkRew = spkRew(:,iSort,:);
spkRewBin = spkRewBin(:,iSort,:);

seshOut = seshIn;
seshOut.spkRewTrial = squeeze(imfs*max(mean(spkRew(end/2+1:end/2+imfs/2,:,:),2)));
seshOut.spkRewTrace = imfs*mean(mean(spkRew(:,:,sesh.iUCO(:,2)),3),2);
%% 1. Previous trial's outcome
spkMat = zeros(size(spkRew,1),size(spkRew,2),3);
xMat = cell(3,1);
vMat = xMat;
colMat = [158 31 99; 0 0 0; 194 117 20]/255;
for iOutcome = 1:3
    spkMat(:,:,iOutcome) = mean(spkRew(:,:,iUCOlast(:,iOutcome)),3);
    xMat{iOutcome,1} = sesh.xRew(:,iUCOlast(:,iOutcome));
    vMat{iOutcome,1} = sesh.vRew(:,iUCOlast(:,iOutcome));
%     seshOut.lastOutcomeSpk(iOutcome,1) = imfs*mean(mean(spkMat(end/2+1:end/2+imfs/5,:,iOutcome),2));
%     seshOut.lastOutcomeSpkCell(:,iOutcome) = imfs*mean(spkMat(end/2+1:end/2+imfs/5,:,iOutcome));
    seshOut.lastOutcomeSpk(iOutcome,1) = imfs*max(mean(spkMat(end/2+1:end/2+imfs/2,:,iOutcome),2));
    seshOut.lastOutcomeSpkCell(:,iOutcome) = imfs*max(spkMat(end/2+1:end/2+imfs/2,:,iOutcome));
end
spkMatInc = cat(2,spkMat(:,:,1),spkMat(:,:,3));
seshOut.lastOutcomeSpk(iOutcome+1,1) = imfs*max(mean(spkMatInc(end/2+1:end/2+imfs/2,:),2));
seshInfo.splitType = 'Outcome';
seshInfo.range = {'Undershoot';'Correct';'Overshoot'};
if plotbool; figOutcome = fig3_rewAnalysis_plot(imfs*spkMat,xMat,vMat,imfs,vrfs,seshInfo,colMat); end

%% 2. Previous trial's outcome
spkMat = zeros(size(spkRew,1),size(spkRew,2),3);
xMat = cell(3,1);
vMat = xMat;
colMat = [158 31 99; 0 0 0; 194 117 20]/255;
for iOutcome = 1:3
    spkMat(:,:,iOutcome) = mean(spkRewBin(:,:,iUCOlast(:,iOutcome)),3);    
    xMat{iOutcome,1} = sesh.xRew(:,iUCOlast(:,iOutcome));
    vMat{iOutcome,1} = sesh.vRew(:,iUCOlast(:,iOutcome));
%     seshOut.lastOutcomeBinSpk(iOutcome,1) = imfs*mean(mean(spkMat(end/2+1:end/2+imfs/5,:,iOutcome),2));
%     seshOut.lastOutcomeBinSpkCell(:,iOutcome) = imfs*mean(spkMat(end/2+1:end/2+imfs/5,:,iOutcome));
    seshOut.lastOutcomeBinSpk(iOutcome,1) = imfs*max(mean(spkMat(end/2+1:end/2+imfs/2,:,iOutcome),2));
    seshOut.lastOutcomeBinSpkCell(:,iOutcome) = imfs*max(spkMat(end/2+1:end/2+imfs/2,:,iOutcome));
end
spkMatInc = cat(2,spkMat(:,:,1),spkMat(:,:,3));
seshOut.lastOutcomeBinSpk(iOutcome+1,1) = imfs*max(mean(spkMatInc(end/2+1:end/2+imfs/2,:),2));
seshInfo.splitType = 'OutcomeBin';
seshInfo.range = {'Undershoot';'Correct';'Overshoot'};
if plotbool; figOutcomeBin = fig3_rewAnalysis_plot(imfs*spkMat,xMat,vMat,imfs,vrfs,seshInfo,colMat); end

%% Save maybe
if savebool && plotbool
    figname = sprintf('%s_01%s_%s',seshInfo.fov,'UCOlast',seshInfo.name);
    savefig(figOutcome, fullfile(figFold,[figname,'.fig']));
    saveas(figOutcome,fullfile(figFold,[figname,'.png']));
    print(figOutcome,fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    
    figname = sprintf('%s_02%s_%s',seshInfo.fov,'UCOlastBin',seshInfo.name);
    savefig(figOutcomeBin, fullfile(figFold,[figname,'.fig']));
    saveas(figOutcomeBin,fullfile(figFold,[figname,'.png']));
    print(figOutcomeBin,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');        
end