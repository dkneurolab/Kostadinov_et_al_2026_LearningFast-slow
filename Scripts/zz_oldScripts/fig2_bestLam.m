function [EVbox_bestmu, bestBoxLi] = fig2_bestLam(GLMdata, DMdataBoxSmall, setFiles, figFold, figFold2, savebool)


xvalFold = setFiles.inits.xvalfold;
roiNames = setFiles.mzones;
lambda = setFiles.inits.lambda;

%% Get best lambdas here

bestBeta_box_xval = zeros(size(GLMdata(1).DMtestbox,2)+1,numel(roiNames),xvalFold);
EVbox = zeros(numel(roiNames),numel(lambda),xvalFold);
for iFold = 1:xvalFold
    [~, ~, ~, varExp] = runRidgeRegression_dk2(GLMdata(iFold),lambda,'box');
    EVbox(:,:,iFold) = varExp;
end
EVbox_Lmu = mean(EVbox,3);
[~, bestBoxLi] = max(EVbox_Lmu, [], 2);

for iFold = 1:xvalFold
    bestBeta_box_xval(:,:,iFold) = runRidgeRegression_dk2best(GLMdata(iFold),lambda,'box',bestBoxLi);
end

EVbox_best = zeros(numel(roiNames),xvalFold);
for j = 1:numel(roiNames)
    EVbox_best(j,:) = EVbox(j,bestBoxLi(j),:);
end
EVbox_bestmu = mean(EVbox_best,2);

clear GLMdata trialBasis

%% Plot regressor weights
fw_box = figure;
nBases = DMdataBoxSmall.nBases;
cumnBases = [1; 1+cumsum(nBases)];
basenames = cat(1,{'Offset'},DMdataBoxSmall.fields);
if strcmpi(setFiles.inits.cellgrping,'zones')
    ijump = 1;
    plotTitle = 'Microzone ';
elseif strcmpi(setFiles.inits.cellgrping,'cells')
    ijump = floor(numel(roiNames)/10);
    nSubplot = ceil(numel(roiNames)/ijump);
    plotTitle = 'PC num ';
end
iSubplot = 0;
for iPred = 1:ijump:numel(roiNames)
    iSubplot = iSubplot + 1;
    wboxtemp = squeeze(bestBeta_box_xval(:,iPred,:));
    xwbox = (0.5:1:size(wboxtemp,1))';
    wboxmean = nanmean(wboxtemp,2);
    wboxsd = nanstd(wboxtemp,[],2);
    wboxylim = max(sum([abs(wboxmean) abs(wboxsd)],2))*1.1;
    if wboxylim < .11; wboxylim = .11; end
    if isnan(wboxylim); wboxylim = 0.5; end
    if strcmpi(setFiles.inits.cellgrping,'zones')
        subplot(numel(roiNames),1,iSubplot);
    elseif strcmpi(setFiles.inits.cellgrping,'cells')
        subplot(ceil(nSubplot/2),2,iSubplot);
    end
    bar(xwbox,wboxmean,'EdgeColor','none'); hold on;
    errorbar(xwbox,wboxmean,wboxsd,'k','CapSize',0,'Marker','none','LineStyle','none');
    ylim([-wboxylim wboxylim]);
    xlim([0 numel(wboxmean)]);
    box off
    plot([cumnBases'; cumnBases'],[-wboxylim*ones(1,numel(cumnBases)); wboxylim*ones(1,numel(cumnBases))],'k-')
    xticks(cumnBases);
    xticklabels({});
    title(sprintf('%s%s',plotTitle,roiNames{iPred}(2:end)));
    if iSubplot == nSubplot
        %             suptitle('Summary of microzone beta values, boxcar filters')
        namelocs = mean(cat(2,[0; cumnBases(1:end-1)], cumnBases),2);
        text(namelocs,repmat(-1.5*wboxylim,size(namelocs)),basenames,'HorizontalAlignment','center');
    end
end

%% Saveout
if savebool
    savefig(fw_box,fullfile(figFold,'Betas_box.fig'),'compact');
    print(fw_box,fullfile(figFold,'Betas_box.eps'), '-depsc', '-painters');
    saveas(fw_box,fullfile(figFold,'Betas_box.png'));
    savefig(fw_box,fullfile(figFold2,'Betas_box.fig'),'compact');
    print(fw_box,fullfile(figFold2,'Betas_box.eps'), '-depsc', '-painters');
    saveas(fw_box,fullfile(figFold2,'Betas_box.png'));
    save(fullfile(figFold,sprintf('GLMoptims_%s.mat',setFiles.inits.cellgrping)), 'bestBeta_box_xval', 'bestBoxLi', 'lambda', 'EVbox', 'EVbox_best', 'EVbox_Lmu','EVbox_bestmu');
end
