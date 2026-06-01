function statStruct = fig2_plotHoldout(trialBasis,roiNames,bestLam,GLMholdout,DMdataSmall,inits,figFold,figFold2,savebool)
% Only take good trials

nShuf = 200;

fitcolor = distinguishable_colors(2e3,[1 1 1; 0 0 0]);
fitcolor(1:6,:) = fitcolor([6,5,2,3,1,4],:);

iTrials = GLMholdout.iTrials.iTrials;
iHoldTrials = GLMholdout.iTrials.iHoldTrials;
% iRews = GLMholdout.iTrials.iRews;
% iHoldRews = GLMholdout.iTrials.iHoldRews;
% iGoodTrials = GLMholdout.iTrials.trialsGOOD;
% trialbasis.expt = trialbasis0.expt(iGoodTrials);
% trialbasis.bs = trialbasis0.bs(iGoodTrials);
% trialbasis.bsbox = trialbasis0.bsbox(iGoodTrials);
% nBases = round(inits.basisdur/inits.basiswid);
plotwid = inits.basisdur/1000*inits.imfs;

% if ~iGoodTrials(1)
%     trialBasis.bs(1).CueBeep = zeros(size(trialBasis.bs(1).TrialON));
%     trialBasis.bs(1).CueReward = zeros(size(trialBasis.bs(1).TrialON));
%     trialBasis.bs(1).CueReward = trialBasis.bs(1).CueReward(:,1:end/2+1);
%     trialBasis.bs(1).RandReward = zeros(size(trialBasis.bs(1).TrialON));
%     trialBasis.bsbox(1).CueBeep = zeros(size(trialBasis.bsbox(1).TrialON));
%     trialBasis.bsbox(1).CueReward = zeros(size(trialBasis.bsbox(1).TrialON));
%     trialBasis.bsbox(1).CueReward = trialBasis.bsbox(1).CueReward(:,1:end/2);
%     trialBasis.bsbox(1).RandReward = zeros(size(trialBasis.bsbox(1).TrialON));
% end

%% Make design and data matrices
Xtrain = makeDesignMatrix_dk(trialBasis, iTrials, 'bsbox');
% DMtrainTrial = DMtrainTrial(:,DMdataSmall.iSmall);
Xtest = makeDesignMatrix_dk(trialBasis, iHoldTrials, 'bsbox');
% DMtestTrial = DMtestTrial(:,DMdataSmall.iSmall);
% DMtrainRew = makeDesignMatrix_dk(rewbasis, iRews, 'bsbox');
% DMtrainRew = DMtrainRew(:,DMdataSmall.iSmall);
% DMtestRew = makeDesignMatrix_dk(rewbasis, iHoldRews, 'bsbox');
% DMtestRew = DMtestRew(:,DMdataSmall.iSmall);

% Xtrain = [DMtrainTrial; DMtrainRew]; Xtrain = [ones(size(Xtrain,1),1) Xtrain];
% Xtest = [DMtestTrial; DMtestRew]; Xtest = [ones(size(Xtest,1),1) Xtest];

Xtrain = [ones(size(Xtrain,1),1) Xtrain];
Xtest = [ones(size(Xtest,1),1) Xtest];

Ytrain = zeros(size(Xtrain,1),numel(roiNames));
Ytest = zeros(size(Xtest,1),numel(roiNames));
bestBeta = zeros(size(Xtrain,2),numel(roiNames));

if ~exist(fullfile(figFold,sprintf('bestBetaShuf_%s.mat',inits.cellgrping)),'file')

    XTrial_shuf = fig2_makeDesignMatrix_dk_shuffle(trialBasis, iTrials, 'bsbox',nShuf,true(numel(fieldnames(trialBasis.bsbox)),1));
    bestBetaShuf = zeros(size(Xtrain,2),numel(roiNames),nShuf);
    tic;
    for j = 1:numel(roiNames)

        YtrialTemp = getSpikes_dk(trialBasis,iTrials,roiNames{j});
    %     YrewTemp = getSpikes_dk(rewbasis,iRews,roiNames{j});
        Ytrain(:,j) = YtrialTemp*inits.imfs;

        YtrialTemp = getSpikes_dk(trialBasis,iHoldTrials,roiNames{j});
    %     YrewTemp = getSpikes_dk(rewbasis,iHoldRews,roiNames{j});
        Ytest(:,j) = YtrialTemp*inits.imfs;
        bestBeta(:,j) = (Xtrain'*Xtrain+bestLam(j)*eye(size(Xtrain,2))) \ (Xtrain' * Ytrain(:,j));
        for iShuf = 1:nShuf
            Xtrain_shuf = XTrial_shuf{iShuf};
            Xtrain_shuf = [ones(size(Xtrain_shuf,1),1) Xtrain_shuf]; %#ok<AGROW>
            bestBetaShuf(:,j,iShuf) = (Xtrain_shuf'*Xtrain_shuf+bestLam(j)*eye(size(Xtrain_shuf,2))) \ (Xtrain_shuf' * Ytrain(:,j));
        end
        if ~mod(j,10)
            toc;
        end
    end

    clear XTrial_shuf Xtrain_shuf
    save(fullfile(figFold,sprintf('bestBetaShuf_%s.mat',inits.cellgrping)),'bestBetaShuf');

else
    bestBetaShuf = load(fullfile(figFold,sprintf('bestBetaShuf_%s.mat',inits.cellgrping)));
    bestBetaShuf = bestBetaShuf.bestBetaShuf;
    nShuf = size(bestBetaShuf,3);
    for j = 1:numel(roiNames)
        YtrialTemp = getSpikes_dk(trialBasis,iTrials,roiNames{j});
        Ytrain(:,j) = YtrialTemp*inits.imfs;
        YtrialTemp = getSpikes_dk(trialBasis,iHoldTrials,roiNames{j});
        Ytest(:,j) = YtrialTemp*inits.imfs;
        bestBeta(:,j) = (Xtrain'*Xtrain+bestLam(j)*eye(size(Xtrain,2))) \ (Xtrain' * Ytrain(:,j));
    end
end

x = [-10:10]';
y = exp(-x.^2/(2*1^2));
y(1:10) = 0; y = y./sum(y);

YtestPlot = convn(Ytest,y,'same');
YtestHat = Xtest*bestBeta;
YtestHatPlot = convn(YtestHat,y,'same');

YtestHat_shuf = zeros(size(Xtest,1),numel(roiNames),nShuf);

for iShuf = 1:nShuf
    YtestHat_shuf(:,:,iShuf) = Xtest*bestBetaShuf(:,:,iShuf);
end

EV = zeros(numel(roiNames),1); Rsq = EV; pRsq = EV;
EV_shuf = zeros(numel(roiNames),nShuf); Rsq_shuf = EV_shuf; pRsq_shuf = Rsq_shuf;
for iCell = 1:numel(roiNames)
    EV(iCell) = calCrossValExpVar(Ytrain(:,iCell), Ytest(:,iCell), YtestHat(:,iCell));
    [Rsq0, Pval0] = corrcoef(Ytest(:,iCell),YtestHat(:,iCell));
    Rsq(iCell) = Rsq0(2,1);
    pRsq(iCell) = Pval0(2,1);
    for iShuf = 1:nShuf
        EV_shuf(iCell,iShuf) = calCrossValExpVar(Ytrain(:,iCell), Ytest(:,iCell), YtestHat_shuf(:,iCell,iShuf));
        [Rsq0, Pval0] = corrcoef(Ytest(:,iCell),YtestHat_shuf(:,iCell,iShuf));
        Rsq_shuf(iCell,iShuf) = Rsq0(2,1);
        pRsq_shuf(iCell,iShuf) = Pval0(2,1);
    end
end

EV_sig = EV - (mean(EV_shuf,2)+2*std(EV_shuf,[],2)) > 0;
Rsq_sig = Rsq - (mean(Rsq_shuf,2)+2*std(Rsq_shuf,[],2)) > 0;

YtestHat_shufMu = mean(YtestHat_shuf,3);

clear YtrialTemp YrewTemp Rsq0 Pval0 YtestHat_shuf bestBetaShuf
%% Make figure!
statStruct = v2struct(EV,Rsq,pRsq, EV_shuf, Rsq_shuf, pRsq_shuf, EV_sig, Rsq_sig);
% meanstruct = struct;
% meanstruct.EV = EV;
% meanstruct.Rsq = Rsq;
% meanstruct.pRsq = pRsq;
% meanstruct.EV_shuf = EV_shuf;
% meanstruct.Rsq_shuf = Rsq_shuf;
% meanstruct.pRsq_shuf = pRsq_shuf;
meanStruct.bestBeta = bestBeta;
meanStruct.Ytest = Ytest;
meanStruct.YtestHat = YtestHat;
testTrialB.expt = trialBasis.expt(iHoldTrials);
% testTrialB.bs = trialBasis.bs(iHoldTrials);
testTrialB.bsbox = trialBasis.bsbox(iHoldTrials);
% testRewB.expt = rewbasis.expt(iHoldRews);
% testRewB.bs = rewbasis.bs(iHoldRews);
% testRewB.bsbox = rewbasis.bsbox(iHoldRews);
% nTrials = numel(iHoldTrials);
maxduri = (inits.prebuffsec+inits.postbuffsec+inits.resptime)*inits.imfs;
trialends = [0; cumsum(cat(1,testTrialB.expt.duration))];
% rewends = [0; cumsum(cat(1,testRewB.expt.duration))];
nTrials = numel(testTrialB.expt);
% nRewards = numel(testRewB.expt);

%% Compute and plot
% Trial onsets - note that this relies on trialbasis and rewbasis to still
% have their full number of predictors
trialONall = cat(1,testTrialB.bsbox.TrialON);
trialONi = threshdetect_dk(trialONall(:,end/2+1),.5,1);
tplot = (-.1:1/inits.imfs:.9-1/inits.imfs)+1/inits.imfs/2;
tplot2 = (-.5:1/inits.imfs:.5-1/inits.imfs)+1/inits.imfs/2;
yOn = zeros(30,numel(roiNames),numel(trialONi)); yOnHat = yOn; yOnPlot = yOn; yOnHatPlot = yOn; yOnHatShuf = yOn;
for i = 1:numel(trialONi)
    yinds = trialONi(i)-3:trialONi(i)+27-1;
    yOn(:,:,i) = Ytest(yinds,:);
    yOnHat(:,:,i) = YtestHat(yinds,:);
    yOnPlot(:,:,i) = YtestPlot(yinds,:);
    yOnHatPlot(:,:,i) = YtestHatPlot(yinds,:);
    yOnHatShuf(:,:,i) = YtestHat_shufMu(yinds,:);
end

% Reward time and movement - can be done together
yREWc = zeros(plotwid,numel(roiNames),nTrials)*nan; yREWcHat = yREWc; yREWcPlot = yREWc; yREWcHatPlot = yREWc; yREWcHatShuf = yREWc;
movi = zeros(numel(cat(1,testTrialB.expt.movonset)),1); movdummy = 0;
for i = 1:nTrials        
    if any(testTrialB.bsbox(i).RewTcorrect(:) > 0) && size(testTrialB.bsbox(i).RewTcorrect,1) < maxduri
        yi = threshdetect_dk(testTrialB.bsbox(i).RewTcorrect(:,1),.5,1)+trialends(i);
        yinds = yi-plotwid/2:yi+plotwid/2-1;
        yREWc(:,:,i) = Ytest(yinds,:);
        yREWcPlot(:,:,i) = YtestPlot(yinds,:);
        yREWcHat(:,:,i) = YtestHat(yinds,:);
        yREWcHatPlot(:,:,i) = YtestHatPlot(yinds,:);
        yREWcHatShuf(:,:,i) = YtestHat_shufMu(yinds,:);
    end     
    if ~isempty(testTrialB.expt(i).movonset)
        movtemp = testTrialB.expt(i).movonset+1;            
        movi(movdummy+1:movdummy+numel(movtemp)) = movtemp+1-testTrialB.expt(i).start+trialends(i);
        movdummy = movdummy+numel(movtemp);
    end
end    
yREW_good = ~isnan(squeeze(yREWc(1,1,:)));
% yREWc = yREWc(:,:,yREW_good(:,1)); yREWcPlot = yREWcPlot(:,:,yREW_good(:,2));
% yREWcHat = yREWcHat(:,:,yREW_good(:,2)); yREWcHatPlot = yREWcHatPlot(:,:,yREW_good(:,2));
% nUCO = sum(yREW_good);

yMov = zeros(plotwid,numel(roiNames),numel(movi)); yMovHat = yMov; yMovPlot = yMov; yMovHatPlot = yMov; yMovHatShuf = yMov;
for j = 1:numel(movi)
    yinds = movi(j)-plotwid/2:movi(j)+plotwid/2-1;
    yMov(:,:,j) = Ytest(yinds,:);
    yMovPlot(:,:,j) = YtestPlot(yinds,:);
    yMovHat(:,:,j) = YtestHat(yinds,:);
    yMovHatPlot(:,:,j) = YtestHatPlot(yinds,:);
    yMovHatShuf(:,:,j) = YtestHat_shufMu(yinds,:);
end

yOnPlotMean = mean(yOnPlot,3);
yOnPlotSEM = std(yOnPlot,[],3)/sqrt(numel(trialONi));

yMovPlotMean = mean(yMovPlot,3);
yMovPlotSEM = std(yMovPlot,[],3)/sqrt(numel(movi));

yREWcPlotMean = nanmean(yREWcPlot,3);
yREWcPlotSEM = nanstd(yREWcPlot,[],3)/sqrt(numel(yREW_good));


j = 0;
for nPlots = 1:5:numel(roiNames)
    f1 = figure;
    for nRois = 1:5
        j = j+1;
        if j <= numel(roiNames)
        nRow = mod(nRois,5); if nRow == 0; nRow = 5; end
        iPlot = 4*(nRow-1)+1;
        subplot(5,4,iPlot); hold on;
        h = boundedline(tplot,yOnPlotMean(:,j),yOnPlotSEM(:,j),'cmap',[0 0 0],'alpha');
        set(h,'linewidth',2);
        plot(tplot,mean(yOnHatPlot(:,j,:),3),'Color',fitcolor(j,:),'LineWidth',2.5);
        plot(tplot,mean(yOnHatShuf(:,j,:),3),'--','Color',fitcolor(j,:)/2,'LineWidth',2.5);
        xline(0,'--','LineWidth',1);
        plotmax = ceil(1.1*(max(yOnPlotMean(:,j)+yOnPlotSEM(:,j))));
        if plotmax < 5; plotmax = 5; end; if plotmax > 50; plotmax = 50; end; if isnan(plotmax); plotmax = 5; end
        axis([-.1 0.9 0 plotmax]);
        xticks([-1:.25:1]); if nRow < 5; xticklabels({}); end
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',18,'TickDir','out');
        if nRow == 3; ylabel('Firing rate (spk/s)'); end
        if nRow == 1
            title(sprintf('Trial reward, PCs %i:%i',nPlots,nPlots+4));
            if j == numel(roiNames); title(sprintf('Trial reward, PCs %i:%i',nPlots,numel(roiNames))); end
        end
        
        subplot(5,4,iPlot+1); hold on;
        h = boundedline(tplot2,yMovPlotMean(:,j),yMovPlotSEM(:,j),'cmap',[0 0 0],'alpha');
        set(h,'linewidth',2);
        plot(tplot2,mean(yMovHatPlot(:,j,:),3),'Color',fitcolor(j,:),'LineWidth',2.5);
        plot(tplot2,mean(yMovHatShuf(:,j,:),3),'--','Color',fitcolor(j,:)/2,'LineWidth',2.5);
        xline(0,'--','LineWidth',1);
        plotmax = ceil(1.1*(max(yMovPlotMean(:,j)+yMovPlotSEM(:,j))));
        if plotmax < 5; plotmax = 5; end; if plotmax > 50; plotmax = 50; end; if isnan(plotmax); plotmax = 5; end
        axis([-.5 0.5 0 plotmax]);
        xticks([-1:.25:1]); if nRow < 5; xticklabels({}); else; xlabel('Time (s)'); end
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',18,'TickDir','out');
        if nRow == 1
            title(sprintf('Trial reward, PCs %i:%i',nPlots,nPlots+4));
            if j == numel(roiNames); title(sprintf('Trial reward, PCs %i:%i',nPlots,numel(roiNames))); end
        end
        
        subplot(5,4,iPlot+2); hold on;
        h = boundedline(tplot2,yREWcPlotMean(:,j),yREWcPlotSEM(:,j),'cmap',[0 0 0],'alpha');
        set(h,'linewidth',2);
        plot(tplot2,mean(yREWcHatPlot(:,j,:),3),'Color',fitcolor(j,:),'LineWidth',2.5);
        plot(tplot2,mean(yREWcHatShuf(:,j,:),3),'--','Color',fitcolor(j,:)/2,'LineWidth',2.5);
        xline(0,'--','LineWidth',1);
        plotmax = ceil(1.1*(max(yREWcPlotMean(:,j)+yREWcPlotSEM(:,j))));
        if plotmax < 5; plotmax = 5; end; if plotmax > 50; plotmax = 50; end;  if isnan(plotmax); plotmax = 5; end
        axis([-.5 0.5 0 plotmax]);
        xticks([-1:.25:1]); if nRow < 5; xticklabels({}); end
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',18,'TickDir','out');
        if nRow == 1
            title(sprintf('Trial reward, PCs %i:%i',nPlots,nPlots+4));
            if j == numel(roiNames); title(sprintf('Trial reward, PCs %i:%i',nPlots,numel(roiNames))); end
        end
        
        
        subplot(5,4,iPlot+3);
        statsummary{1,1} = sprintf('Exp. var. = %.1f',100*EV(j));
        statsummary{3,1} = sprintf('Rsq. = %.2f',Rsq(j));
        if EV_sig(j)
            statsummary{2,1} = sprintf('Exp. var. p < 0.05');
        else
            statsummary{2,1} = sprintf('Exp. var. p > 0.05');
        end
        if Rsq_sig(j)
            statsummary{4,1} = sprintf('Rsq. p < 0.05');
        else
            statsummary{4,1} = sprintf('Rsq. p > 0.05');
        end
        
        text(0.5,0.5,statsummary,'FontSize',18,'HorizontalAlignment','center'); axis off;
        
        end
    end
    
    if savebool
        figname = sprintf('00_FullModel_PCs_%03i_%03i',nPlots,nPlots+4);
        figFold_fig = fullfile(figFold,'figsMat'); if ~exist(figFold_fig,'dir'); mkdir(figFold_fig); end
        figFold_eps = fullfile(figFold,'figsEps'); if ~exist(figFold_eps,'dir'); mkdir(figFold_eps); end
        figFold_png = fullfile(figFold,'figsPng'); if ~exist(figFold_png,'dir'); mkdir(figFold_png); end
        figFold2_fig = fullfile(figFold2,'figsMat'); if ~exist(figFold2_fig,'dir'); mkdir(figFold2_fig); end
        figFold2_eps = fullfile(figFold2,'figsEps'); if ~exist(figFold2_eps,'dir'); mkdir(figFold2_eps); end
        figFold2_png = fullfile(figFold2,'figsPng'); if ~exist(figFold2_png,'dir'); mkdir(figFold2_png); end
        savefig(f1, fullfile(figFold_fig,[figname,'.fig']));
        saveas(f1,fullfile(figFold_png,[figname,'.png']));
        print(f1,fullfile(figFold_eps,[figname,'.eps']), '-depsc', '-painters');
        savefig(f1, fullfile(figFold2_fig,[figname,'.fig']));
        saveas(f1,fullfile(figFold2_eps,[figname,'.png']));
        print(f1,fullfile(figFold2_png,[figname,'.eps']), '-depsc', '-painters');
    end
end

% plot(tplot,mean(yOnPlot,3),'k','LineWidth',2); hold on;
% plot(tplot,mean(yOnHatPlot,3),'LineWidth',1.5); title('Trial ONSET');
meanStruct.yOn = yOn;
meanStruct.yOnHat = yOnHat;
meanStruct.yOnPlot = yOnPlot;
meanStruct.yOnHatPlot = yOnHatPlot;
meanStruct.yOnHatShuf = yOnHatShuf;
meanStruct.yMov = yMov;
meanStruct.yMovPlot = yMovPlot;
meanStruct.yMovHat = yMovHat;
meanStruct.yMovHatPlot = yMovHatPlot;
meanStruct.yMovHatShuf = yMovHatShuf;
meanStruct.yREWc = yREWc;
meanStruct.yREWcPlot = yREWcPlot;
meanStruct.yREWcHat = yREWcHat;
meanStruct.yREWcHatPlot = yREWcHatPlot;
meanStruct.yREWcHatShuf = yREWcHatShuf;

clear yi yinds licktemp lickdummy licki movtemp movdummy movi yREWu yREWu_hat yREWc yREWc_hat yREWo yREWo_hat yLick yLick_hat yMov yMov_hat 


if savebool
   save(fullfile(figFold,'plotStructs.mat'),'meanStruct','statStruct'); 
end


end