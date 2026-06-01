function fig2_GLMdistDep(sessions,dataFold,paperFold,savebool)
%%
figFold = fullfile(paperFold,'Figures','Fig2','Sub4_DistDep');
if ~exist(figFold,'dir'); mkdir(figFold); end
xDistMax = 500;
xDistBin = 20;

%% Loop through and extract data
distSummary = struct;
for iSesh = 1:numel(sessions)
    dataPath = fullfile(dataFold,sprintf('Version_%i',sessions(iSesh).version),sessions(iSesh).name,...
        sprintf('%s_%s',sessions(iSesh).date,sessions(iSesh).fov),'GLM','glmRfiles');
    glmLocal = load(fullfile(dataPath,'GLMsig+VisMovRewLick.mat')); glmLocal = glmLocal.GLMoutput;
    fovParams = load(fullfile(dataFold,sprintf('Version_%i',sessions(iSesh).version),sessions(iSesh).name,...
        sprintf('%s_%s',sessions(iSesh).date,sessions(iSesh).fov),'fovparams.mat'));
    fovSets = load(fullfile(dataFold,sprintf('Version_%i',sessions(iSesh).version),sessions(iSesh).name,...
        sprintf('%s_%s',sessions(iSesh).date,sessions(iSesh).fov),'setsonly.mat'));
    fovParams = fovParams.fovparams.final; fovSets = fovSets.sets;
    
    % Find cells with any model at all and trim other arrays to size
    iCells = false(size(glmLocal))';
    for j = 1:numel(glmLocal)
        if ~isempty(glmLocal(j).GLMrankPCA) && glmLocal(j).sigBool
            iCells(j) = true;
        end
    end
    fovParams = fovParams(iCells);
    fovSets.Fspont = fovSets.Fspont(:,iCells);
    glmLocal = glmLocal(iCells);
    
    % Correct for DK105 0223
    if strcmpi(sessions(iSesh).name,'dk105') && strcmpi(sessions(iSesh).date,'2018-02-23') && ...
            strcmpi(sessions(iSesh).fov,'Xnor')
        for ifix = 1:numel(fovParams)
            fovParams(ifix).rois = fliplr(fovParams(ifix).rois);
        end
    end
    
    % Calculate intercellular distance in x and y    
    [xDists,yDists] = fig2_GLMdistCalc(fovParams);
    xDistum = xDists*670/512; yDistum = yDists*670/512; clear xDists yDists;
    
    % Calculate paired logic gates
    visAnd = and(cat(1,glmLocal.sigBool_shufVis_proj),cat(1,glmLocal.sigBool_shufVis_proj)');
    visXnor = and(~cat(1,glmLocal.sigBool_shufVis_proj),~cat(1,glmLocal.sigBool_shufVis_proj)');
    visXnor = or(visAnd,visXnor);
    movAnd = and(cat(1,glmLocal.sigBool_shufMov_proj),cat(1,glmLocal.sigBool_shufMov_proj)');
    movXnor = and(~cat(1,glmLocal.sigBool_shufMov_proj),~cat(1,glmLocal.sigBool_shufMov_proj)');
    movXnor = or(movAnd,movXnor);    
    rewAnd = and(cat(1,glmLocal.sigBool_shufRew_proj),cat(1,glmLocal.sigBool_shufRew_proj)');
    rewXnor = and(~cat(1,glmLocal.sigBool_shufRew_proj),~cat(1,glmLocal.sigBool_shufRew_proj)');
    rewXnor = or(rewAnd,rewXnor);
    lickAnd = and(cat(1,glmLocal.sigBool_shufLick_proj),cat(1,glmLocal.sigBool_shufLick_proj)');
    lickXnor = and(~cat(1,glmLocal.sigBool_shufLick_proj),~cat(1,glmLocal.sigBool_shufLick_proj)');
    lickXnor = or(lickAnd,lickXnor);
    
    % Calculate correlations and put everyone into upper triangle
    Rspont = corrcoef(fovSets.Fspont);
    iTriu = triu(Rspont,1); iTriu = logical(iTriu(:));
    Rspont2 = Rspont(iTriu);     xDistum2 = xDistum(iTriu); yDistum2 = yDistum(iTriu);     
    visAnd = visAnd(iTriu); visXnor = visXnor(iTriu);
    movAnd = movAnd(iTriu); movXnor = movXnor(iTriu);
    rewAnd = rewAnd(iTriu); rewXnor = rewXnor(iTriu);
    lickAnd = lickAnd(iTriu); lickXnor = lickXnor(iTriu);
     
    [hCounts,hEdges,hBinid] = histcounts(xDistum2,[0:xDistBin:xDistMax]);
    hCounts = hCounts'; hEdges = hEdges';
    
    RbinMu = zeros(numel(hCounts,1));
    RbinSD = RbinMu;
    visAndMu = RbinMu; visAndSD = RbinMu;
    visXnorMu = RbinMu; visXnorSD = RbinMu;
    movAndMu = RbinMu; movAndSD = RbinMu;
    movXnorMu = RbinMu; movXnorSD = RbinMu;
    rewAndMu = RbinMu; rewAndSD = RbinMu;
    rewXnorMu = RbinMu; rewXnorSD = RbinMu;
    lickAndMu = RbinMu; lickAndSD = RbinMu;
    lickXnorMu = RbinMu; lickXnorSD = RbinMu;
    for iBin = 1:numel(hCounts)
        RbinMu(iBin,1) = mean(Rspont2(hBinid == iBin));
        RbinSD(iBin,1) = std(Rspont2(hBinid == iBin));
        
        visAndMu(iBin,1) = mean(visAnd(hBinid == iBin));
        visAndSD(iBin,1) = std(visAnd(hBinid == iBin));
        visXnorMu(iBin,1) = mean(visXnor(hBinid == iBin));
        visXnorSD(iBin,1) = std(visXnor(hBinid == iBin));
        
        movAndMu(iBin,1) = mean(movAnd(hBinid == iBin));
        movAndSD(iBin,1) = std(movAnd(hBinid == iBin));
        movXnorMu(iBin,1) = mean(movXnor(hBinid == iBin));
        movXnorSD(iBin,1) = std(movXnor(hBinid == iBin));
        
        rewAndMu(iBin,1) = mean(rewAnd(hBinid == iBin));
        rewAndSD(iBin,1) = std(rewAnd(hBinid == iBin));
        rewXnorMu(iBin,1) = mean(rewXnor(hBinid == iBin));
        rewXnorSD(iBin,1) = std(rewXnor(hBinid == iBin));
        
        lickAndMu(iBin,1) = mean(lickAnd(hBinid == iBin));
        lickAndSD(iBin,1) = std(lickAnd(hBinid == iBin));
        lickXnorMu(iBin,1) = mean(lickXnor(hBinid == iBin));
        lickXnorSD(iBin,1) = std(lickXnor(hBinid == iBin));
    end
    
    distSummary(iSesh).name = sessions(iSesh).name;
    distSummary(iSesh).date = sessions(iSesh).date;
    distSummary(iSesh).fov = sessions(iSesh).fov;
    distSummary(iSesh).iLobV = strcmpi(distSummary(iSesh).fov,'lobv');
    distSummary(iSesh).iSim2 = strcmpi(distSummary(iSesh).fov,'sim2');
    distSummary(iSesh).iCells = iCells;
    
    distSummary(iSesh).distAll = v2struct(hCounts,hEdges,RbinMu,RbinSD,...
        visAndMu,visAndSD,visXnorMu,visXnorSD,...
        movAndMu,movAndSD,movXnorMu,movXnorSD,...
        rewAndMu,rewAndSD,rewXnorMu,rewXnorSD,...
        lickAndMu,lickAndSD,lickXnorMu,lickXnorSD);
    
    distSummary(iSesh).RbinMu = RbinMu;
    distSummary(iSesh).visAndMu = visAndMu;
    distSummary(iSesh).visXnorMu = visXnorMu;
    distSummary(iSesh).movAndMu = movAndMu;
    distSummary(iSesh).movXnorMu = movXnorMu;
    distSummary(iSesh).rewAndMu = rewAndMu;
    distSummary(iSesh).rewXnorMu = rewXnorMu;
    distSummary(iSesh).lickAndMu = lickAndMu;
    distSummary(iSesh).lickXnorMu = lickXnorMu;
    
    
%     figure; 
%     subplot(1,2,1); hold on;
%     xplot = hEdges(1:end-1)+xDistBin/2;
%     plot(xplot,RbinMu,'-k','LineWidth',2);
%     plot(xplot,visAndMu,'-b','LineWidth',1.5);
%     plot(xplot,movAndMu,'-g','LineWidth',1.5);
%     plot(xplot,rewAndMu,'-c','LineWidth',1.5);
%     plot(xplot,lickAndMu,'-m','LineWidth',1.5);
%     xlim([0 1]);
%     
%     subplot(1,2,2); hold on;
%     xplot = hEdges(1:end-1)+xDistBin/2;
%     plot(xplot,RbinMu,'-k','LineWidth',2);
%     plot(xplot,visXnorMu,'-b','LineWidth',1.5);
%     plot(xplot,movXnorMu,'-g','LineWidth',1.5);
%     plot(xplot,rewXnorMu,'-c','LineWidth',1.5);
%     plot(xplot,lickXnorMu,'-m','LineWidth',1.5);
%     xlim([0 1]);
%     
%     suptitle(sprintf('%s, %s, %s',sessions(iSesh).name,sessions(iSesh).date,sessions(iSesh).fov));
    
end

%% Summary plots

xplot = hEdges(1:end-1)+xDistBin/2;
iLobV = cat(1,distSummary.iLobV); distLobV = distSummary(iLobV);
iSim2 = cat(1,distSummary.iSim2); distSim2  = distSummary(iSim2);
lobVRbinMu = mean(cat(2,distLobV.RbinMu),2,'omitnan');
lobVRbinSD = std(cat(2,distLobV.RbinMu),[],2,'omitnan');
lobVvisAndMu = mean(cat(2,distLobV.visAndMu),2,'omitnan');
lobVvisAndSD = std(cat(2,distLobV.visAndMu),[],2,'omitnan');
lobVmovAndMu = mean(cat(2,distLobV.movAndMu),2,'omitnan');
lobVmovAndSD = std(cat(2,distLobV.movAndMu),[],2,'omitnan');
lobVrewAndMu = mean(cat(2,distLobV.rewAndMu),2,'omitnan');
lobVrewAndSD = std(cat(2,distLobV.rewAndMu),[],2,'omitnan');
lobVlickAndMu = mean(cat(2,distLobV.lickAndMu),2,'omitnan');
lobVlickAndSD = std(cat(2,distLobV.lickAndMu),[],2,'omitnan');

lobVvisXnorMu = mean(cat(2,distLobV.visXnorMu),2,'omitnan');
lobVvisXnorSD = std(cat(2,distLobV.visXnorMu),[],2,'omitnan');
lobVmovXnorMu = mean(cat(2,distLobV.movXnorMu),2,'omitnan');
lobVmovXnorSD = std(cat(2,distLobV.movXnorMu),[],2,'omitnan');
lobVrewXnorMu = mean(cat(2,distLobV.rewXnorMu),2,'omitnan');
lobVrewXnorSD = std(cat(2,distLobV.rewXnorMu),[],2,'omitnan');
lobVlickXnorMu = mean(cat(2,distLobV.lickXnorMu),2,'omitnan');
lobVlickXnorSD = std(cat(2,distLobV.lickXnorMu),[],2,'omitnan');

sim2RbinMu = mean(cat(2,distSim2.RbinMu),2,'omitnan');
sim2RbinSD = std(cat(2,distSim2.RbinMu),[],2,'omitnan');
sim2visAndMu = mean(cat(2,distSim2.visAndMu),2,'omitnan');
sim2visAndSD = std(cat(2,distSim2.visAndMu),[],2,'omitnan');
sim2movAndMu = mean(cat(2,distSim2.movAndMu),2,'omitnan');
sim2movAndSD = std(cat(2,distSim2.movAndMu),[],2,'omitnan');
sim2rewAndMu = mean(cat(2,distSim2.rewAndMu),2,'omitnan');
sim2rewAndSD = std(cat(2,distSim2.rewAndMu),[],2,'omitnan');
sim2lickAndMu = mean(cat(2,distSim2.lickAndMu),2,'omitnan');
sim2lickAndSD = std(cat(2,distSim2.lickAndMu),[],2,'omitnan');

sim2visXnorMu = mean(cat(2,distSim2.visXnorMu),2,'omitnan');
sim2visXnorSD = std(cat(2,distSim2.visXnorMu),[],2,'omitnan');
sim2movXnorMu = mean(cat(2,distSim2.movXnorMu),2,'omitnan');
sim2movXnorSD = std(cat(2,distSim2.movXnorMu),[],2,'omitnan');
sim2rewXnorMu = mean(cat(2,distSim2.rewXnorMu),2,'omitnan');
sim2rewXnorSD = std(cat(2,distSim2.rewXnorMu),[],2,'omitnan');
sim2lickXnorMu = mean(cat(2,distSim2.lickXnorMu),2,'omitnan');
sim2lickXnorSD = std(cat(2,distSim2.lickXnorMu),[],2,'omitnan');


fSumAnd = figure;
subplot(1,2,1);
yyaxis left
h1 = boundedline(xplot,lobVvisAndMu,lobVvisAndSD/sqrt(sum(iLobV)),'cmap',[.1 .1 .8],'alpha'); set(h1,'linewidth',1.5);
h1 = boundedline(xplot,lobVmovAndMu,lobVmovAndSD/sqrt(sum(iLobV)),'cmap',[.1 .8 .1],'alpha'); set(h1,'linewidth',1.5);
h1 = boundedline(xplot,lobVrewAndMu,lobVrewAndSD/sqrt(sum(iLobV)),'cmap',[.1 .8 .8],'alpha'); set(h1,'linewidth',1.5);
h1 = boundedline(xplot,lobVlickAndMu,lobVlickAndSD/sqrt(sum(iLobV)),'cmap',[.5 .1 .8],'alpha'); set(h1,'linewidth',1.5);
ylim([0 1]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',16,'TickDir','out');
yyaxis right
h1 = boundedline(xplot,lobVRbinMu,lobVRbinSD/sqrt(sum(iLobV)),'cmap',[0 0 0],'alpha'); set(h1,'linewidth',1.5);
ylim([0 .2]);xlim([0 300]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',16,'TickDir','out');
title('Lobule V, AND logic');
axis square;

subplot(1,2,2);
yyaxis left
h1 = boundedline(xplot,sim2visAndMu,sim2visAndSD/sqrt(sum(iLobV)),'cmap',[.1 .1 .8],'alpha'); set(h1,'linewidth',1.5);
h1 = boundedline(xplot,sim2movAndMu,sim2movAndSD/sqrt(sum(iLobV)),'cmap',[.1 .8 .1],'alpha'); set(h1,'linewidth',1.5);
h1 = boundedline(xplot,sim2rewAndMu,sim2rewAndSD/sqrt(sum(iLobV)),'cmap',[.1 .8 .8],'alpha'); set(h1,'linewidth',1.5);
h1 = boundedline(xplot,sim2lickAndMu,sim2lickAndSD/sqrt(sum(iLobV)),'cmap',[.5 .1 .8],'alpha'); set(h1,'linewidth',1.5);
ylim([0 1]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',16,'TickDir','out');
yyaxis right
h1 = boundedline(xplot,sim2RbinMu,sim2RbinSD/sqrt(sum(iLobV)),'cmap',[0 0 0],'alpha'); set(h1,'linewidth',1.5);
ylim([0 .2]);xlim([0 300]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',16,'TickDir','out');
title('Lobule simplex, AND logic');
axis square;

fSumXnor = figure;
subplot(1,2,1);
yyaxis left
h1 = boundedline(xplot,lobVvisXnorMu,lobVvisXnorSD/sqrt(sum(iLobV)),'cmap',[.1 .1 .8],'alpha'); set(h1,'linewidth',1.5);
h1 = boundedline(xplot,lobVmovXnorMu,lobVmovXnorSD/sqrt(sum(iLobV)),'cmap',[.1 .8 .1],'alpha'); set(h1,'linewidth',1.5);
h1 = boundedline(xplot,lobVrewXnorMu,lobVrewXnorSD/sqrt(sum(iLobV)),'cmap',[.1 .8 .8],'alpha'); set(h1,'linewidth',1.5);
h1 = boundedline(xplot,lobVlickXnorMu,lobVlickXnorSD/sqrt(sum(iLobV)),'cmap',[.5 .1 .8],'alpha'); set(h1,'linewidth',1.5);
ylim([0 1]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',16,'TickDir','out');
yyaxis right
h1 = boundedline(xplot,lobVRbinMu,lobVRbinSD/sqrt(sum(iLobV)),'cmap',[0 0 0],'alpha'); set(h1,'linewidth',1.5);
ylim([0 .2]);xlim([0 300]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',16,'TickDir','out');
title('Lobule V, XNOR logic');
axis square;

subplot(1,2,2);
yyaxis left
h1 = boundedline(xplot,sim2visXnorMu,sim2visXnorSD/sqrt(sum(iLobV)),'cmap',[.1 .1 .8],'alpha'); set(h1,'linewidth',1.5);
h1 = boundedline(xplot,sim2movXnorMu,sim2movXnorSD/sqrt(sum(iLobV)),'cmap',[.1 .8 .1],'alpha'); set(h1,'linewidth',1.5);
h1 = boundedline(xplot,sim2rewXnorMu,sim2rewXnorSD/sqrt(sum(iLobV)),'cmap',[.1 .8 .8],'alpha'); set(h1,'linewidth',1.5);
h1 = boundedline(xplot,sim2lickXnorMu,sim2lickXnorSD/sqrt(sum(iLobV)),'cmap',[.5 .1 .8],'alpha'); set(h1,'linewidth',1.5);
ylim([0 1]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',16,'TickDir','out');
yyaxis right
h1 = boundedline(xplot,sim2RbinMu,sim2RbinSD/sqrt(sum(iLobV)),'cmap',[0 0 0],'alpha'); set(h1,'linewidth',1.5);
ylim([0 .2]);xlim([0 300]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',16,'TickDir','out');
title('Lobule simplex, XNOR logic');
axis square;

if savebool
    fig1name = 'distDepEncodeAND';
    savefig(fSumAnd,fullfile(figFold,[fig1name,'.fig']));
    saveas(fSumAnd,fullfile(figFold,[fig1name,'.png']));
    print(fSumAnd,fullfile(figFold,[fig1name,'.eps']), '-depsc', '-painters');
    
    fig2name = 'distDepEncodeXNOR';
    savefig(fSumXnor,fullfile(figFold,[fig2name,'.fig']));
    saveas(fSumXnor,fullfile(figFold,[fig2name,'.png']));
    print(fSumXnor,fullfile(figFold,[fig2name,'.eps']), '-depsc', '-painters');

end

