function fig2_GLM_plotTraces(sesh,GLMsig,eventsOut0,figFold0,savebool)

%% Smooth everybody!!
eventsOut = eventsOut0;
nSmooth = 3;
for iEv = 1:numel(eventsOut)
    eventsOut(iEv).yDataMu = smoothdata(eventsOut(iEv).yDataMu,1,'gaussian',nSmooth);
    eventsOut(iEv).yDataSD = smoothdata(eventsOut(iEv).yDataSD,1,'gaussian',nSmooth);
    eventsOut(iEv).yHatMu = smoothdata(eventsOut(iEv).yHatMu,1,'gaussian',nSmooth);
    eventsOut(iEv).yHatSD = smoothdata(eventsOut(iEv).yHatSD,1,'gaussian',nSmooth);
    eventsOut(iEv).yShufVisMu = smoothdata(eventsOut(iEv).yShufVisMu,1,'gaussian',nSmooth);
    eventsOut(iEv).yShufVisSD = smoothdata(eventsOut(iEv).yShufVisSD,1,'gaussian',nSmooth);
    eventsOut(iEv).yShufMovMu = smoothdata(eventsOut(iEv).yShufMovMu,1,'gaussian',nSmooth);
    eventsOut(iEv).yShufMovSD = smoothdata(eventsOut(iEv).yShufMovSD,1,'gaussian',nSmooth);
    eventsOut(iEv).yShufRewMu = smoothdata(eventsOut(iEv).yShufRewMu,1,'gaussian',nSmooth);
    eventsOut(iEv).yShufRewSD = smoothdata(eventsOut(iEv).yShufRewSD,1,'gaussian',nSmooth);
    eventsOut(iEv).yShufLickMu = smoothdata(eventsOut(iEv).yShufLickMu,1,'gaussian',nSmooth);
    eventsOut(iEv).yShufLickSD = smoothdata(eventsOut(iEv).yShufLickSD,1,'gaussian',nSmooth);
end

clear eventsOut0

%%
iCells = eventsOut(1).iCells;
GLMsig = GLMsig(iCells);

for iFig = 1:4:numel(iCells)
    iFigSub = iFig:iFig+3; iFigSub(iFigSub > numel(iCells)) = [];
    fPlot = figure;
    suptitle(sprintf('%s: %s - expert GLM fits',sesh.name,sesh.fov));
    for jFig = 1:numel(iFigSub)
        
        % Trial Onset aligned
        subplot(4,4,(jFig-1)*4+1)
        hold on;
        plot(eventsOut(1).tEvent,eventsOut(1).yHatMu(:,iFigSub(jFig)),'-','Color',[.5 .1 .1],'LineWidth',1);
        plot(eventsOut(1).tEvent,eventsOut(1).yShufVisMu(:,iFigSub(jFig)),'-','Color',[.1 .1 .8],'LineWidth',1);
        plot(eventsOut(1).tEvent,eventsOut(1).yShufMovMu(:,iFigSub(jFig)),'-','Color',[.1 .8 .1],'LineWidth',1);
        plot(eventsOut(1).tEvent,eventsOut(1).yShufRewMu(:,iFigSub(jFig)),'-','Color',[.1 .8 .8],'LineWidth',1);
        plot(eventsOut(1).tEvent,eventsOut(1).yShufLickMu(:,iFigSub(jFig)),'-','Color',[.5 .1 .8],'LineWidth',1);        
        h1 = boundedline(eventsOut(1).tEvent,eventsOut(1).yDataMu(:,iFigSub(jFig)),eventsOut(1).yDataSD(:,iFigSub(jFig))/sqrt(eventsOut(1).nEvents),'-k','alpha');
        set(h1,'linewidth',2);        
        xlim([-.2 .8]); ylim([0 15]);
        xticks([0:.2:.8]); yticks([0:3:15]);
        xline(0,'-','Color',[1 1 1]*.5);
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',16,'TickDir','out');
        title(sprintf('PC#%i,Trial ON',iCells(iFigSub(jFig))),'FontSize',11,'FontWeight','normal');
        t = text([.6 .6],[15 13],{sprintf('Full fit sig: %i',GLMsig(iFigSub(jFig)).sigBool);...
            sprintf('Vis fit sig: %i',GLMsig(iFigSub(jFig)).sigBool_shufVis_proj)},...
            'VerticalAlignment','top');
        t(2).Color = [.1 .1 .8];
        
        subplot(4,4,(jFig-1)*4+2)
        hold on;
        plot(eventsOut(2).tEvent,eventsOut(2).yHatMu(:,iFigSub(jFig)),'-','Color',[.5 .1 .1],'LineWidth',1);
        plot(eventsOut(2).tEvent,eventsOut(2).yShufVisMu(:,iFigSub(jFig)),'-','Color',[.1 .1 .8],'LineWidth',1);
        plot(eventsOut(2).tEvent,eventsOut(2).yShufMovMu(:,iFigSub(jFig)),'-','Color',[.1 .8 .1],'LineWidth',1);
        plot(eventsOut(2).tEvent,eventsOut(2).yShufRewMu(:,iFigSub(jFig)),'-','Color',[.1 .8 .8],'LineWidth',1);
        plot(eventsOut(2).tEvent,eventsOut(2).yShufLickMu(:,iFigSub(jFig)),'-','Color',[.5 .1 .8],'LineWidth',1);        
        h1 = boundedline(eventsOut(2).tEvent,eventsOut(2).yDataMu(:,iFigSub(jFig)),eventsOut(2).yDataSD(:,iFigSub(jFig))/sqrt(eventsOut(2).nEvents),'-k','alpha');
        set(h1,'linewidth',2);        
        xlim([-.5 .5]); ylim([0 15]);
        xticks([-.4:.2:.4]); yticks([0:3:15]);
        xline(0,'-','Color',[1 1 1]*.5);
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',16,'TickDir','out');
        title('Movement','FontSize',11,'FontWeight','normal');
        text(.3,15,{sprintf('Mov fit sig: %i',GLMsig(iFigSub(jFig)).sigBool_shufMov_proj)},...
            'VerticalAlignment','top','Color',[.1 .8 .1]);
        
        subplot(4,4,(jFig-1)*4+3)
        hold on;
        plot(eventsOut(3).tEvent,eventsOut(3).yHatMu(:,iFigSub(jFig)),'-','Color',[.5 .1 .1],'LineWidth',1);
        plot(eventsOut(3).tEvent,eventsOut(3).yShufVisMu(:,iFigSub(jFig)),'-','Color',[.1 .1 .8],'LineWidth',1);
        plot(eventsOut(3).tEvent,eventsOut(3).yShufMovMu(:,iFigSub(jFig)),'-','Color',[.1 .8 .1],'LineWidth',1);
        plot(eventsOut(3).tEvent,eventsOut(3).yShufRewMu(:,iFigSub(jFig)),'-','Color',[.1 .8 .8],'LineWidth',1);
        plot(eventsOut(3).tEvent,eventsOut(3).yShufLickMu(:,iFigSub(jFig)),'-','Color',[.5 .1 .8],'LineWidth',1);        
        h1 = boundedline(eventsOut(3).tEvent,eventsOut(3).yDataMu(:,iFigSub(jFig)),eventsOut(3).yDataSD(:,iFigSub(jFig))/sqrt(eventsOut(3).nEvents),'-k','alpha');
        set(h1,'linewidth',2);        
        xlim([-.5 .5]); ylim([0 15]);
        xticks([-.4:.2:.4]); yticks([0:3:15]);
        xline(0,'-','Color',[1 1 1]*.5);
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',16,'TickDir','out');
        title('Trial Reward','FontSize',11,'FontWeight','normal');
        text(.3,15,{sprintf('Rew fit sig: %i',GLMsig(iFigSub(jFig)).sigBool_shufRew_proj)},...
            'VerticalAlignment','top','Color',[.1 .8 .8]);
        
        subplot(4,4,(jFig-1)*4+4)
        hold on;
        plot(eventsOut(4).tEvent,eventsOut(4).yHatMu(:,iFigSub(jFig)),'-','Color',[.5 .1 .1],'LineWidth',1);
        plot(eventsOut(4).tEvent,eventsOut(4).yShufVisMu(:,iFigSub(jFig)),'-','Color',[.1 .1 .8],'LineWidth',1);
        plot(eventsOut(4).tEvent,eventsOut(4).yShufMovMu(:,iFigSub(jFig)),'-','Color',[.1 .8 .1],'LineWidth',1);
        plot(eventsOut(4).tEvent,eventsOut(4).yShufRewMu(:,iFigSub(jFig)),'-','Color',[.1 .8 .8],'LineWidth',1);
        plot(eventsOut(4).tEvent,eventsOut(4).yShufLickMu(:,iFigSub(jFig)),'-','Color',[.5 .1 .8],'LineWidth',1);        
        h1 = boundedline(eventsOut(4).tEvent,eventsOut(4).yDataMu(:,iFigSub(jFig)),eventsOut(4).yDataSD(:,iFigSub(jFig))/sqrt(eventsOut(4).nEvents),'-k','alpha');
        set(h1,'linewidth',2);        
        xlim([-.5 .5]); ylim([0 15]);
        xticks([-.4:.2:.4]); yticks([0:3:15]);
        xline(0,'-','Color',[1 1 1]*.5);
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',16,'TickDir','out');
        title('Licking','FontSize',11,'FontWeight','normal');
        text(.3,15,{sprintf('Lick fit sig: %i',GLMsig(iFigSub(jFig)).sigBool_shufLick_proj)},...
            'VerticalAlignment','top','Color',[.5 .1 .8]);
        
        
    end
    if savebool
        figName = sprintf('%s_%s_PCs_%03.f_to_%03.f',sesh.name,sesh.fov,iCells(iFigSub(1)),iCells(iFigSub(end)));
        figFold = fullfile(figFold0,sprintf('%s_%s',sesh.name,sesh.fov));
        if ~exist(figFold,'dir'); mkdir(figFold); end        
        savefig(fPlot,fullfile(figFold,[figName,'_Traces.fig']));
        saveas(fPlot,fullfile(figFold,[figName,'_Traces.png']));
        print(fPlot,fullfile(figFold,[figName,'_Traces.eps']), '-depsc', '-painters');
    end
end
