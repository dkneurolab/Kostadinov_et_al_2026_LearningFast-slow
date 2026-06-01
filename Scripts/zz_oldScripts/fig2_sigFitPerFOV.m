function fig2_sigFitPerFOV(glmSessions,dataFold_sub,paperFold_sub,savebool)


iLobv = false(numel(glmSessions),1);
iSim2 = false(numel(glmSessions),1);
percEV = zeros(numel(glmSessions),1);
percRsq = zeros(numel(glmSessions),1);
for i = 1:numel(glmSessions)
    if strcmpi(glmSessions(i).fov, 'lobv')
        iLobv(i) = true;
    else
        iSim2(i) = true;
    end
    
    percEV(i) = sum(glmSessions(i).statStruct.EV_sig)/numel(glmSessions(i).statStruct.EV_sig);
    percRsq(i) = sum(glmSessions(i).statStruct.Rsq_sig)/numel(glmSessions(i).statStruct.Rsq_sig);
    
    
    
    
end

f1 = figure; 
subplot(2,2,1); hold on;
bar(1,mean(percEV(iLobv)),0.5,'EdgeColor','none','FaceColor',[.6 .0 .6],'FaceAlpha',0.5);
errorbar(1, mean(percEV(iLobv)),std(percEV(iLobv)),'k','CapSize',0,'LineWidth',2);
bar(3,mean(percEV(iSim2)),0.5,'EdgeColor','none','FaceColor',[.0 .6 .6],'FaceAlpha',0.5);
errorbar(3, mean(percEV(iSim2)),std(percEV(iSim2)),'k','CapSize',0,'LineWidth',2);
ylim([0 1.1]); yticks([0:0.2:1]);
xlim([0 4]); xticks([1,3]); xticklabels({'Lobule V', 'Simplex'})
ylabel({'Proportion of PC dendrites with'; 'significant EV fit'});
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
axis square;

subplot(2,2,2); hold on;
bar(1,mean(percRsq(iLobv)),0.5,'EdgeColor','none','FaceColor',[.6 .0 .6],'FaceAlpha',0.5);
errorbar(1, mean(percRsq(iLobv)),std(percRsq(iLobv)),'k','CapSize',0,'LineWidth',2);
bar(3,mean(percRsq(iSim2)),0.5,'EdgeColor','none','FaceColor',[.0 .6 .6],'FaceAlpha',0.5);
errorbar(3, mean(percRsq(iSim2)),std(percRsq(iSim2)),'k','CapSize',0,'LineWidth',2);
ylim([0 1.1]); yticks([0:0.2:1]);
xlim([0 4]); xticks([1,3]); xticklabels({'Lobule V', 'Simplex'})
ylabel({'Proportion of PC dendrites with'; 'significant Rsq fit'});
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
axis square;

if savebool
    savefig(f1,fullfile(dataFold_sub,'sigFitPerFOV.fig'));
    saveas(f1,fullfile(dataFold_sub,'sigFitPerFOV.png'));
    print(f1,fullfile(dataFold_sub,'sigFitPerFOV.eps'),'-depsc','-painters')
    
    savefig(f1,fullfile(paperFold_sub,'sigFitPerFOV.fig'));
    saveas(f1,fullfile(paperFold_sub,'sigFitPerFOV.png'));
    print(f1,fullfile(paperFold_sub,'sigFitPerFOV.eps'),'-depsc','-painters')

end
    
