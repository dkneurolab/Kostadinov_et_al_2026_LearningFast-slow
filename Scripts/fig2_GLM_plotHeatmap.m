function figOut = fig2_GLM_plotHeatmap(session,eventsOut,nCells,t0,figName)


figOut = figure;
subplot(1,6,1);
imagesc(eventsOut(1).yDataMu',[0 8]); box off
yticks(.5:50:nCells); xticks(0.5:3:30.5); xticklabels({}); yticklabels({});
xline(t0,'--w');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
title('Data');
subplot(1,6,2);
imagesc(eventsOut(1).yHatMu',[0 8]); box off
yticks(.5:50:nCells); xticks(0.5:3:30.5); xticklabels({}); yticklabels({});
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
title('Full model');
xline(t0,'--w');
subplot(1,6,3);
imagesc(eventsOut(1).yShufVisMu',[0 8]); box off
yticks(.5:50:nCells); xticks(0.5:3:30.5); xticklabels({}); yticklabels({});
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
xline(t0,'--w');
title('Vis shuf');
subplot(1,6,4);
imagesc(eventsOut(1).yShufMovMu',[0 8]); box off
yticks(.5:50:nCells); xticks(0.5:3:30.5); xticklabels({}); yticklabels({});
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
xline(t0,'--w');
title('Mov shuf');
subplot(1,6,5);
imagesc(eventsOut(1).yShufRewMu',[0 8]); box off
yticks(.5:50:nCells); xticks(0.5:3:30.5); xticklabels({}); yticklabels({});
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
xline(t0,'--w');
title('Rew shuf');
subplot(1,6,6);
imagesc(eventsOut(1).yShufLickMu',[0 8]); box off
yticks(.5:50:nCells); xticks(0.5:3:30.5); xticklabels({}); yticklabels({});
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
xline(t0,'--w');
title('Lick shuf');
colormap viridis
suptitle(sprintf('%s: %s, %s',session.name,session.fov,figName));