figFold = '/Users/dimitar/Desktop/CF_learning_paper/paperFigs/Fig2/Sub1_GLM/DK105_lobv_v4_03-03';

plotStructs = load(fullfile(figFold,'plotStructs.mat'));

meanStruct = plotStructs.meanStruct;
nRois = size(meanStruct.Ytest,2);

f1 = figure;
subplot(1,3,1);
imagesc([mean(meanStruct.yOnPlot,3); mean(meanStruct.yOnHatPlot,3)]',[0 10]); colormap viridis;
hold on;
xline(3.5,'--w','LineWidth',2);
xline(33.5,'--w','LineWidth',2);
xline(30.5,'w','LineWidth',4);
xticks([3.5 18.5 33.5 48.5]);
xticklabels([0 0.5 0 0.5]);
yticks([1,40:40:nRois]);
xlabel('Time (s)');
ylabel('PC number');
title('Trial On: data (left), model (right)');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out'); box off;
subplot(1,3,2);
imagesc([mean(meanStruct.yMovPlot,3); mean(meanStruct.yMovHatPlot,3)]',[0 10]); colormap viridis;
hold on;
xline(15.5,'--w','LineWidth',2);
xline(45.5,'--w','LineWidth',2);
xline(30.5,'w','LineWidth',4);
xticks([0.5 15.5 30.5 45.5 60.5]);
xticklabels({'-0.5','0','(-)0.5','0','0.5'});
yticks([1,40:40:nRois]);
xlabel('Time (s)');
title('Movement: data (left), model (right)');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out'); box off;
subplot(1,3,3);
imagesc([mean(meanStruct.yREWcPlot,3); mean(meanStruct.yREWcHatPlot,3)]',[0 10]); colormap viridis;
hold on;
xline(15.5,'--w','LineWidth',2);
xline(45.5,'--w','LineWidth',2);
xline(30.5,'w','LineWidth',4);
xticks([0.5 15.5 30.5 45.5 60.5]);
xticklabels({'-0.5','0','(-)0.5','0','0.5'});
yticks([1,40:40:nRois]);
xlabel('Time (s)');
title('Trial reward: data (left), model (right)');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out'); box off;

figname = 'summaryHeatmap';
savefig(f1, fullfile(figFold,[figname,'.fig']));
saveas(f1,fullfile(figFold,[figname,'.png']));
print(f1,fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');

