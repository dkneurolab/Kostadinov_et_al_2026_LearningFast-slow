function figOut = fig3_matchFOVs_PlotFig(plot1,plot2,plot3,fov,plotTypes,seshTypes)
%% Make fig and populate one by one;
nRois = size(cat(2,plot1.spkMapSort),2);
t2p = ((1:120)/30)-1/60;
tvr = ((1:400)/100)-1/200;
figOut = figure;

%% Heatmap plots
% Plot 1
subplot(6,3,[1,4,7]);
spkMapSort = cat(2,plot1.spkMapSort);
coeff = pca(spkMapSort);
[~, iSort] = sort(coeff(:,1),'descend');
spkMap1 = cat(2,plot1.spkMap); spkMap1 = spkMap1(:,iSort);
imagesc(t2p,[1:nRois],spkMap1',[0 8]);
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks(0:.5:4); xticklabels({[],[],seshTypes{1},[],[],[],seshTypes{2},[],[],}); xtickangle(0);
yticks(0.5:50:nRois); yticklabels([1,50:50:nRois]);
hold on;
xline(1,'--w','LineWidth',1); xline(3,'--w','LineWidth',1); xline(2,'-w','LineWidth',2);
title(plotTypes(1));
ylabel([fov,' PC dendrite number']);

% Plot 2
subplot(6,3,[1,4,7]+1);
spkMapSort = cat(2,plot2.spkMapSort);
coeff = pca(spkMapSort);
[~, iSort] = sort(coeff(:,1),'descend');
spkMap2 = cat(2,plot2.spkMap); spkMap2 = spkMap2(:,iSort);
imagesc(t2p,[1:nRois],spkMap2',[0 8]);
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks(0:.5:4); xticklabels({[],[],seshTypes{1},[],[],[],seshTypes{2},[],[],}); xtickangle(0);
yticks(0.5:50:nRois); yticklabels([1,50:50:nRois]);
hold on;
xline(1,'--w','LineWidth',1); xline(3,'--w','LineWidth',1); xline(2,'-w','LineWidth',2);
title(plotTypes(2));
% ylabel([fov,' PC dendrite number']);

% Plot 3
subplot(6,3,[1,4,7]+2);
spkMapSort = cat(2,plot3.spkMapSort);
coeff = pca(spkMapSort);
[~, iSort] = sort(coeff(:,1),'descend');
spkMap3 = cat(2,plot3.spkMap); spkMap3 = spkMap3(:,iSort);
imagesc(t2p,[1:nRois],spkMap3',[0 8]);
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks(0:.5:4); xticklabels({[],[],seshTypes{1},[],[],[],seshTypes{2},[],[],}); xtickangle(0);
yticks(0.5:50:nRois); yticklabels([1,50:50:nRois]);
hold on;
xline(1,'--w','LineWidth',1); xline(3,'--w','LineWidth',1); xline(2,'-w','LineWidth',2);
title(plotTypes(3));
% ylabel([fov,' PC dendrite number']);

%% Average plots

subplot(6,3,10);
h = boundedline(t2p,mean(cat(2,plot1.spkMu),2), std(cat(2,plot1.spkMu),[],2)/sqrt(numel(plot1)),'-k','alpha');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(2,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 4 0 12]);
xticks([0:0.5:4]); xticklabels({-1,[],0,[],'(-)1',[],0,[],1,});
yticks([0:2:12]);
ylabel({'Pop. response','(spk/s)'});

subplot(6,3,11);
h = boundedline(t2p,mean(cat(2,plot2.spkMu),2), std(cat(2,plot2.spkMu),[],2)/sqrt(numel(plot2)),'-k','alpha');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(2,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 4 0 12]);
xticks([0:0.5:4]); xticklabels({-1,[],0,[],'(-)1',[],0,[],1,});
yticks([0:2:12]);

subplot(6,3,12);
h = boundedline(t2p,mean(cat(2,plot3.spkMu),2), std(cat(2,plot3.spkMu),[],2)/sqrt(numel(plot3)),'-k','alpha');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(2,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 4 0 12]);
xticks([0:0.5:4]); xticklabels({-1,[],0,[],'(-)1',[],0,[],1,});
yticks([0:2:12]);

%% Movement plots

subplot(6,3,13);
h = boundedline(tvr,mean(cat(2,plot1.movMu),2), std(cat(2,plot1.movMu),[],2)/sqrt(numel(plot1)),'-k','alpha');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(2,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 4 -2.5 10]);
xticks([0:0.5:4]); xticklabels({-1,[],0,[],'(-)1',[],0,[],1,});
yticks([-2.5:2.5:10]);
ylabel({'Wheel movement','(cm/s)'});

subplot(6,3,14);
h = boundedline(tvr,mean(cat(2,plot2.movMu),2), std(cat(2,plot2.movMu),[],2)/sqrt(numel(plot2)),'-k','alpha');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(2,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 4 -2.5 10]);
xticks([0:0.5:4]); xticklabels({-1,[],0,[],'(-)1',[],0,[],1,});
yticks([-2.5:2.5:10]);

subplot(6,3,15);
h = boundedline(tvr,mean(cat(2,plot3.movMu),2), std(cat(2,plot3.movMu),[],2)/sqrt(numel(plot3)),'-k','alpha');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(2,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 4 -2.5 10]);
xticks([0:0.5:4]); xticklabels({-1,[],0,[],'(-)1',[],0,[],1,});
yticks([-2.5:2.5:10]);

%% Licking plots

subplot(6,3,16);
h = boundedline(tvr,mean(cat(2,plot1.lickMu),2), std(cat(2,plot1.lickMu),[],2)/sqrt(numel(plot1)),'-k','alpha');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(2,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 4 0 1]);
xticks([0:0.5:4]); xticklabels({-1,[],0,[],'(-)1',[],0,[],1,});
yticks([0:.25:1]);
ylabel('Licking');
xlabel('Time (s)');

subplot(6,3,17);
h = boundedline(tvr,mean(cat(2,plot2.lickMu),2), std(cat(2,plot2.lickMu),[],2)/sqrt(numel(plot2)),'-k','alpha');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(2,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 4 0 1]);
xticks([0:0.5:4]); xticklabels({-1,[],0,[],'(-)1',[],0,[],1,});
yticks([0:.25:1]);
xlabel('Time (s)');

subplot(6,3,18);
h = boundedline(tvr,mean(cat(2,plot3.lickMu),2), std(cat(2,plot3.lickMu),[],2)/sqrt(numel(plot3)),'-k','alpha');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(2,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 4 0 1]);
xticks([0:0.5:4]); xticklabels({-1,[],0,[],'(-)1',[],0,[],1,});
yticks([0:.25:1]);

end