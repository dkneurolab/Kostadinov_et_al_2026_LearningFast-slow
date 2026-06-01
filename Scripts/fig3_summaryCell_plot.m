function figsOut = fig3_summaryCell_plot(seshSummary,ilobv,isim2)

%% First figure: most movement stuff
figsOut(1) = figure;
% UCO
colMat = [158 31 99; 0 0 0; 194 117 20]/255;
dataMat = cat(1,seshSummary(ilobv).outcomeSpkCell);
subplot(2,5,1); pbaspect([1 1 1]); hold on;
plot([1,2,3],dataMat','Color',[0 0 0 .05],'LineWidth',.5);
e1 = errorbar(1:3,mean(dataMat),std(dataMat)/sqrt(size(dataMat,1)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 3.5 0 20]); box off;
xticks(1:3); xticklabels({'Under','Correct','Over'}); xtickangle(0)
yticks(0:2:20);
ylabel('Lobule V field rate (spk/s)');
xlabel('Outcome');

colMat = [158 31 99; 0 0 0; 194 117 20]/255;
dataMat = cat(1,seshSummary(isim2).outcomeSpkCell);
subplot(2,5,6); pbaspect([1 1 1]); hold on;
plot([1,2,3],dataMat','Color',[0 0 0 .05],'LineWidth',.5);
e1 = errorbar(1:3,mean(dataMat),std(dataMat)/sqrt(size(dataMat,1)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 3.5 0 20]); box off;
xticks(1:3); xticklabels({'Under','Correct','Over'}); xtickangle(0)
yticks(0:2:20);
ylabel('Simplex field rate (spk/s)');
xlabel('Outcome');

% xMax
colMat = repmat(linspace(.8,0,5),3,1)';
dataMat = cat(1,seshSummary(ilobv).xMaxSpkCell);
subplot(2,5,2); pbaspect([1 1 1]); hold on;
plot(1:5,dataMat','Color',[0 0 0 .05],'LineWidth',.5);
e1 = errorbar(1:5,mean(dataMat),std(dataMat)/sqrt(size(dataMat,1)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 5.5 0 20]); box off;
xticks(1:5); xticklabels({'0-20','20-40','40-60','60-80','80-100'}); xtickangle(30)
yticks(0:2:20);
% ylabel('Lobule V field rate (spk/s)');
xlabel('xMax percentile');

colMat = repmat(linspace(.8,0,5),3,1)';
dataMat = cat(1,seshSummary(isim2).xMaxSpkCell);
subplot(2,5,7); pbaspect([1 1 1]); hold on;
plot(1:5,dataMat','Color',[0 0 0 .05],'LineWidth',.5);
e1 = errorbar(1:5,mean(dataMat),std(dataMat)/sqrt(size(dataMat,1)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 5.5 0 20]); box off;
xticks(1:5); xticklabels({'0-20','20-40','40-60','60-80','80-100'}); xtickangle(30)
yticks(0:2:20);
% ylabel('Simplex field rate (spk/s)');
xlabel('xMax percentile');

% vMax
colMat = repmat(linspace(.8,0,5),3,1)';
dataMat = cat(1,seshSummary(ilobv).vMaxSpkCell);
subplot(2,5,3); pbaspect([1 1 1]); hold on;
plot(1:5,dataMat','Color',[0 0 0 .05],'LineWidth',.5);
e1 = errorbar(1:5,mean(dataMat),std(dataMat)/sqrt(size(dataMat,1)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 5.5 0 20]); box off;
xticks(1:5); xticklabels({'0-20','20-40','40-60','60-80','80-100'}); xtickangle(30)
yticks(0:2:20);
% ylabel('Lobule V field rate (spk/s)');
xlabel('vMax percentile');

colMat = repmat(linspace(.8,0,5),3,1)';
dataMat = cat(1,seshSummary(isim2).vMaxSpkCell);
subplot(2,5,8); pbaspect([1 1 1]); hold on;
plot(1:5,dataMat','Color',[0 0 0 .05],'LineWidth',.5);
e1 = errorbar(1:5,mean(dataMat),std(dataMat)/sqrt(size(dataMat,1)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 5.5 0 20]); box off;
xticks(1:5); xticklabels({'0-20','20-40','40-60','60-80','80-100'}); xtickangle(30)
yticks(0:2:20);
% ylabel('Simplex field rate (spk/s)');
xlabel('vMax percentile');

% aMax
colMat = repmat(linspace(.8,0,5),3,1)';
dataMat = cat(1,seshSummary(ilobv).aMaxSpkCell);
subplot(2,5,4); pbaspect([1 1 1]); hold on;
plot(1:5,dataMat','Color',[0 0 0 .05],'LineWidth',.5);
e1 = errorbar(1:5,mean(dataMat),std(dataMat)/sqrt(size(dataMat,1)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 5.5 0 20]); box off;
xticks(1:5); xticklabels({'0-20','20-40','40-60','60-80','80-100'}); xtickangle(30)
yticks(0:2:20);
% ylabel('Lobule V field rate (spk/s)');
xlabel('aMax percentile');

colMat = repmat(linspace(.8,0,5),3,1)';
dataMat = cat(1,seshSummary(isim2).aMaxSpkCell);
subplot(2,5,9); pbaspect([1 1 1]); hold on;
plot(1:5,dataMat','Color',[0 0 0 .05],'LineWidth',.5);
e1 = errorbar(1:5,mean(dataMat),std(dataMat)/sqrt(size(dataMat,1)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 5.5 0 20]); box off;
xticks(1:5); xticklabels({'0-20','20-40','40-60','60-80','80-100'}); xtickangle(30)
yticks(0:2:20);
% ylabel('Simplex field rate (spk/s)');
xlabel('aMax percentile');

% Rcorr
colMat = repmat(linspace(.8,0,5),3,1)';
dataMat = cat(1,seshSummary(ilobv).RcorrSpkCell);
subplot(2,5,5); pbaspect([1 1 1]); hold on;
plot(1:5,dataMat','Color',[0 0 0 .05],'LineWidth',.5);
e1 = errorbar(1:5,mean(dataMat),std(dataMat)/sqrt(size(dataMat,1)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 5.5 0 20]); box off;
xticks(1:5); xticklabels({'0-20','20-40','40-60','60-80','80-100'}); xtickangle(30)
yticks(0:2:20);
% ylabel('Lobule V field rate (spk/s)');
xlabel('Rcorr percentile');

colMat = repmat(linspace(.8,0,5),3,1)';
dataMat = cat(1,seshSummary(isim2).RcorrSpkCell);
subplot(2,5,10); pbaspect([1 1 1]); hold on;
plot(1:5,dataMat','Color',[0 0 0 .05],'LineWidth',.5);
e1 = errorbar(1:5,mean(dataMat),std(dataMat)/sqrt(size(dataMat,1)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 5.5 0 20]); box off;
xticks(1:5); xticklabels({'0-20','20-40','40-60','60-80','80-100'}); xtickangle(30)
yticks(0:2:20);
% ylabel('Simplex field rate (spk/s)');
xlabel('Rcorr percentile');

%% Figure 2: most reward stuff

figsOut(2) = figure;
% UCOLast
colMat = [158 31 99; 0 0 0; 194 117 20]/255;
dataMat = cat(1,seshSummary(ilobv).lastOutcomeSpkCell);
subplot(2,5,1); pbaspect([1 1 1]); hold on;
plot([1,2,3],dataMat','Color',[0 0 0 .05],'LineWidth',.5);
e1 = errorbar(1:3,mean(dataMat),std(dataMat)/sqrt(size(dataMat,1)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 3.5 0 12]); box off;
xticks(1:3); xticklabels({'Under','Correct','Over'}); xtickangle(0)
yticks(0:2:12);
ylabel('Lobule V reward response (spk/s)');
xlabel('Outcome of last trial');

colMat = [158 31 99; 0 0 0; 194 117 20]/255;
dataMat = cat(1,seshSummary(isim2).lastOutcomeSpkCell);
subplot(2,5,6); pbaspect([1 1 1]); hold on;
plot([1,2,3],dataMat','Color',[0 0 0 .05],'LineWidth',.5);
e1 = errorbar(1:3,mean(dataMat),std(dataMat)/sqrt(size(dataMat,1)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 3.5 0 12]); box off;
xticks(1:3); xticklabels({'Under','Correct','Over'}); xtickangle(0)
yticks(0:2:12);
ylabel('Simplex reward response (spk/s)');
xlabel('Outcome of last trial');

% UCOLast - binary
colMat = [158 31 99; 0 0 0; 194 117 20]/255;
dataMat = cat(1,seshSummary(ilobv).lastOutcomeBinSpkCell);
subplot(2,5,2); pbaspect([1 1 1]); hold on;
plot([1,2,3],dataMat','Color',[0 0 0 .05],'LineWidth',.5);
e1 = errorbar(1:3,mean(dataMat),std(dataMat)/sqrt(size(dataMat,1)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 3.5 0 20]); box off;
xticks(1:3); xticklabels({'Under','Correct','Over'}); xtickangle(0)
yticks(0:2:20);
% ylabel('Lobule V reward response (spk/s)');
xlabel('Outcome of last trial (binary)');

colMat = [158 31 99; 0 0 0; 194 117 20]/255;
dataMat = cat(1,seshSummary(isim2).lastOutcomeBinSpkCell);
subplot(2,5,7); pbaspect([1 1 1]); hold on;
plot([1,2,3],dataMat','Color',[0 0 0 .05],'LineWidth',.5);
e1 = errorbar(1:3,mean(dataMat),std(dataMat)/sqrt(size(dataMat,1)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 3.5 0 20]); box off;
xticks(1:3); xticklabels({'Under','Correct','Over'}); xtickangle(0)
yticks(0:2:20);
% ylabel('Simplex reward response (spk/s)');
xlabel('Outcome of last trial (binary)');

% spkMax - plot peak velocity
colMat = repmat(linspace(.8,0,5),3,1)';
dataMat = cat(2,seshSummary(ilobv).spkMaxvMax)';
subplot(2,5,4); pbaspect([1 1 1]); hold on;
plot(1:5,dataMat','Color',[0 0 0 .05],'LineWidth',.5);
e1 = errorbar(1:5,mean(dataMat),std(dataMat)/sqrt(size(dataMat,1)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 5.5 -2.5 15]); box off;
xticks(1:5); xticklabels({'0-20','20-40','40-60','60-80','80-100'}); xtickangle(30)
yticks(-2.5:2.5:15);
ylabel('Lobule V Peak velocity (cm/s)');
xlabel('Max activity percentile');

colMat = repmat(linspace(.8,0,5),3,1)';
dataMat = cat(2,seshSummary(isim2).spkMaxvMax)';
subplot(2,5,9); pbaspect([1 1 1]); hold on;
plot(1:5,dataMat','Color',[0 0 0 .05],'LineWidth',.5);
e1 = errorbar(1:5,mean(dataMat),std(dataMat)/sqrt(size(dataMat,1)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 5.5 -2.5 15]); box off;
xticks(1:5); xticklabels({'0-20','20-40','40-60','60-80','80-100'}); xtickangle(30)
yticks(-2.5:2.5:15);
ylabel('Simplex Peak velocity (cm/s)');
xlabel('Max activity percentile');

% syncMax - plot peak velocity
colMat = repmat(linspace(.8,0,5),3,1)';
dataMat = cat(2,seshSummary(ilobv).syncMaxvMax)';
subplot(2,5,5); pbaspect([1 1 1]); hold on;
plot(1:5,dataMat','Color',[0 0 0 .05],'LineWidth',.5);
e1 = errorbar(1:5,mean(dataMat),std(dataMat)/sqrt(size(dataMat,1)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 5.5 -2.5 15]); box off;
xticks(1:5); xticklabels({'0-20','20-40','40-60','60-80','80-100'}); xtickangle(30)
yticks(-2.5:2.5:15);
% ylabel('Lobule V Peak velocity (cm/s)');
xlabel('Max synchrony percentile');

colMat = repmat(linspace(.8,0,5),3,1)';
dataMat = cat(2,seshSummary(isim2).syncMaxvMax)';
subplot(2,5,10); pbaspect([1 1 1]); hold on;
plot(1:5,dataMat','Color',[0 0 0 .05],'LineWidth',.5);
e1 = errorbar(1:5,mean(dataMat),std(dataMat)/sqrt(size(dataMat,1)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 5.5 -2.5 15]); box off;
xticks(1:5); xticklabels({'0-20','20-40','40-60','60-80','80-100'}); xtickangle(30)
yticks(-2.5:2.5:15);
% ylabel('Simplex Peak velocity (cm/s)');
xlabel('Max synchrony percentile');

end

