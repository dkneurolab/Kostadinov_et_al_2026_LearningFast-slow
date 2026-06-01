function [figsOut, statsOut] = fig3_summary_plot(seshSummary,ilobv,isim2)

lobvMovMax = max(cat(2,seshSummary(ilobv).spkMovTrace));
sim2MovMax = max(cat(2,seshSummary(isim2).spkMovTrace));
lobvRewMax = max(cat(2,seshSummary(ilobv).spkRewTrace));
sim2RewMax = max(cat(2,seshSummary(isim2).spkRewTrace));

statsOut = struct;

%% First figure: most movement stuff
figsOut(1) = figure;
% UCO
colMat = [158 31 99; 0 0 0; 194 117 20]/255;
dataMat = cat(2,seshSummary(ilobv).outcomeSpk)';
[p,tbl,stats] = friedman(dataMat);
c = multcompare(stats,'CType', 'bonferroni');
statsOut.lobv_outcomeSpk = v2struct(p,tbl,stats,c);
subplot(2,5,1); pbaspect([1 1 1]); hold on;
plot([1,2,3],dataMat','Color',[0 0 0 .5],'LineWidth',.5);
e1 = errorbar(1:3,mean(dataMat),std(dataMat)/sqrt(numel(ilobv)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 3.5 0 8]); box off;
xticks(1:3); xticklabels({'Under','Correct','Over'}); xtickangle(0)
yticks(0:2:8);
ylabel('Lobule V field rate (spk/s)');
xlabel('Outcome');

colMat = [158 31 99; 0 0 0; 194 117 20]/255;
dataMat = cat(2,seshSummary(isim2).outcomeSpk)';
[p,tbl,stats] = friedman(dataMat);
c = multcompare(stats,'CType', 'bonferroni');
statsOut.sim2_outcomeSpk = v2struct(p,tbl,stats,c);
subplot(2,5,6); pbaspect([1 1 1]); hold on;
plot([1,2,3],dataMat','Color',[0 0 0 .5],'LineWidth',.5);
e1 = errorbar(1:3,mean(dataMat),std(dataMat)/sqrt(numel(isim2)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 3.5 0 8]); box off;
xticks(1:3); xticklabels({'Under','Correct','Over'}); xtickangle(0)
yticks(0:2:8);
ylabel('Simplex field rate (spk/s)');
xlabel('Outcome');

% xMax
colMat = repmat(linspace(.8,0,5),3,1)';
dataMat = cat(2,seshSummary(ilobv).xMaxSpk)';
[p,tbl,stats] = friedman(dataMat);
c = multcompare(stats,'CType', 'bonferroni');
statsOut.lobv_xMaxSpk = v2struct(p,tbl,stats,c);
subplot(2,5,2); pbaspect([1 1 1]); hold on;
plot(1:5,dataMat','Color',[0 0 0 .5],'LineWidth',.5);
e1 = errorbar(1:5,mean(dataMat),std(dataMat)/sqrt(numel(ilobv)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 5.5 0 8]); box off;
xticks(1:5); xticklabels({'0-20','20-40','40-60','60-80','80-100'}); xtickangle(30)
yticks(0:2:8);
% ylabel('Lobule V field rate (spk/s)');
xlabel('xMax percentile');

colMat = repmat(linspace(.8,0,5),3,1)';
dataMat = cat(2,seshSummary(isim2).xMaxSpk)';
[p,tbl,stats] = friedman(dataMat);
c = multcompare(stats,'CType', 'bonferroni');
statsOut.sim2_xMaxSpk = v2struct(p,tbl,stats,c);
subplot(2,5,7); pbaspect([1 1 1]); hold on;
plot(1:5,dataMat','Color',[0 0 0 .5],'LineWidth',.5);
e1 = errorbar(1:5,mean(dataMat),std(dataMat)/sqrt(numel(isim2)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 5.5 0 8]); box off;
xticks(1:5); xticklabels({'0-20','20-40','40-60','60-80','80-100'}); xtickangle(30)
yticks(0:2:8);
% ylabel('Simplex field rate (spk/s)');
xlabel('xMax percentile');

% vMax
colMat = repmat(linspace(.8,0,5),3,1)';
dataMat = cat(2,seshSummary(ilobv).vMaxSpk)';
[p,tbl,stats] = friedman(dataMat);
c = multcompare(stats,'CType', 'bonferroni');
statsOut.lobv_vMaxSpk = v2struct(p,tbl,stats,c);
subplot(2,5,3); pbaspect([1 1 1]); hold on;
plot(1:5,dataMat','Color',[0 0 0 .5],'LineWidth',.5);
e1 = errorbar(1:5,mean(dataMat),std(dataMat)/sqrt(numel(ilobv)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 5.5 0 8]); box off;
xticks(1:5); xticklabels({'0-20','20-40','40-60','60-80','80-100'}); xtickangle(30)
yticks(0:2:8);
% ylabel('Lobule V field rate (spk/s)');
xlabel('vMax percentile');

colMat = repmat(linspace(.8,0,5),3,1)';
dataMat = cat(2,seshSummary(isim2).vMaxSpk)';
[p,tbl,stats] = friedman(dataMat);
c = multcompare(stats,'CType', 'bonferroni');
statsOut.sim2_vMaxSpk = v2struct(p,tbl,stats,c);
subplot(2,5,8); pbaspect([1 1 1]); hold on;
plot(1:5,dataMat','Color',[0 0 0 .5],'LineWidth',.5);
e1 = errorbar(1:5,mean(dataMat),std(dataMat)/sqrt(numel(isim2)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 5.5 0 8]); box off;
xticks(1:5); xticklabels({'0-20','20-40','40-60','60-80','80-100'}); xtickangle(30)
yticks(0:2:8);
% ylabel('Simplex field rate (spk/s)');
xlabel('vMax percentile');

% aMax
colMat = repmat(linspace(.8,0,5),3,1)';
dataMat = cat(2,seshSummary(ilobv).aMaxSpk)';
[p,tbl,stats] = friedman(dataMat);
c = multcompare(stats,'CType', 'bonferroni');
statsOut.lobv_aMaxSpk = v2struct(p,tbl,stats,c);
subplot(2,5,4); pbaspect([1 1 1]); hold on;
plot(1:5,dataMat','Color',[0 0 0 .5],'LineWidth',.5);
e1 = errorbar(1:5,mean(dataMat),std(dataMat)/sqrt(numel(ilobv)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 5.5 0 8]); box off;
xticks(1:5); xticklabels({'0-20','20-40','40-60','60-80','80-100'}); xtickangle(30)
yticks(0:2:8);
% ylabel('Lobule V field rate (spk/s)');
xlabel('aMax percentile');

colMat = repmat(linspace(.8,0,5),3,1)';
dataMat = cat(2,seshSummary(isim2).aMaxSpk)';
[p,tbl,stats] = friedman(dataMat);
c = multcompare(stats,'CType', 'bonferroni');
statsOut.sim2_aMaxSpk = v2struct(p,tbl,stats,c);
subplot(2,5,9); pbaspect([1 1 1]); hold on;
plot(1:5,dataMat','Color',[0 0 0 .5],'LineWidth',.5);
e1 = errorbar(1:5,mean(dataMat),std(dataMat)/sqrt(numel(isim2)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 5.5 0 8]); box off;
xticks(1:5); xticklabels({'0-20','20-40','40-60','60-80','80-100'}); xtickangle(30)
yticks(0:2:8);
% ylabel('Simplex field rate (spk/s)');
xlabel('aMax percentile');

% Rcorr
colMat = repmat(linspace(.8,0,5),3,1)';
dataMat = cat(2,seshSummary(ilobv).RcorrSpkB)';
[p,tbl,stats] = friedman(dataMat(:,1:3));
c = multcompare(stats,'CType', 'bonferroni');
statsOut.lobv_RcorrSpkB = v2struct(p,tbl,stats,c);
subplot(2,5,5); pbaspect([1 1 1]); hold on;
plot(1:5,dataMat','Color',[0 0 0 .5],'LineWidth',.5);
e1 = errorbar(1:5,mean(dataMat,'omitnan'),std(dataMat,'omitnan')/sqrt(numel(ilobv)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData),'omitnan'),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 5.5 0 8]); box off;
xticks(1:5); xticklabels({'0-20','20-40','40-60','60-80','80-100'}); xtickangle(30)
yticks(0:2:8);
% ylabel('Lobule V field rate (spk/s)');
xlabel('Rcorr percentile');

colMat = repmat(linspace(.8,0,5),3,1)';
dataMat = cat(2,seshSummary(isim2).RcorrSpkB)';
[p,tbl,stats] = friedman(dataMat(:,1:3));
c = multcompare(stats,'CType', 'bonferroni');
statsOut.sim2_RcorrSpkB = v2struct(p,tbl,stats,c);
subplot(2,5,10); pbaspect([1 1 1]); hold on;
plot(1:5,dataMat','Color',[0 0 0 .5],'LineWidth',.5);
e1 = errorbar(1:5,mean(dataMat),std(dataMat)/sqrt(numel(isim2)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 5.5 0 8]); box off;
xticks(1:5); xticklabels({'0-20','20-40','40-60','60-80','80-100'}); xtickangle(30)
yticks(0:2:8);
% ylabel('Simplex field rate (spk/s)');
xlabel('Rcorr percentile');

%% Figure 2: most reward stuff

figsOut(2) = figure;
% UCOLast
colMat = [158 31 99; 0 0 0; 194 117 20]/255;
dataMat = cat(2,seshSummary(ilobv).lastOutcomeSpk)'; dataMat = dataMat(:,1:3);
[p,tbl,stats] = friedman(dataMat);
c = multcompare(stats,'CType', 'bonferroni');
statsOut.lobv_lastOutcomeSpk = v2struct(p,tbl,stats,c);
subplot(2,5,1); pbaspect([1 1 1]); hold on;
plot([1,2,3],dataMat','Color',[0 0 0 .5],'LineWidth',.5);
e1 = errorbar(1:3,mean(dataMat),std(dataMat)/sqrt(numel(ilobv)));
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
dataMat = cat(2,seshSummary(isim2).lastOutcomeSpk)';  dataMat = dataMat(:,1:3);
[p,tbl,stats] = friedman(dataMat);
c = multcompare(stats,'CType', 'bonferroni');
statsOut.sim2_lastOutcomeSpk = v2struct(p,tbl,stats,c);
subplot(2,5,6); pbaspect([1 1 1]); hold on;
plot([1,2,3],dataMat','Color',[0 0 0 .5],'LineWidth',.5);
e1 = errorbar(1:3,mean(dataMat),std(dataMat)/sqrt(numel(isim2)));
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
dataMat = cat(2,seshSummary(ilobv).lastOutcomeBinSpk)'; dataMat = dataMat(:,1:3);
[p,tbl,stats] = friedman(dataMat);
c = multcompare(stats,'CType', 'bonferroni');
statsOut.lobv_lastOutcomeBinSpk = v2struct(p,tbl,stats,c);
subplot(2,5,2); pbaspect([1 1 1]); hold on;
plot([1,2,3],dataMat','Color',[0 0 0 .5],'LineWidth',.5);
e1 = errorbar(1:3,mean(dataMat),std(dataMat)/sqrt(numel(ilobv)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 3.5 0 12]); box off;
xticks(1:3); xticklabels({'Under','Correct','Over'}); xtickangle(0)
yticks(0:2:12);
% ylabel('Lobule V reward response (spk/s)');
xlabel('Outcome of last trial (binary)');

colMat = [158 31 99; 0 0 0; 194 117 20]/255;
dataMat = cat(2,seshSummary(isim2).lastOutcomeBinSpk)';  dataMat = dataMat(:,1:3);
[p,tbl,stats] = friedman(dataMat);
c = multcompare(stats,'CType', 'bonferroni');
statsOut.sim2_lastOutcomeBinSpk = v2struct(p,tbl,stats,c);
subplot(2,5,7); pbaspect([1 1 1]); hold on;
plot([1,2,3],dataMat','Color',[0 0 0 .5],'LineWidth',.5);
e1 = errorbar(1:3,mean(dataMat),std(dataMat)/sqrt(numel(isim2)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 3.5 0 12]); box off;
xticks(1:3); xticklabels({'Under','Correct','Over'}); xtickangle(0)
yticks(0:2:12);
% ylabel('Simplex reward response (spk/s)');
xlabel('Outcome of last trial (binary)');

% spkMax - plot peak velocity
colMat = repmat(linspace(.8,0,5),3,1)';
dataMat = cat(2,seshSummary(ilobv).spkMaxvMax)';
[p,tbl,stats] = friedman(dataMat);
c = multcompare(stats,'CType', 'bonferroni');
statsOut.lobv_spkMaxvMax = v2struct(p,tbl,stats,c);
subplot(2,5,4); pbaspect([1 1 1]); hold on;
plot(1:5,dataMat','Color',[0 0 0 .5],'LineWidth',.5);
e1 = errorbar(1:5,mean(dataMat),std(dataMat)/sqrt(numel(ilobv)));
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
[p,tbl,stats] = friedman(dataMat);
c = multcompare(stats,'CType', 'bonferroni');
statsOut.sim2_spkMaxvMax = v2struct(p,tbl,stats,c);
subplot(2,5,9); pbaspect([1 1 1]); hold on;
plot(1:5,dataMat','Color',[0 0 0 .5],'LineWidth',.5);
e1 = errorbar(1:5,mean(dataMat),std(dataMat)/sqrt(numel(ilobv)));
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
[p,tbl,stats] = friedman(dataMat);
c = multcompare(stats,'CType', 'bonferroni');
statsOut.lobv_syncMaxvMax = v2struct(p,tbl,stats,c);
subplot(2,5,5); pbaspect([1 1 1]); hold on;
plot(1:5,dataMat','Color',[0 0 0 .5],'LineWidth',.5);
e1 = errorbar(1:5,mean(dataMat),std(dataMat)/sqrt(numel(ilobv)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 5.5 -2.5 15]); box off;
xticks(1:5); xticklabels({'0-20','20-40','40-60','60-80','80-100'}); xtickangle(30)
yticks(-2.5:2.5:15);
ylabel('Lobule V Peak velocity (cm/s)');
xlabel('Max synchrony percentile');

colMat = repmat(linspace(.8,0,5),3,1)';
dataMat = cat(2,seshSummary(isim2).syncMaxvMax)';
[p,tbl,stats] = friedman(dataMat);
c = multcompare(stats,'CType', 'bonferroni');
statsOut.sim2_syncMaxvMax = v2struct(p,tbl,stats,c);
subplot(2,5,10); pbaspect([1 1 1]); hold on;
plot(1:5,dataMat','Color',[0 0 0 .5],'LineWidth',.5);
e1 = errorbar(1:5,mean(dataMat),std(dataMat)/sqrt(numel(ilobv)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 5.5 -2.5 15]); box off;
xticks(1:5); xticklabels({'0-20','20-40','40-60','60-80','80-100'}); xtickangle(30)
yticks(-2.5:2.5:15);
ylabel('Simplex Peak velocity (cm/s)');
xlabel('Max synchrony percentile');

%% Fig3 - normalized Rcorr and reward history

figsOut(3) = figure;
% Rcorr
colMat = repmat(linspace(.8,0,5),3,1)';
dataMat = (cat(2,seshSummary(ilobv).RcorrSpkB)')./lobvMovMax';
[p,tbl,stats] = friedman(dataMat(:,1:3));
c = multcompare(stats,'CType', 'bonferroni');
statsOut.lobv_RcorrSpkB_norm = v2struct(p,tbl,stats,c);
subplot(2,2,1); pbaspect([.7 1 1]); hold on;
plot(1:5,dataMat','Color',[0 0 0 .5],'LineWidth',.5);
e1 = errorbar(1:5,mean(dataMat,'omitnan'),std(dataMat,'omitnan')/sqrt(numel(ilobv)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData),'omitnan'),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 3.5 .5 1.5]); box off;
xticks(1:3); xticklabels({'0-0.5','0.5-0.9','0.9-1'}); xtickangle(30)
yticks(.5:.25:1.5);
ylabel('Lobule V field rate (norm)');
% xlabel('Rcorr grouping');

colMat = repmat(linspace(.8,0,5),3,1)';
dataMat = (cat(2,seshSummary(isim2).RcorrSpkB)')./sim2MovMax';
[p,tbl,stats] = friedman(dataMat(:,1:3));
c = multcompare(stats,'CType', 'bonferroni');
statsOut.sim2_RcorrSpkB_norm = v2struct(p,tbl,stats,c);
subplot(2,2,3); pbaspect([.7 1 1]); hold on;
plot(1:5,dataMat','Color',[0 0 0 .5],'LineWidth',.5);
e1 = errorbar(1:5,mean(dataMat,'omitnan'),std(dataMat,'omitnan')/sqrt(numel(isim2)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData),'omitnan'),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 3.5 .5 1.5]); box off;
xticks(1:3); xticklabels({'0-0.5','0.5-0.9','0.9-1'}); xtickangle(30)
yticks(.5:.25:1.5);
ylabel('Simplex field rate (norm)');
xlabel('Rcorr grouping');

% UCOLast
colMat = [158 31 99; 0 0 0; 194 117 20]/255;
dataMat = cat(2,seshSummary(ilobv).lastOutcomeSpk)'; dataMat = (dataMat(:,1:3)./lobvRewMax');
[p,tbl,stats] = friedman(dataMat(:,1:3));
c = multcompare(stats,'CType', 'bonferroni');
statsOut.lobv_lastOutcomeSpk_norm = v2struct(p,tbl,stats,c);
subplot(2,2,2); pbaspect([.7 1 1]); hold on;
plot([1,2,3],dataMat','Color',[0 0 0 .5],'LineWidth',.5);
e1 = errorbar(1:3,mean(dataMat),std(dataMat)/sqrt(numel(ilobv)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 3.5 .5 1.5]); box off;
xticks(1:3); xticklabels({'Under','Correct','Over'}); xtickangle(0)
yticks(.5:.25:1.5);
ylabel('Lobule V reward response (spk/s)');
% xlabel('Outcome of last trial');

colMat = [158 31 99; 0 0 0; 194 117 20]/255;
dataMat = cat(2,seshSummary(isim2).lastOutcomeSpk)'; dataMat = (dataMat(:,1:3)./sim2RewMax');
[p,tbl,stats] = friedman(dataMat(:,1:3));
c = multcompare(stats,'CType', 'bonferroni');
statsOut.sim2_lastOutcomeSpk_norm = v2struct(p,tbl,stats,c);
subplot(2,2,4); pbaspect([.7 1 1]); hold on;
plot([1,2,3],dataMat','Color',[0 0 0 .5],'LineWidth',.5);
e1 = errorbar(1:3,mean(dataMat),std(dataMat)/sqrt(numel(isim2)));
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
for iData = 1:size(dataMat,2)
    scatter(iData,mean(dataMat(:,iData)),60,'MarkerEdgeColor','none','MarkerFaceColor',colMat(iData,:));
end
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([.5 3.5 .5 1.5]); box off;
xticks(1:3); xticklabels({'Under','Correct','Over'}); xtickangle(0)
yticks(.5:.25:1.5);
ylabel('Simplex reward response (spk/s)');
xlabel('Outcome of last trial');

end

