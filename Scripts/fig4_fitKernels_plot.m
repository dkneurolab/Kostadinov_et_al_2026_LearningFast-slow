function [figOut,statsOut] = fig4_fitKernels_plot(GLMdata,plotType)
%% Make data
% v1
kernelMatv1 = cat(2,GLMdata.(sprintf('%sKernel_v1',plotType)));
iGood = ~isnan(kernelMatv1(1,:));
kernelMatv1 = kernelMatv1(:,iGood);

% v41
kernelMatv41 = cat(2,GLMdata.(sprintf('%sKernel_v41',plotType)));
kernelMatv41 = kernelMatv41(:,iGood);

% v4n
kernelMatv4n = cat(2,GLMdata.(sprintf('%sKernel_v4n',plotType)));
kernelMatv4n = kernelMatv4n(:,iGood);

% Sort based on v4n
pcaCoeff = pca(kernelMatv4n);
[~, iSort] = sort(pcaCoeff(:,1),'descend');
kernelMatv1 = kernelMatv1(:,iSort);
kernelMatv41 = kernelMatv41(:,iSort);
kernelMatv4n = kernelMatv4n(:,iSort);

xmap = [1:size(kernelMatv1,1)]';

% For naming later
fov = inputname(1);
fov = fov(end-3:end);
nPCs = size(kernelMatv1,2);

% Calculate correlations
Rmat = corrcoef([kernelMatv1 kernelMatv41 kernelMatv4n]);
Rv1v41All = Rmat(1:nPCs,nPCs+1:2*nPCs); % all v1-v41
Rv1v41Matched = Rv1v41All(1:size(Rv1v41All,1)+1:end); % matched cells v1-v41
Rv1v4nAll = Rmat(1:nPCs,2*nPCs+1:3*nPCs); % all v1-v4n
Rv1v4nMatched = Rv1v4nAll(1:size(Rv1v4nAll,1)+1:end); % matched cells v1-v4n
Rv41v4nAll = Rmat(nPCs+1:2*nPCs,2*nPCs+1:3*nPCs); % all v41-v4n
Rv41v4nMatched = Rv41v4nAll(1:size(Rv41v4nAll,1)+1:end); % matched cells v41-v4n


% Rmat = pdist2([kernelMatv1 kernelMatv41 kernelMatv4n],[kernelMatv1 kernelMatv41 kernelMatv4n],'cosine')

%% Plot stuff
figOut = figure;
% Heatmaps
subplot(4,6,[1:6:13]);
imagesc(kernelMatv1',[0 .5]); colormap viridis;
axis image; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
yticks([.5:20:nPCs]);
xticks([.5 15.5 30.5]); xtickangle(0);
title(sprintf('%s-%s-v1',fov,plotType));

subplot(4,6,[1:6:13]+1);
imagesc(kernelMatv41',[0 .5]); colormap viridis;
axis image; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
yticks([.5:20:nPCs]);
xticks([.5 15.5 30.5]); xtickangle(0);
title(sprintf('%s-%s-v41',fov,plotType));

subplot(4,6,[1:6:13]+2);
imagesc(kernelMatv4n',[0 .5]); colormap viridis;
axis image; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
yticks([.5:20:nPCs]);
xticks([.5 15.5 30.5]); xtickangle(0);
title(sprintf('%s-%s-v4n',fov,plotType));

% Mean traces
subplot(4,6,19);
boundedline(xmap,mean(kernelMatv1,2,'omitnan'),std(kernelMatv1,[],2,'omitnan')/sqrt(sum(~isnan(kernelMatv1(1,:)))),'cmap',[0 0 0]);
ylim([-.25 1]);
yticks([-.25:.25:1]);
xlim([0 size(kernelMatv1,1)]+.5);
xticks([.5 15.5 30.5]); xtickangle(0);
axis square; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');

subplot(4,6,20);
boundedline(xmap,mean(kernelMatv41,2,'omitnan'),std(kernelMatv41,[],2,'omitnan')/sqrt(sum(~isnan(kernelMatv41(1,:)))),'cmap',[0 0 0]);
ylim([-.25 1]);
yticks([-.25:.25:1]);
xlim([0 size(kernelMatv41,1)]+.5);
xticks([.5 15.5 30.5]); xtickangle(0);
axis square; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');

subplot(4,6,21);
boundedline(xmap,mean(kernelMatv4n,2,'omitnan'),std(kernelMatv4n,[],2,'omitnan')/sqrt(sum(~isnan(kernelMatv4n(1,:)))),'cmap',[0 0 0]);
ylim([-.25 1]);
yticks([-.25:.25:1]);
xlim([0 size(kernelMatv4n,1)]+.5);
xticks([.5 15.5 30.5]); xtickangle(0);
axis square; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');

% Big correlation matrix
subplot(2,4,3); imagesc(Rmat,[0 1]);
axis square
box off
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([.5:nPCs:3*nPCs+.5]);
yticks([.5:nPCs:3*nPCs+.5]);
xline(nPCs+.5,'--w');
xline(2*nPCs+.5,'--w');
yline(nPCs+.5,'--w');
yline(2*nPCs+.5,'--w');

% Histgrams of indvidual comparisons
subplot(2,6,10);
h1 = histogram(Rv1v41Matched(~isnan(Rv1v41Matched)),[-1:.05:1],'Normalization','probability');
axis square
box off
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('v1-v41 Rmatched')
xticks([-1:.5:1]); xtickangle(0)
xline(median(Rv1v41Matched,'omitnan'),'--k');
xlim([-1.025 1.025]);
ylim([0 .3]);

subplot(2,6,11);
h2 = histogram(Rv1v4nMatched(~isnan(Rv1v4nMatched)),[-1:.05:1],'Normalization','probability');
axis square
box off
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('v1-v4n Rmatched')
xticks([-1:.5:1]); xtickangle(0)
xline(median(Rv1v4nMatched,'omitnan'),'--k');
xlim([-1.025 1.025]);
ylim([0 .3]);

subplot(2,6,12);
h1 = histogram(Rv41v4nMatched(~isnan(Rv41v4nMatched)),[-1:.05:1],'Normalization','probability');
axis square
box off
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('v41-v4n Rmatched')
xline(median(Rv41v4nMatched,'omitnan'),'--k');
xticks([-1:.5:1]); xtickangle(0)
xlim([-1.025 1.025]);
ylim([0 .3]);

% Cumulative histograms
v1v41Cum = sort(Rv1v41Matched(~isnan(Rv1v41Matched)),'ascend');
v1v4nCum = sort(Rv1v4nMatched(~isnan(Rv1v4nMatched)),'ascend');
v41v4nCum = sort(Rv41v4nMatched(~isnan(Rv41v4nMatched)),'ascend');

subplot(2,4,4);
plot([-1 v1v41Cum(1) v1v41Cum 1],[0,0:numel(v1v41Cum),numel(v1v41Cum)]/numel(v1v41Cum),'r','LineWidth',1.5); hold on;
plot([-1 v1v4nCum(1) v1v4nCum 1],[0,0:numel(v1v4nCum),numel(v1v4nCum)]/numel(v1v4nCum),'Color',[.4 .4 .4],'LineWidth',1.5);
plot([-1 v41v4nCum(1) v41v4nCum 1],[0,0:numel(v41v4nCum),numel(v41v4nCum)]/numel(v41v4nCum),'k','LineWidth',1.5);
axis square
box off
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1.025 1.025 0 1.025]);
legend('naive-learning','naive-expert','learning-expert','Location','northwest');

%% Crunch stats
% [~,p] = kstest2(v1v41Cum,v1v4nCum);
% ksOut(1) = p;
% [~,p] = kstest2(v1v41Cum,v41v4nCum);
% ksOut(2) = p;
% [~,p] = kstest2(v1v4nCum,v41v4nCum);
% ksOut(3) = p;

Rmatched0 = [Rv1v41Matched; Rv1v4nMatched; Rv41v4nMatched]';
Rmatched_iGood = sum(isnan(Rmatched0),2) == 0;
Rmatched = Rmatched0(Rmatched_iGood,:);

[p,tbl,stats] = friedman(Rmatched);
c = multcompare(stats,'CType', 'bonferroni');

statsOut = v2struct(p,tbl,stats,c);

end
