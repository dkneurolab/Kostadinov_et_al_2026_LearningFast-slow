function fig4_kernelRvsDist_v4only(seshv4, paperFold, savebool)


%% Load GLMv144:

glmFold = fullfile(paperFold,'Figures','Fig4','Sub5_fitKernels');

GLMv144 = load(fullfile(glmFold,'kernelsGLMv144.mat'));
GLMv144 = GLMv144.GLMv144;

%% Loop through and get matching xCoordPx

for isesh = 1:numel(GLMv144)
    xCtrPx = [];
    for jseshv4 = 1:numel(seshv4)
        if (strcmpi(GLMv144(isesh).name,seshv4(jseshv4).name) && GLMv144(isesh).ilobv && strcmpi(seshv4(jseshv4).fov,'lobv')) || ...
                strcmpi(GLMv144(isesh).name,seshv4(jseshv4).name) && GLMv144(isesh).isim2 && strcmpi(seshv4(jseshv4).fov,'sim2')
            xCtrPx = seshv4(jseshv4).xCtrsPx;
        end
    end
    iPC_v4n_v41_all = cat(2,GLMv144(isesh).GLMmatchv4n_v41.iPC_v4);

    GLMv144(isesh).movMatchInd = [];
    GLMv144(isesh).movMatch_xCtr = [];
    GLMv144(isesh).movMatch_Dist = [];
    GLMv144(isesh).movMatch_Similarity = [];
    GLMv144(isesh).rewMatchInd = [];
    GLMv144(isesh).rewMatch_xCtr = [];
    GLMv144(isesh).rewMatch_Dist = [];
    GLMv144(isesh).rewMatch_Similarity = [];

    % Deal with movement first:
    plotType = 'mov';
    kernelMatv41 = cat(2,GLMv144(isesh).(sprintf('%sKernel_v41',plotType)));
    iGood = ~isnan(kernelMatv41(1,:));
    if sum(iGood) > 0
        % Compute correlations again (same math as in fig4_fitKernels_plot.m)
        kernelMatv41 = kernelMatv41(:,iGood);
        kernelMatv4n = cat(2,GLMv144(isesh).(sprintf('%sKernel_v4n',plotType)));
        kernelMatv4n = kernelMatv4n(:,iGood);
        nPCs = size(kernelMatv41,2);
        Rmat = corrcoef([kernelMatv41 kernelMatv4n]);
        Rmat_compare = Rmat(1:nPCs,nPCs+1:end);
        Rmat_matched0 = Rmat_compare(1:nPCs+1:end)';

        % Compute cell IDs and locations for same cells
        sigFits = cat(2,GLMv144(isesh).GLMmatchv4n_v41.sigBool_shufMov);
        xCtrInds0 = iPC_v4n_v41_all(~isnan(sigFits));
        xCtrInds0 = xCtrInds0(iGood);
        xCtr0 = xCtrPx(xCtrInds0)*670/512;

        % Check again for nans and make final all to all computations
        matchedGood = ~isnan(Rmat_matched0);
        Rmat_matched = Rmat_matched0(matchedGood);
        xCtrInds = xCtrInds0(matchedGood);
        xCtrMatch = xCtr0(matchedGood);

        matchDist0 = xCtrMatch - xCtrMatch.';
        matchSimilarity0 = Rmat_matched - Rmat_matched.';

        matchDist = abs(matchDist0(triu(true(size(matchDist0)),1)));
        matchSimilarity = abs(matchSimilarity0(triu(true(size(matchSimilarity0)),1)));

        GLMv144(isesh).movMatchInd = xCtrInds';
        GLMv144(isesh).movMatch_xCtr = xCtrMatch;
        GLMv144(isesh).movMatch_Dist = matchDist;
        GLMv144(isesh).movMatch_Similarity = matchSimilarity;
    end

    % Deal with reward second
    plotType = 'rew';
    kernelMatv41 = cat(2,GLMv144(isesh).(sprintf('%sKernel_v41',plotType)));
    iGood = ~isnan(kernelMatv41(1,:));
    if sum(iGood) > 0
        % Compute correlations again (same math as in fig4_fitKernels_plot.m)
        kernelMatv41 = kernelMatv41(:,iGood);
        kernelMatv4n = cat(2,GLMv144(isesh).(sprintf('%sKernel_v4n',plotType)));
        kernelMatv4n = kernelMatv4n(:,iGood);
        nPCs = size(kernelMatv41,2);
        Rmat = corrcoef([kernelMatv41 kernelMatv4n]);
        Rmat_compare = Rmat(1:nPCs,nPCs+1:end);
        Rmat_matched0 = Rmat_compare(1:nPCs+1:end)';

        % Compute cell IDs and locations for same cells
        sigFits = cat(2,GLMv144(isesh).GLMmatchv4n_v41.sigBool_shufRew);
        xCtrInds0 = iPC_v4n_v41_all(~isnan(sigFits));
        xCtrInds0 = xCtrInds0(iGood);
        xCtr0 = xCtrPx(xCtrInds0)*670/512;

        % Check again for nans and make final all to all computations
        matchedGood = ~isnan(Rmat_matched0);
        Rmat_matched = Rmat_matched0(matchedGood);
        xCtrInds = xCtrInds0(matchedGood);
        xCtrMatch = xCtr0(matchedGood);

        matchDist0 = xCtrMatch - xCtrMatch.';
        matchSimilarity0 = Rmat_matched - Rmat_matched.';

        matchDist = abs(matchDist0(triu(true(size(matchDist0)),1)));
        matchSimilarity = abs(matchSimilarity0(triu(true(size(matchSimilarity0)),1)));

        GLMv144(isesh).rewMatchInd = xCtrInds';
        GLMv144(isesh).rewMatch_xCtr = xCtrMatch;
        GLMv144(isesh).rewMatch_Dist = matchDist;
        GLMv144(isesh).rewMatch_Similarity = matchSimilarity;
    end

end

%% Concatenate and do stats:
movDistAll = cat(1,GLMv144.movMatch_Dist);
movSimAll = cat(1,GLMv144.movMatch_Similarity);
rewDistAll = cat(1,GLMv144.rewMatch_Dist);
rewSimAll = cat(1,GLMv144.rewMatch_Similarity);

f1 = figure;
subplot(2,2,1);
scatter(movDistAll,movSimAll, 'filled', MarkerEdgeColor='none',MarkerFaceColor = 'b', MarkerFaceAlpha = 0.3); axis square; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlim([-25 650]); ylim([-0.05 1.55])
ylabel('Difference in v41-v4n corr.')
xlabel('Intercellular distance (um)')
title('Slow learning: movement')
subplot(2,2,2);
scatter(rewDistAll,rewSimAll, 'filled', MarkerEdgeColor='none',MarkerFaceColor = 'b', MarkerFaceAlpha = 0.3); axis square; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlim([-25 650]); ylim([-0.05 1.55])
ylabel('Difference in v41-v4n corr.')
xlabel('Intercellular distance (um)')
title('Slow learning: reward')

iMov_under100 = find(movDistAll < 100);
iMov_over100 = find(movDistAll >= 100);
iRew_under100 = find(rewDistAll < 100);
iRew_over100 = find(rewDistAll >= 100);

subplot(2,2,3); hold on;
xplot = randn(numel(iMov_under100),1)*.2;
xplot(abs(xplot) > 0.5) = 0;
scatter(xplot+1,movSimAll(iMov_under100),'filled', MarkerEdgeColor='none',MarkerFaceColor = 'k', MarkerFaceAlpha = 0.3); axis square; box off;
xplot = randn(numel(iMov_over100),1)*.2;
xplot(abs(xplot) > 0.5) = 0;
scatter(xplot+3,movSimAll(iMov_over100),'filled', MarkerEdgeColor='none',MarkerFaceColor = [0.5 .1 .1], MarkerFaceAlpha = 0.3); axis square; box off;
xlim([0 4]); ylim([-0.05 1.55])
axis square; box off
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([1,3]); xticklabels({'<100','>=100'})
xlabel('Intercellular distance (um)')
ylabel('Difference in v41-v4n corr.')
a = ranksum(movSimAll(iMov_under100),movSimAll(iMov_over100));
text(2,1.5,sprintf('p = %.3f',a));

subplot(2,2,4); hold on;
xplot = randn(numel(iRew_under100),1)*.2;
xplot(abs(xplot) > 0.5) = 0;
scatter(xplot+1,rewSimAll(iRew_under100),'filled', MarkerEdgeColor='none',MarkerFaceColor = 'k', MarkerFaceAlpha = 0.3); axis square; box off;
xplot = randn(numel(iRew_over100),1)*.2;
xplot(abs(xplot) > 0.5) = 0;
scatter(xplot+3,rewSimAll(iRew_over100),'filled', MarkerEdgeColor='none',MarkerFaceColor = [0.5 .1 .1], MarkerFaceAlpha = 0.3); axis square; box off;
xlim([0 4]); ylim([-0.05 1.55])
axis square; box off
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([1,3]); xticklabels({'<100','>=100'})
xlabel('Intercellular distance (um)')
ylabel('Difference in v41-v4n corr.')
a = ranksum(rewSimAll(iRew_under100),rewSimAll(iRew_over100));
text(2,1.5,sprintf('p = %.3f',a));

suptitle('Distance-dependence of kernel correlations: Learning-expert');

%% Concatenate and do stats: lobule V
ilobv = cat(1,GLMv144.ilobv);
movDistAll_lobv = cat(1,GLMv144(ilobv).movMatch_Dist);
movSimAll_lobv = cat(1,GLMv144(ilobv).movMatch_Similarity);
rewDistAll_lobv = cat(1,GLMv144(ilobv).rewMatch_Dist);
rewSimAll_lobv = cat(1,GLMv144(ilobv).rewMatch_Similarity);

f2 = figure;
subplot(2,2,1);
scatter(movDistAll_lobv,movSimAll_lobv, 'filled', MarkerEdgeColor='none',MarkerFaceColor = 'b', MarkerFaceAlpha = 0.3); axis square; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlim([-25 650]); ylim([-0.05 1.55])
ylabel('Difference in v41-v4n corr.')
xlabel('Intercellular distance (um)')
title('Slow learning: movement (lobV)')
subplot(2,2,2);
scatter(rewDistAll_lobv,rewSimAll_lobv, 'filled', MarkerEdgeColor='none',MarkerFaceColor = 'b', MarkerFaceAlpha = 0.3); axis square; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlim([-25 650]); ylim([-0.05 1.55])
ylabel('Difference in v41-v4n corr.')
xlabel('Intercellular distance (um)')
title('Slow learning: reward (lobV)')

iMov_under100_lobv = find(movDistAll_lobv < 100);
iMov_over100_lobv = find(movDistAll_lobv >= 100);
iRew_under100_lobv = find(rewDistAll_lobv < 100);
iRew_over100_lobv = find(rewDistAll_lobv >= 100);

subplot(2,2,3); hold on;
xplot = randn(numel(iMov_under100_lobv),1)*.2;
xplot(abs(xplot) > 0.5) = 0;
scatter(xplot+1,movSimAll_lobv(iMov_under100_lobv),'filled', MarkerEdgeColor='none',MarkerFaceColor = 'k', MarkerFaceAlpha = 0.3); axis square; box off;
xplot = randn(numel(iMov_over100_lobv),1)*.2;
xplot(abs(xplot) > 0.5) = 0;
scatter(xplot+3,movSimAll_lobv(iMov_over100_lobv),'filled', MarkerEdgeColor='none',MarkerFaceColor = [0.5 .1 .1], MarkerFaceAlpha = 0.3); axis square; box off;
xlim([0 4]); ylim([-0.05 1.55])
axis square; box off
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([1,3]); xticklabels({'<100','>=100'})
xlabel('Intercellular distance (um)')
ylabel('Difference in v41-v4n corr.')
a = ranksum(movSimAll_lobv(iMov_under100_lobv),movSimAll_lobv(iMov_over100_lobv));
text(2,1.5,sprintf('p = %.3f',a));

subplot(2,2,4); hold on;
xplot = randn(numel(iRew_under100_lobv),1)*.2;
xplot(abs(xplot) > 0.5) = 0;
scatter(xplot+1,rewSimAll_lobv(iRew_under100_lobv),'filled', MarkerEdgeColor='none',MarkerFaceColor = 'k', MarkerFaceAlpha = 0.3); axis square; box off;
xplot = randn(numel(iRew_over100_lobv),1)*.2;
xplot(abs(xplot) > 0.5) = 0;
scatter(xplot+3,rewSimAll_lobv(iRew_over100_lobv),'filled', MarkerEdgeColor='none',MarkerFaceColor = [0.5 .1 .1], MarkerFaceAlpha = 0.3); axis square; box off;
xlim([0 4]); ylim([-0.05 1.55])
axis square; box off
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([1,3]); xticklabels({'<100','>=100'})
xlabel('Intercellular distance (um)')
ylabel('Difference in v41-v4n corr.')
a = ranksum(rewSimAll_lobv(iRew_under100_lobv),rewSimAll_lobv(iRew_over100_lobv));
text(2,1.5,sprintf('p = %.3f',a));

suptitle('Distance-dependence of kernel correlations: Learning-expert');

%% Concatenate and do stats:
isim2 = cat(1,GLMv144.isim2);
movDistAll_sim2 = cat(1,GLMv144(isim2).movMatch_Dist);
movSimAll_sim2 = cat(1,GLMv144(isim2).movMatch_Similarity);
rewDistAll_sim2 = cat(1,GLMv144(isim2).rewMatch_Dist);
rewSimAll_sim2 = cat(1,GLMv144(isim2).rewMatch_Similarity);

f3 = figure;
subplot(2,2,1);
scatter(movDistAll_sim2,movSimAll_sim2, 'filled', MarkerEdgeColor='none',MarkerFaceColor = 'b', MarkerFaceAlpha = 0.3); axis square; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlim([-25 650]); ylim([-0.05 1.55])
ylabel('Difference in v41-v4n corr.')
xlabel('Intercellular distance (um)')
title('Slow learning: movement (sim2)')
subplot(2,2,2);
scatter(rewDistAll_sim2,rewSimAll_sim2, 'filled', MarkerEdgeColor='none',MarkerFaceColor = 'b', MarkerFaceAlpha = 0.3); axis square; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlim([-25 650]); ylim([-0.05 1.55])
ylabel('Difference in v41-v4n corr.')
xlabel('Intercellular distance (um)')
title('Slow learning: reward (sim2)')

iMov_under100_sim2 = find(movDistAll_sim2 < 100);
iMov_over100_sim2 = find(movDistAll_sim2 >= 100);
iRew_under100_sim2 = find(rewDistAll_sim2 < 100);
iRew_over100_sim2 = find(rewDistAll_sim2 >= 100);

subplot(2,2,3); hold on;
xplot = randn(numel(iMov_under100_sim2),1)*.2;
xplot(abs(xplot) > 0.5) = 0;
scatter(xplot+1,movSimAll_sim2(iMov_under100_sim2),'filled', MarkerEdgeColor='none',MarkerFaceColor = 'k', MarkerFaceAlpha = 0.3); axis square; box off;
xplot = randn(numel(iMov_over100_sim2),1)*.2;
xplot(abs(xplot) > 0.5) = 0;
scatter(xplot+3,movSimAll_sim2(iMov_over100_sim2),'filled', MarkerEdgeColor='none',MarkerFaceColor = [0.5 .1 .1], MarkerFaceAlpha = 0.3); axis square; box off;
xlim([0 4]); ylim([-0.05 1.55])
axis square; box off
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([1,3]); xticklabels({'<100','>=100'})
xlabel('Intercellular distance (um)')
ylabel('Difference in v41-v4n corr.')
a = ranksum(movSimAll_sim2(iMov_under100_sim2),movSimAll_sim2(iMov_over100_sim2));
text(2,1.5,sprintf('p = %.3f',a));

subplot(2,2,4); hold on;
xplot = randn(numel(iRew_under100_sim2),1)*.2;
xplot(abs(xplot) > 0.5) = 0;
scatter(xplot+1,rewSimAll_sim2(iRew_under100_sim2),'filled', MarkerEdgeColor='none',MarkerFaceColor = 'k', MarkerFaceAlpha = 0.3); axis square; box off;
xplot = randn(numel(iRew_over100_sim2),1)*.2;
xplot(abs(xplot) > 0.5) = 0;
scatter(xplot+3,rewSimAll_sim2(iRew_over100_sim2),'filled', MarkerEdgeColor='none',MarkerFaceColor = [0.5 .1 .1], MarkerFaceAlpha = 0.3); axis square; box off;
xlim([0 4]); ylim([-0.05 1.55])
axis square; box off
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([1,3]); xticklabels({'<100','>=100'})
xlabel('Intercellular distance (um)')
ylabel('Difference in v41-v4n corr.')
a = ranksum(rewSimAll_sim2(iRew_under100_sim2),rewSimAll_sim2(iRew_over100_sim2));
text(2,1.5,sprintf('p = %.3f',a));

suptitle('Distance-dependence of kernel correlations: Learning-expert');

%% Summary bar graphs only

f4 = figure;
subplot(2,3,1)
e1 = errorbar([1,3],[mean(movSimAll(iMov_under100)), mean(movSimAll(iMov_over100))], [std(movSimAll(iMov_under100))/sqrt(numel(iMov_under100)), std(movSimAll(iMov_over100))/sqrt(numel(iMov_over100))]);
% e1 = errorbar([1,3],[mean(movSimAll(iMov_under100)), mean(movSimAll(iMov_over100))], [std(movSimAll(iMov_under100)), std(movSimAll(iMov_over100))]);
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
hold on; axis square; box off
scatter([1,3],[mean(movSimAll(iMov_under100)), mean(movSimAll(iMov_over100))],60,'MarkerEdgeColor','none','MarkerFaceColor','k');
xlim([0 4]); ylim([0 0.6])
% xlim([0 4]); ylim([-0.05 1.05])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([1,3]); xticklabels({'<100','>=100'})
xlabel('Intercellular distance (um)')
ylabel('Difference in v1-v4n corr.')
title('Combined FOVs, movement')
a = ranksum(movSimAll(iMov_under100),movSimAll(iMov_over100));
text(2,.5,sprintf('p = %.3f',a));

subplot(2,3,4)
e1 = errorbar([1,3],[mean(rewSimAll(iRew_under100)), mean(rewSimAll(iRew_over100))], [std(rewSimAll(iRew_under100))/sqrt(numel(iRew_under100)), std(rewSimAll(iRew_over100))/sqrt(numel(iRew_over100))]);
% e1 = errorbar([1,3],[mean(rewSimAll(iRew_under100)), mean(rewSimAll(iRew_over100))], [std(rewSimAll(iRew_under100)), std(rewSimAll(iRew_over100))]);
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
hold on; axis square; box off
scatter([1,3],[mean(rewSimAll(iRew_under100)), mean(rewSimAll(iRew_over100))],60,'MarkerEdgeColor','none','MarkerFaceColor','k');
xlim([0 4]); ylim([0.0 0.6])
% xlim([0 4]); ylim([-0.05 1.05])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([1,3]); xticklabels({'<100','>=100'})
xlabel('Intercellular distance (um)')
ylabel('Difference in v1-v4n corr.')
title('Combined FOVs, reward')
a = ranksum(rewSimAll(iRew_under100),rewSimAll(iRew_over100));
text(2,.5,sprintf('p = %.3f',a));

subplot(2,3,2)
e1 = errorbar([1,3],[mean(movSimAll_lobv(iMov_under100_lobv)), mean(movSimAll_lobv(iMov_over100_lobv))], [std(movSimAll_lobv(iMov_under100_lobv))/sqrt(numel(iMov_under100_lobv)), std(movSimAll_lobv(iMov_over100_lobv))/sqrt(numel(iMov_over100_lobv))]);
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
hold on; axis square; box off
scatter([1,3],[mean(movSimAll_lobv(iMov_under100_lobv)), mean(movSimAll_lobv(iMov_over100_lobv))],60,'MarkerEdgeColor','none','MarkerFaceColor','k');
xlim([0 4]); ylim([0.0 0.6])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([1,3]); xticklabels({'<100','>=100'})
xlabel('Intercellular distance (um)')
ylabel('Difference in v1-v4n corr.')
title('Lobule V, movement')
a = ranksum(movSimAll_lobv(iMov_under100_lobv),movSimAll_lobv(iMov_over100_lobv));
text(2,.5,sprintf('p = %.3f',a));

subplot(2,3,5)
e1 = errorbar([1,3],[mean(rewSimAll_lobv(iRew_under100_lobv)), mean(rewSimAll_lobv(iRew_over100_lobv))], [std(rewSimAll_lobv(iRew_under100_lobv))/sqrt(numel(iRew_over100_lobv)), std(rewSimAll_lobv(iRew_over100_lobv))/sqrt(numel(iRew_over100_lobv))]);
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
hold on; axis square; box off
scatter([1,3],[mean(rewSimAll_lobv(iRew_under100_lobv)), mean(rewSimAll_lobv(iRew_over100_lobv))],60,'MarkerEdgeColor','none','MarkerFaceColor','k');
xlim([0 4]); ylim([0.0 0.6])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([1,3]); xticklabels({'<100','>=100'})
xlabel('Intercellular distance (um)')
ylabel('Difference in v1-v4n corr.')
title('Lobule V, reward')
a = ranksum(rewSimAll_lobv(iRew_under100_lobv),rewSimAll_lobv(iRew_over100_lobv));
text(2,.5,sprintf('p = %.3f',a));

subplot(2,3,3)
e1 = errorbar([1,3],[mean(movSimAll_sim2(iMov_under100_sim2)), mean(movSimAll_sim2(iMov_over100_sim2))], [std(movSimAll_sim2(iMov_under100_sim2))/sqrt(numel(iMov_under100_sim2)), std(movSimAll_sim2(iMov_over100_sim2))/sqrt(numel(iMov_over100_sim2))]);
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
hold on; axis square; box off
scatter([1,3],[mean(movSimAll_sim2(iMov_under100_sim2)), mean(movSimAll_sim2(iMov_over100_sim2))],60,'MarkerEdgeColor','none','MarkerFaceColor','k');
xlim([0 4]); ylim([0.0 0.6])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([1,3]); xticklabels({'<100','>=100'})
xlabel('Intercellular distance (um)')
ylabel('Difference in v1-v4n corr.')
title('Simplex, movement')
a = ranksum(movSimAll_sim2(iMov_under100_sim2),movSimAll_sim2(iMov_over100_sim2));
text(2,.5,sprintf('p = %.3f',a));

subplot(2,3,6)
e1 = errorbar([1,3],[mean(rewSimAll_sim2(iRew_under100_sim2)), mean(rewSimAll_sim2(iRew_over100_sim2))], [std(rewSimAll_sim2(iRew_under100_sim2))/sqrt(numel(iRew_under100_sim2)), std(rewSimAll_sim2(iRew_over100_sim2))/sqrt(numel(iRew_over100_sim2))]);
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
hold on; axis square; box off
scatter([1,3],[mean(rewSimAll_sim2(iRew_under100_sim2)), mean(rewSimAll_sim2(iRew_over100_sim2))],60,'MarkerEdgeColor','none','MarkerFaceColor','k');
xlim([0 4]); ylim([0.0 0.6])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([1,3]); xticklabels({'<100','>=100'})
xlabel('Intercellular distance (um)')
ylabel('Difference in v1-v4n corr.')
title('Simplex, reward')
a = ranksum(rewSimAll_sim2(iRew_under100_sim2),rewSimAll_sim2(iRew_over100_sim2));
text(2,.5,sprintf('p = %.3f',a));

suptitle('Distance-dependence of kernel correlations: Learning-expert');

%% Save maybe
if savebool
    figname = 'kernelRvsDist_v41_v4n_all';
    savefig(f1, fullfile(glmFold,[figname,'.fig']));
    saveas(f1,fullfile(glmFold,[figname,'.png']));
    print(f1,fullfile(glmFold,[figname,'.eps']), '-depsc', '-vector');

    figname = 'kernelRvsDist_v41_v4n_lobv';
    savefig(f2, fullfile(glmFold,[figname,'.fig']));
    saveas(f2,fullfile(glmFold,[figname,'.png']));
    print(f2,fullfile(glmFold,[figname,'.eps']), '-depsc', '-vector');

    figname = 'kernelRvsDist_v41_v4n_sim2';
    savefig(f3, fullfile(glmFold,[figname,'.fig']));
    saveas(f3,fullfile(glmFold,[figname,'.png']));
    print(f3,fullfile(glmFold,[figname,'.eps']), '-depsc', '-vector');

    figname = 'kernelRvsDist_v41_v4n_summary';
    savefig(f4, fullfile(glmFold,[figname,'.fig']));
    saveas(f4,fullfile(glmFold,[figname,'.png']));
    print(f4,fullfile(glmFold,[figname,'.eps']), '-depsc', '-vector');
    
    GLMdistance_v41_v4 = rmfield(GLMv144,{'GLMmatchv4n_v41', 'GLMmatchv4n_v41'});

    save(fullfile(glmFold,'GLMdistance_v41_v4.mat'),'GLMdistance_v41_v4');

end