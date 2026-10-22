function [roiStatsAll, f1] = fig2_Supp_roiCorrs_summary_rev1(roiStatsAll)
% Make summary plots for ROI matching

%% Loop through and tabulate stats

for iroi = 1:numel(roiStatsAll)
    roiStatsAll(iroi).delDays = abs(days(datetime(roiStatsAll(iroi).date2, 'InputFormat', 'yyyy-MM-dd') - datetime(roiStatsAll(iroi).date1, 'InputFormat', 'yyyy-MM-dd')));
    
    % rmean
    roiStatsAll(iroi).r_real_mu = mean(roiStatsAll(iroi).r_real,'omitnan');
    roiStatsAll(iroi).r_real_sem = std(roiStatsAll(iroi).r_real,'omitnan') / sqrt(sum(~isnan(roiStatsAll(iroi).r_real)));
    roiStatsAll(iroi).r_next_mu = mean(roiStatsAll(iroi).r_next,'omitnan');
    roiStatsAll(iroi).r_next_sem = std(roiStatsAll(iroi).r_next,'omitnan') / sqrt(sum(~isnan(roiStatsAll(iroi).r_next)));
    
    % perc overlap
    roiStatsAll(iroi).foot_real_mu = mean(roiStatsAll(iroi).percFoot_real,'omitnan');
    roiStatsAll(iroi).foot_real_sem = std(roiStatsAll(iroi).percFoot_real,'omitnan') / sqrt(sum(~isnan(roiStatsAll(iroi).percFoot_real)));
    roiStatsAll(iroi).foot_next_mu = mean(roiStatsAll(iroi).percFoot_next,'omitnan');
    roiStatsAll(iroi).foot_next_sem = std(roiStatsAll(iroi).percFoot_next,'omitnan') / sqrt(sum(~isnan(roiStatsAll(iroi).percFoot_next)));

    % ctr dist
    roiStatsAll(iroi).ctr_real_mu = mean(roiStatsAll(iroi).ctrDist_real,'omitnan');
    roiStatsAll(iroi).ctr_real_sem = std(roiStatsAll(iroi).ctrDist_real,'omitnan') / sqrt(sum(~isnan(roiStatsAll(iroi).ctrDist_real)));
    roiStatsAll(iroi).ctr_next_mu = mean(roiStatsAll(iroi).ctrDist_next,'omitnan');
    roiStatsAll(iroi).ctr_next_sem = std(roiStatsAll(iroi).ctrDist_next,'omitnan') / sqrt(sum(~isnan(roiStatsAll(iroi).ctrDist_next)));

end

%% Plot stuff

delDays = cat(1,roiStatsAll.delDays);
pairplotX = [ones(numel(roiStatsAll),1) 2*ones(numel(roiStatsAll),1)];

f1 = figure;
% ROI corrs
subplot(2,3,1); hold on;
yData1 = cat(1, roiStatsAll.r_real_mu);
yError1 = cat(1, roiStatsAll.r_real_sem);
yData2 = cat(1, roiStatsAll.r_next_mu);
yError2 = cat(1, roiStatsAll.r_next_sem);
oData1 = regress_hd(delDays,yData1,1);
oData2 = regress_hd(delDays,yData2,1);
s1 = scatter(delDays, yData1, 36, 'k','filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5);
e1 = errorbar(delDays, yData1, yError1,'k', 'LineStyle','none', 'CapSize',0, 'LineWidth',2);
s2 = scatter(delDays, yData2, 36, 'b','filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5);
e2 = errorbar(delDays, yData2, yError2,'k', 'LineStyle','none', 'CapSize',0, 'LineWidth',2);
h1 = boundedline(oData1.x_fit,oData1.y_fit,oData1.delta,'k','alpha');
h2 = boundedline(oData2.x_fit,oData2.y_fit,oData2.delta,'b','alpha');
axis square; box off;
xlim([-1.5 75]); ylim([-.02 1]);
xticks([0:15:75]); yticks([0:0.2:1]);
ylabel('ROI footprint correlation');
xlabel('Days between sessions');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
text(45,0.2,{sprintf('rsqReal = %0.3f',oData1.rsq); sprintf('pReal = %0.3f',oData1.corr_p); sprintf('rsqNext = %0.3f',oData2.rsq); sprintf('pNext = %0.3f',oData2.corr_p)});

subplot(2,3,4);
plot(pairplotX', [yData1, yData2]', 'k', LineWidth = 2);
axis square; box off;
xlim([0.5 2.5]); ylim([-.02 1]);
xticks([1,2]); yticks([0:0.2:1]); xticklabels({'Real','Neighbour'})
ylabel('ROI footprint correlation');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
[a,b] = signrank(yData1,yData2);
text(1.5, 0.2, sprintf('p = %0.4f', a))

% Footprint
subplot(2,3,2); hold on;
yData1 = cat(1, roiStatsAll.foot_real_mu);
yError1 = cat(1, roiStatsAll.foot_real_sem);
yData2 = cat(1, roiStatsAll.foot_next_mu);
yError2 = cat(1, roiStatsAll.foot_next_sem);
oData1 = regress_hd(delDays,yData1,1);
oData2 = regress_hd(delDays,yData2,1);
s1 = scatter(delDays, yData1, 36, 'k','filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5);
e1 = errorbar(delDays, yData1, yError1,'k', 'LineStyle','none', 'CapSize',0, 'LineWidth',2);
s2 = scatter(delDays, yData2, 36, 'b','filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5);
e2 = errorbar(delDays, yData2, yError2,'k', 'LineStyle','none', 'CapSize',0, 'LineWidth',2);
h1 = boundedline(oData1.x_fit,oData1.y_fit,oData1.delta,'k','alpha');
h2 = boundedline(oData2.x_fit,oData2.y_fit,oData2.delta,'b','alpha');
axis square; box off;
xlim([-1.5 75]); ylim([-.02 1]);
xticks([0:15:75]); yticks([0:0.2:1]);
ylabel({'Pixel overlap per', 'matched PC dendrite (%)'});
xlabel('Days between sessions');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
text(45,.2,{sprintf('rsqReal = %0.3f',oData1.rsq); sprintf('pReal = %0.3f',oData1.corr_p); sprintf('rsqNext = %0.3f',oData2.rsq); sprintf('pNext = %0.3f',oData2.corr_p)});

subplot(2,3,5);
plot(pairplotX', [yData1, yData2]', 'k', LineWidth = 2);
axis square; box off;
xlim([0.5 2.5]); ylim([-.02 1]);
xticks([1,2]); yticks([0:0.2:1]); xticklabels({'Real','Neighbour'})
ylabel({'Pixel overlap per', 'matched PC dendrite (%)'});
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
[a,b] = signrank(yData1,yData2);
text(1.5, 0.2, sprintf('p = %0.4f', a))

% Inter-centroid distance
subplot(2,3,3); hold on;
yData1 = cat(1, roiStatsAll.ctr_real_mu);
yError1 = cat(1, roiStatsAll.ctr_real_sem);
yData2 = cat(1, roiStatsAll.ctr_next_mu);
yError2 = cat(1, roiStatsAll.ctr_next_sem);
oData1 = regress_hd(delDays,yData1,1);
oData2 = regress_hd(delDays,yData2,1);
s1 = scatter(delDays, yData1, 36, 'k','filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5);
e1 = errorbar(delDays, yData1, yError1,'k', 'LineStyle','none', 'CapSize',0, 'LineWidth',2);
s2 = scatter(delDays, yData2, 36, 'b','filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5);
e2 = errorbar(delDays, yData2, yError2,'k', 'LineStyle','none', 'CapSize',0, 'LineWidth',2);
h1 = boundedline(oData1.x_fit,oData1.y_fit,oData1.delta,'k','alpha');
h2 = boundedline(oData2.x_fit,oData2.y_fit,oData2.delta,'b','alpha');
axis square; box off;
xlim([-1.5 75]); ylim([-.02 1]*250);
xticks([0:15:75]); yticks([0:0.2:1]*250);
ylabel('Inter-centroid distance (um)');
xlabel('Days between sessions');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
text(60,125,{sprintf('rsqReal = %0.3f',oData1.rsq); sprintf('pReal = %0.3f',oData1.corr_p); sprintf('rsqNext = %0.3f',oData2.rsq); sprintf('pNext = %0.3f',oData2.corr_p)});

subplot(2,3,6);
plot(pairplotX', [yData1, yData2]', 'k', LineWidth = 2);
axis square; box off;
xlim([0.5 2.5]); ylim([-.02 1]*250);
xticks([1,2]); yticks([0:0.2:1]*250); xticklabels({'Real','Neighbour'})
ylabel('Inter-centroid distance (um)');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
[a,b] = signrank(yData1,yData2);
text(2, 50, sprintf('p = %0.4f', a))


end