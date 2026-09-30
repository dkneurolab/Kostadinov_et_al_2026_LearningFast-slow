function fig1_v1to4learning_rev1(seshV1, seshV4,nSesh,dataFold,paperFold,savebool)
%% Note that pretty much everything is hard-coded with the expectation that nSesh == 6

% Set directories for behavioural data and output figures
figFold = fullfile(paperFold,'Figures','Fig1','Sub1_v1Learning_rev1');



%% Loop through
mouseNames = unique(cat(1,seshV4.name),'rows');
mouseStr4 = cat(2,seshV4.name);
mouseStr1 = cat(2,seshV1.name);

learnMat = struct;
for iMouse = 1:size(mouseNames,1)
    learnMat(iMouse).name = mouseNames(iMouse,:);
    learnMat(iMouse).iSeshv1 = (regexp(mouseStr1,mouseNames(iMouse,:))+4)/5;
    learnMat(iMouse).iSeshv4 = (regexp(mouseStr4,mouseNames(iMouse,:))+4)/5;
    seshTemp = seshV4(learnMat(iMouse).iSeshv4);
    
    learnMat(iMouse).percMov = nan; learnMat(iMouse).movAllv1 = nan;
    learnMat(iMouse).RmovSesh = nan; learnMat(iMouse).RmovCmu = nan; 
    learnMat(iMouse).rewIntervals = nan; learnMat(iMouse).movmu = nan(size(seshTemp(1).movvC));
    learnMat(iMouse).movsd = nan(size(seshTemp(1).movvC)); learnMat(iMouse).nMov = nan;
    % May not have v1 sessions for everyone
    if ~isempty(learnMat(iMouse).iSeshv1)
        learnMat(iMouse).percMov = mean(cat(1,seshV1(learnMat(iMouse).iSeshv1).percMov));
        learnMat(iMouse).movAllv1 = cat(2,seshV1(learnMat(iMouse).iSeshv1).movAll);
        RmovSesh0 = corrcoef(learnMat(iMouse).movAllv1);
        iUpper = triu(true(size(RmovSesh0)), 1);
        learnMat(iMouse).RmovSesh = mean(RmovSesh0(iUpper));
        learnMat(iMouse).rewIntervals = median(cat(1,seshV1(learnMat(iMouse).iSeshv1).rewIntervals));
        movvC = -mean(cat(2,seshTemp(end-1:end).movvC),2);
        RmovC0 = corrcoef([learnMat(iMouse).movAllv1 movvC]);
        learnMat(iMouse).RmovCmu = mean(RmovC0(1:end-1,end));        
        learnMat(iMouse).movmu = mean(learnMat(iMouse).movAllv1,2);
        learnMat(iMouse).movsd = std(learnMat(iMouse).movAllv1,[],2);
        learnMat(iMouse).nMov = size(learnMat(iMouse).movAllv1,2);
    end
    learnMat(iMouse).RmovCmu = [learnMat(iMouse).RmovCmu; cat(1,seshTemp.RmovCmu)];
    learnMat(iMouse).percMov = [learnMat(iMouse).percMov; cat(1,seshTemp.percMov)];
    for jSesh = 1:numel(seshTemp)
        learnMat(iMouse).rewIntervals(jSesh+1,1) = median(seshTemp(jSesh).rewIntervals);
        learnMat(iMouse).movmu(:,jSesh+1) = -mean(seshTemp(jSesh).movv,2);
        learnMat(iMouse).movsd(:,jSesh+1) = std(seshTemp(jSesh).movv,[],2);
        learnMat(iMouse).nMov(jSesh+1,1) = size(seshTemp(jSesh).movv,2);
        RmovSesh0 = corrcoef(seshTemp(jSesh).movv);
        iUpper = triu(true(size(RmovSesh0)), 1);
        learnMat(iMouse).RmovSesh(jSesh+1) = mean(RmovSesh0(iUpper));
    end
    
    % mean movement correlation across sessions
    learnMat(iMouse).movmu_corr = corrcoef(learnMat(iMouse).movmu);


    tPlot = (0.5:1:numel(seshTemp(1).movvC))'/100-.5;
    f1 = figure;
    for iPlot = 1:9
        if ~isnan(learnMat(iMouse).movmu(1,iPlot))
            subplot(3,3,iPlot);
            h1 = boundedline(tPlot,learnMat(iMouse).movmu(:,iPlot),learnMat(iMouse).movsd(:,iPlot),'cmap',[0 0 0],'alpha');
%             h1 = boundedline(tPlot,learnMat(iMouse).movmu(:,iPlot),learnMat(iMouse).movsd(:,iPlot)/learnMat(iMouse).nMov(iPlot),'cmap',[0 0 0],'alpha');
            set(h1,'linewidth',1.5);
            ylim([-5 15]); xlim([-.5 1]); xticks([-.5:.25:1]); yticks([-5:2.5:15]); xtickangle(0);
            if iPlot <=6; xticklabels({}); end
            if ~any(iPlot == [1,4,7]); yticklabels({}); end
            if iPlot == 4; ylabel('Velocity (cm/s)'); end
            hold on; xline(0,'--b');
            if iPlot == 2; title(learnMat(iMouse).name); end
            set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
        end
    end
    if savebool
        if ~exist(figFold,'dir'); mkdir(figFold); end
        figname = sprintf('%s_v1_v4_movTraces',learnMat(iMouse).name);
        savefig(f1, fullfile(figFold,[figname,'.fig']));
        saveas(f1,fullfile(figFold,[figname,'.png']));
        print(f1,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');        
    end
        
end

%% Summary stuff
% Proportion of trials with movement:
percMov_v1 = cat(1,seshV1.percMov)*100;
percMov_v4 = cat(1,seshV4.percMov)*100;
[pPercMov, ~] = ranksum(percMov_v4, percMov_v1);

% Mean movements across age
movmuAll = cat(3,learnMat.movmu);
movmuAll_mu = mean(movmuAll,3,'omitnan');
movmuAll_sd = std(movmuAll,[],3,'omitnan'); %/sqrt(size(movmuAll,3));
movmu_naive_mu = movmuAll_mu(:,1);
movmu_naive_sd = movmuAll_sd(:,1);
movmu_inter_mu = mean(reshape(movmuAll(:,2:3,:), size(movmuAll,1), 2*size(movmuAll,3)), 2, 'omitnan');
movmu_inter_sd = std(reshape(movmuAll(:,2:3,:), size(movmuAll,1), 2*size(movmuAll,3)), [], 2, 'omitnan');
movmu_exp_mu = mean(reshape(movmuAll(:,end-1:end,:), size(movmuAll,1), 2*size(movmuAll,3)),2,'omitnan');
movmu_exp_sd = std(reshape(movmuAll(:,end-1:end,:), size(movmuAll,1), 2*size(movmuAll,3)), [], 2, 'omitnan');

% Movement correlations - within session:
RmovSeshAll = cat(1,learnMat.RmovSesh)';
RmovSesh_v1 = RmovSeshAll(1,:)';
RmovSesh_v4_1 = reshape(RmovSeshAll(2:3,:),2*(nSesh+1),[]);
RmovSesh_v4_2 = reshape(RmovSeshAll(4:5,:),2*(nSesh+1),[]);
RmovSesh_v4_3 = reshape(RmovSeshAll(6:7,:),2*(nSesh+1),[]);
RmovSesh_v4_n = reshape(RmovSeshAll(8:9,:),2*(nSesh+1),[]);

% Movement correlations - to trained movement:
RmovCAll = cat(2,learnMat.RmovCmu);
RmovC_v1 = RmovCAll(1,:)';
RmovC_v4_1 = reshape(RmovCAll(2:3,:),2*(nSesh+1),[]);
RmovC_v4_2 = reshape(RmovCAll(4:5,:),2*(nSesh+1),[]);
RmovC_v4_3 = reshape(RmovCAll(6:7,:),2*(nSesh+1),[]);
RmovC_v4_n = reshape(RmovCAll(8:9,:),2*(nSesh+1),[]);
v1_grp = ones(numel(RmovC_v1),1);
v4_1_grp = ones(numel(RmovC_v4_1),1)*2;
v4_2_grp = ones(numel(RmovC_v4_2),1)*3;
v4_3_grp = ones(numel(RmovC_v4_3),1)*4;
v4_n_grp = ones(numel(RmovC_v4_n),1)*5;

RmovSesh_vector = cat(1,RmovSesh_v1,RmovSesh_v4_1,RmovSesh_v4_2,RmovSesh_v4_3,RmovSesh_v4_n);
RmovC_vector = cat(1,RmovC_v1,RmovC_v4_1,RmovC_v4_2,RmovC_v4_3,RmovC_v4_n);
Rmov_grps = cat(1,v1_grp,v4_1_grp,v4_2_grp,v4_3_grp,v4_n_grp);
[movSesh_P, movSesh_Table, movSesh_Stats] = kruskalwallis(RmovSesh_vector, Rmov_grps);
RmovSesh_stats = multcompare(movSesh_Stats, 'Ctype','bonferroni');
[movC_P, movC_Table, movC_Stats] = kruskalwallis(RmovC_vector, Rmov_grps);
RmovC_stats = multcompare(movC_Stats, 'Ctype','bonferroni');

statsAll = v2struct(pPercMov, movSesh_P,movSesh_Table,movSesh_Stats, RmovC_stats, movC_P,movC_Table,movC_Stats, RmovC_stats);

RmovSesh_mu(1) = mean(RmovSesh_v1,'omitnan');
RmovSesh_sem(1) = std(RmovSesh_v1,'omitnan')/sqrt(sum(~isnan(RmovSesh_v1)));
RmovSesh_mu(2) = mean(RmovSesh_v4_1,'omitnan');
RmovSesh_sem(2) = std(RmovSesh_v4_1,'omitnan')/sqrt(2*(nSesh+1));
RmovSesh_mu(3) = mean(RmovSesh_v4_2,'omitnan');
RmovSesh_sem(3) = std(RmovSesh_v4_2,'omitnan')/sqrt(2*(nSesh+1));
RmovSesh_mu(4) = mean(RmovSesh_v4_3,'omitnan');
RmovSesh_sem(4) = std(RmovSesh_v4_3,'omitnan')/sqrt(2*(nSesh+1));
RmovSesh_mu(5) = mean(RmovSesh_v4_n,'omitnan');
RmovSesh_sem(5) = std(RmovSesh_v4_n,'omitnan')/sqrt(2*(nSesh+1));


RmovC_mu(1) = mean(RmovC_v1,'omitnan');
RmovC_sem(1) = std(RmovC_v1,'omitnan')/sqrt(sum(~isnan(RmovC_v1)));
RmovC_mu(2) = mean(RmovC_v4_1,'omitnan');
RmovC_sem(2) = std(RmovC_v4_1,'omitnan')/sqrt(2*(nSesh+1));
RmovC_mu(3) = mean(RmovC_v4_2,'omitnan');
RmovC_sem(3) = std(RmovC_v4_2,'omitnan')/sqrt(2*(nSesh+1));
RmovC_mu(4) = mean(RmovC_v4_3,'omitnan');
RmovC_sem(4) = std(RmovC_v4_3,'omitnan')/sqrt(2*(nSesh+1));
RmovC_mu(5) = mean(RmovC_v4_n,'omitnan');
RmovC_sem(5) = std(RmovC_v4_n,'omitnan')/sqrt(2*(nSesh+1));

rewIntAll0 = cat(2,learnMat.rewIntervals)';
rewIntAll = rewIntAll0(:);
Rgrp = zeros(81,1);
Rgrp(1:9) = 1; Rgrp(10:27) = 2; Rgrp(28:45) = 3; Rgrp(46:63) = 4; Rgrp(64:81) = 5;

% Plot mean movements across mice
f2 = figure;
for iPlot = 1:9
    subplot(3,3,iPlot);
    h1 = boundedline(tPlot,movmuAll_mu(:,iPlot),movmuAll_sd(:,iPlot),'cmap',[0 0 0],'alpha');    
    set(h1,'linewidth',1.5);
    ylim([-5 15]); xlim([-.5 1]); xticks([-.5:.25:1]); yticks([-5:2.5:15]); xtickangle(0);
    if iPlot <=6; xticklabels({}); end
    if ~any(iPlot == [1,4,7]); yticklabels({}); end
    if iPlot == 4; ylabel('Velocity (cm/s)'); end
    hold on; xline(0,'--b');
    if iPlot == 2; title(learnMat(iMouse).name); end
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
end

% Plot summary stuff
f3 = figure;
subplot(2,2,1); hold on;
hNaive = boundedline(tPlot, movmu_naive_mu, movmu_naive_sd, 'cmap',[0,0,0]); set(hNaive,'linewidth',1.5)
hInter = boundedline(tPlot, movmu_inter_mu, movmu_inter_sd, 'cmap',[0,0,0], 'alpha'); set(hInter,'linewidth',1.5)
hExp = boundedline(tPlot, movmu_exp_mu, movmu_exp_sd, 'cmap',[.5,0,0], 'alpha'); set(hExp,'linewidth',1.5)
ylim([-5 15]); xlim([-.5 1]); xticks([-.5:.25:1]); yticks([-5:2.5:15]); xtickangle(0);
xline(0,'--b');
xlabel('Time from movement (s)'); ylabel('Velocity (cm/s)');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');

subplot(2,2,2); hold on;
iDays = (1:5)';
h1 = boundedline(iDays,RmovSesh_mu,RmovSesh_sem,'cmap',[0 0 0],'alpha');
set(h1,'linewidth',2);
scatter(iDays,RmovSesh_mu,60,'ok','filled');
errorbar(iDays,RmovSesh_mu, RmovSesh_sem,'k','CapSize',0,'LineWidth',2);
xlim([0.5 5.5]); ylim([.25 .75]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
axis square; box off;
ylabel('Corr with trained movement');

subplot(2,2,4); hold on;
iDays = (1:5)';
h1 = boundedline(iDays,RmovC_mu,RmovC_sem,'cmap',[0 0 0],'alpha');
set(h1,'linewidth',2);
scatter(iDays,RmovC_mu,60,'ok','filled');
errorbar(iDays,RmovC_mu, RmovC_sem,'k','CapSize',0,'LineWidth',2);
xlim([0.5 5.5]); ylim([.25 .75]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
axis square; box off;
ylabel('Corr with trained movement');

subplot(2,2,3); hold on;
boxplot(rewIntAll(~isnan(rewIntAll)),Rgrp(~isnan(rewIntAll)),'BoxStyle','filled');
xlim([0.5 5.5]); ylim([0 85]); yticks([0:10:85]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
axis square; box off;
ylabel('Median reward interval');

if savebool
    figname = 'AllMice_v1_v4_movTraces';
    savefig(f2, fullfile(figFold,[figname,'.fig']));
    saveas(f2,fullfile(figFold,[figname,'.png']));
    print(f2,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');    

    figname = 'Summary_v1_v4_movTraces';
    savefig(f3, fullfile(figFold,[figname,'.fig']));
    saveas(f3,fullfile(figFold,[figname,'.png']));
    print(f3,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');    
    save(fullfile(figFold,'v1_v4_learnMat.mat'),'learnMat','statsAll');
end

end