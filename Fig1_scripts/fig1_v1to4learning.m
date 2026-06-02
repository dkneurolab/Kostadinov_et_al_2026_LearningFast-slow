function fig1_v1to4learning(seshV1, seshV4,nSesh,dataFold,paperFold,savebool)
% Set directories for behavioural data and output figures
figFold = fullfile(paperFold,'Figures','Fig1','Sub1_v1Learning');

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
    learnMat(iMouse).RmovCmu = nan; learnMat(iMouse).rewIntervals = nan;
    learnMat(iMouse).movmu = nan(size(seshTemp(1).movvC));
    learnMat(iMouse).movsd = nan(size(seshTemp(1).movvC));
    learnMat(iMouse).nMov = nan;
    % May not have v1 sessions for everyone
    if ~isempty(learnMat(iMouse).iSeshv1)
        learnMat(iMouse).percMov = mean(cat(1,seshV1(learnMat(iMouse).iSeshv1).percMov));
        learnMat(iMouse).movAllv1 = cat(2,seshV1(learnMat(iMouse).iSeshv1).movAll);
%         learnMat(iMouse).rewIntervals{1} = cat(1,seshV1(learnMat(iMouse).iSeshv1).rewIntervals);
        learnMat(iMouse).rewIntervals = median(cat(1,seshV1(learnMat(iMouse).iSeshv1).rewIntervals));
        movvC = -mean(cat(2,seshTemp(end-1:end).movvC),2);
        Rmov0 = corrcoef([learnMat(iMouse).movAllv1 movvC]);
        learnMat(iMouse).RmovCmu = mean(Rmov0(1:end-1,end));        
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
    end
    
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
        if ~exist(figFold2,'dir'); mkdir(figFold); end
        figname = sprintf('%s_v1_v4_movTraces',learnMat(iMouse).name);
        savefig(f1, fullfile(figFold,[figname,'.fig']));
        saveas(f1,fullfile(figFold,[figname,'.png']));
        print(f1,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');
        savefig(f1, fullfile(figFold2,[figname,'.fig']));
        saveas(f1,fullfile(figFold2,[figname,'.png']));
        print(f1,fullfile(figFold2,[figname,'.eps']), '-depsc', '-vector');
    end
        
end

%% Summary stuff

% Proportion of trials with movement:
percMov_v1 = cat(1,seshV1.percMov)*100;
percMov_v4 = cat(1,seshV4.percMov)*100;
[pPercMov, ~] = ranksum(percMov_v4, percMov_v1);

% Movement correlations:
RmovAll = cat(2,learnMat.RmovCmu);
Rmov_v1 = RmovAll(1,:)';
Rmov_v4_1 = reshape(RmovAll(2:3,:),2*(nSesh+1),[]);
Rmov_v4_2 = reshape(RmovAll(4:5,:),2*(nSesh+1),[]);
Rmov_v4_3 = reshape(RmovAll(6:7,:),2*(nSesh+1),[]);
Rmov_v4_n = reshape(RmovAll(8:9,:),2*(nSesh+1),[]);
v1_grp = ones(numel(Rmov_v1),1);
v4_1_grp = ones(numel(Rmov_v4_1),1)*2;
v4_2_grp = ones(numel(Rmov_v4_2),1)*3;
v4_3_grp = ones(numel(Rmov_v4_3),1)*4;
v4_n_grp = ones(numel(Rmov_v4_n),1)*5;

Rmov_vector = cat(1,Rmov_v1,Rmov_v4_1,Rmov_v4_2,Rmov_v4_3,Rmov_v4_n);
Rmov_grps = cat(1,v1_grp,v4_1_grp,v4_2_grp,v4_3_grp,v4_n_grp);
[movP, movTable, movStats] = kruskalwallis(Rmov_vector, Rmov_grps);
Rmov_stats = multcompare(movStats, 'Ctype','bonferroni');


RmovMu(1) = mean(Rmov_v1,'omitnan');
RmovSem(1) = std(Rmov_v1,'omitnan')/sqrt(sum(~isnan(Rmov_v1)));
RmovMu(2) = mean(Rmov_v4_1,'omitnan');
RmovSem(2) = std(Rmov_v4_1,'omitnan')/sqrt(2*(nSesh+1));
RmovMu(3) = mean(Rmov_v4_2,'omitnan');
RmovSem(3) = std(Rmov_v4_2,'omitnan')/sqrt(2*(nSesh+1));
RmovMu(4) = mean(Rmov_v4_3,'omitnan');
RmovSem(4) = std(Rmov_v4_3,'omitnan')/sqrt(2*(nSesh+1));
RmovMu(5) = mean(Rmov_v4_n,'omitnan');
RmovSem(5) = std(Rmov_v4_n,'omitnan')/sqrt(2*(nSesh+1));

rewIntAll0 = cat(2,learnMat.rewIntervals)';
rewIntAll = rewIntAll0(:);
Rgrp = zeros(81,1);
Rgrp(1:9) = 1; Rgrp(10:27) = 2; Rgrp(28:45) = 3; Rgrp(46:63) = 4; Rgrp(64:81) = 5;



f2 = figure;
subplot(1,2,1); hold on;
iDays = (1:5)';
h1 = boundedline(iDays,RmovMu,RmovSem,'cmap',[0 0 0],'alpha');
set(h1,'linewidth',2);
scatter(iDays,RmovMu,60,'ok','filled');
errorbar(iDays,RmovMu, RmovSem,'k','CapSize',0,'LineWidth',2);
xlim([0.5 5.5]); ylim([.25 .75]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
axis square; box off;
ylabel('Corr with trained movement');

subplot(1,2,2); hold on;
boxplot(rewIntAll(~isnan(rewIntAll)),Rgrp(~isnan(rewIntAll)),'BoxStyle','filled');
xlim([0.5 5.5]); ylim([0 85]); yticks([0:10:85]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
axis square; box off;
ylabel('Median reward interval');

if savebool
    figname = 'Summary_v1_v4_movTraces';
    savefig(f2, fullfile(figFold,[figname,'.fig']));
    saveas(f2,fullfile(figFold,[figname,'.png']));
    print(f2,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');
    savefig(f2, fullfile(figFold2,[figname,'.fig']));
    saveas(f2,fullfile(figFold2,[figname,'.png']));
    print(f2,fullfile(figFold2,[figname,'.eps']), '-depsc', '-vector');
    save(fullfile(figFold,'v1_v4_learnat.mat'),'learnMat');
    save(fullfile(figFold2,'v1_v4_learnat.mat'),'learnMat');
end

end