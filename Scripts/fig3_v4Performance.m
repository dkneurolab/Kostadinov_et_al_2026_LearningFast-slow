function fig3_v4Performance(seshV4,dataFold,paperFold,savebool)

figFold = fullfile(paperFold,'Figures','Fig3','Sub1_behavSummary');
figFold2 = fullfile(dataFold,'paperFigs','Fig3','Sub1_behavSummary');



%% Loop through and get performance from behavstability .mat files
RmovUCO = zeros(numel(seshV4),3);
percUCO = RmovUCO;
grpUCO = RmovUCO; % Need 3 groups: Lob V ==> blue, Sim2 ==> magenta, no imaging ==> gray
grpRGB = [0 0 0.5; 0.5 0 0.5; .5 .5 .5];

for iSesh = 1:numel(seshV4)
    behavStab = load(fullfile(dataFold,'04_Behav-processed',seshV4(iSesh).name,seshV4(iSesh).date,'Behaviour_4s',sprintf('%s_%s_stability.mat',seshV4(iSesh).name,seshV4(iSesh).date)));
    behavStab = behavStab.behavstability;
    RmovUCO(iSesh,:) = [behavStab.params.rhoUmu behavStab.params.rhoCmu behavStab.params.rhoOmu];
    percUCO(iSesh,:) = [numel(behavStab.params.trialsU) numel(behavStab.params.trialsC) numel(behavStab.params.trialsO)]/behavStab.params.trialnum*100;
    if ~seshV4(iSesh).data2p
        grpUCO(iSesh,3) = 1;
    else
        if strcmpi(seshV4(iSesh).fov,'lobv')
            grpUCO(iSesh,1) = 1;
        else
            grpUCO(iSesh,2) = 1;
        end
    end
end

%% Now plot

% All data together
f1 = figure;
subplot(1,2,1); pbaspect([.5 .75 1]); hold on;
iDays = (1:3)';
h1 = boundedline(iDays,mean(percUCO),std(percUCO)/sqrt(size(percUCO,1)),'cmap',[0 0 0],'alpha');
set(h1,'linewidth',2);
scatter(iDays,mean(percUCO),60,'ok','filled');
errorbar(iDays,mean(percUCO), std(percUCO)/sqrt(size(percUCO,1)),'k','CapSize',0,'LineWidth',2);
xlim([0.5 3.5]); xticks([1 2 3]); xticklabels({'Under.','Corr.','Over.'});
ylim([0 80]); yticks([0:10:80]); ylabel('Percent correct');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
axis square; box off;

subplot(1,2,2); pbaspect([.5 .75 1]); hold on;
iDays = (1:3)';
h1 = boundedline(iDays,mean(RmovUCO),std(RmovUCO)/sqrt(size(RmovUCO,1)),'cmap',[0 0 0],'alpha');
set(h1,'linewidth',2);
scatter(iDays,mean(RmovUCO),60,'ok','filled');
errorbar(iDays,mean(RmovUCO), std(RmovUCO)/sqrt(size(RmovUCO,1)),'k','CapSize',0,'LineWidth',2);
xlim([0.5 3.5]); xticks([1 2 3]); xticklabels({'Under.','Corr.','Over.'});
ylim([.6 .8]); yticks([.6:.05:.8]); ylabel('Corr with correct');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
axis square; box off;

% Split by v4 recording sessions
f2 = figure;
subplot(1,2,1); pbaspect([.5 .75 1]); hold on;
iDays = (1:3)';
% plot(iDays-.1,percUCO(logical(grpUCO(:,1)),:),'Color',[grpRGB(1,:) 0.5]);
% plot(iDays+.1,percUCO(logical(grpUCO(:,2)),:),'Color',[grpRGB(2,:) 0.5]);
% plot(iDays,percUCO(logical(grpUCO(:,3)),:),'Color',[grpRGB(3,:) 0.5]);
scatter(iDays-.1,mean(percUCO(logical(grpUCO(:,1)),:)),60,'ok','filled','MarkerFaceColor',grpRGB(1,:));
errorbar(iDays-.1,mean(percUCO(logical(grpUCO(:,1)),:)), std(percUCO(logical(grpUCO(:,1)),:))/sqrt(size(percUCO(logical(grpUCO(:,1)),:),1)),'Color',grpRGB(1,:),'CapSize',0,'LineWidth',2,'LineStyle','none');
scatter(iDays+.1,mean(percUCO(logical(grpUCO(:,2)),:)),60,'ok','filled','MarkerFaceColor',grpRGB(2,:));
errorbar(iDays+.1,mean(percUCO(logical(grpUCO(:,2)),:)), std(percUCO(logical(grpUCO(:,2)),:))/sqrt(size(percUCO(logical(grpUCO(:,2)),:),1)),'Color',grpRGB(2,:),'CapSize',0,'LineWidth',2,'LineStyle','none');
xlim([0.5 3.5]); xticks([1 2 3]); xticklabels({'Under.','Corr.','Over.'});
ylim([0 80]); yticks([0:10:80]); ylabel('Percent correct');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
axis square; box off;

subplot(1,2,2); pbaspect([.5 .75 1]); hold on;
iDays = (1:3)';
scatter(iDays-.1,mean(RmovUCO(logical(grpUCO(:,1)),:)),60,'ok','filled','MarkerFaceColor',grpRGB(1,:));
errorbar(iDays-.1,mean(RmovUCO(logical(grpUCO(:,1)),:)), std(RmovUCO(logical(grpUCO(:,1)),:))/sqrt(size(RmovUCO(logical(grpUCO(:,1)),:),1)),'Color',grpRGB(1,:),'CapSize',0,'LineWidth',2,'LineStyle','none');
scatter(iDays+.1,mean(RmovUCO(logical(grpUCO(:,2)),:)),60,'ok','filled','MarkerFaceColor',grpRGB(2,:));
errorbar(iDays+.1,mean(RmovUCO(logical(grpUCO(:,2)),:)), std(RmovUCO(logical(grpUCO(:,2)),:))/sqrt(size(RmovUCO(logical(grpUCO(:,2)),:),1)),'Color',grpRGB(2,:),'CapSize',0,'LineWidth',2,'LineStyle','none');
xlim([0.5 3.5]); xticks([1 2 3]); xticklabels({'Under.','Corr.','Over.'});
ylim([.55 .85]); yticks([.55:.05:.85]); ylabel('Corr with correct');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
axis square; box off;

%% Stats:
RmovUCO_lobv = RmovUCO(logical(grpUCO(:,1)),:);
RmovUCO_sim2 = RmovUCO(logical(grpUCO(:,2)),:);

[a,b] = signrank(RmovUCO_lobv(:,1), RmovUCO_sim2(:,1));
RmovUCO_pvals(1,1) = a;
[a,b] = signrank(RmovUCO_lobv(:,2), RmovUCO_sim2(:,2));
RmovUCO_pvals(2,1) = a;
[a,b] = signrank(RmovUCO_lobv(:,3), RmovUCO_sim2(:,3));
RmovUCO_pvals(3,1) = a;

percUCO_lobv = percUCO(logical(grpUCO(:,1)),:);
percUCO_sim2 = percUCO(logical(grpUCO(:,2)),:);

[a,b] = signrank(percUCO_lobv(:,1), percUCO_sim2(:,1));
percUCO_pvals(1,1) = a;
[a,b] = signrank(percUCO_lobv(:,2), percUCO_sim2(:,2));
percUCO_pvals(2,1) = a;
[a,b] = signrank(percUCO_lobv(:,3), percUCO_sim2(:,3));
percUCO_pvals(3,1) = a;



%% Save?
if savebool
    if ~exist(figFold,'dir'); mkdir(figFold); end
    if ~exist(figFold2,'dir'); mkdir(figFold2); end
    save(fullfile(figFold,'v4Performance.mat'),'seshV4','RmovUCO','percUCO', 'RmovUCO_pvals', 'percUCO_pvals');
    save(fullfile(figFold2,'v4Performance.mat'),'seshV4','RmovUCO','percUCO', 'RmovUCO_pvals', 'percUCO_pvals');
    figname = 'v4Performance';
    savefig(f1, fullfile(figFold,[figname,'.fig']));
    saveas(f1,fullfile(figFold,[figname,'.png']));
    print(f1,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');
    figname = 'v4Performance_lobv_sim';
    savefig(f2, fullfile(figFold,[figname,'.fig']));
    saveas(f2,fullfile(figFold,[figname,'.png']));
    print(f2,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');
end

end