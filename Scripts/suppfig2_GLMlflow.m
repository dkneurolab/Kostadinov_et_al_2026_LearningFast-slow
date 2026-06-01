% function suppfig2_GLMflow(sessionsv4_n, iExample, savebool)

%% Choose example dataset:
iExample = 6; % DK169, lobule V
glmRoot = fullfile(dataFold,'01_Behav+imaging',sprintf('Version_%i',sessionsv4_n(iExample).version),sessionsv4_n(iExample).name,...
        sprintf('%s_%s',sessionsv4_n(iExample).date,sessionsv4_n(iExample).fov),'GLM','glmRfiles');
load(fullfile(glmRoot,'GLMfull_small.mat'))

%% Reindexing design matrix:

cumBases = [0; cumsum(DMdata.nBases)];

baseRange = [cumBases(1:end-1)+1 cumBases(2:end)];

baseReind = [1,10,9,11,12,2,3,4,5,6,7,8];

DM2 = zeros(size(DM));

a = 0;
for i = 1:size(baseRange,1)
    iRange = baseRange(baseReind(i),1):baseRange(baseReind(i),2);
    DM2(:,a+1:a+numel(iRange)) = DM(:,iRange);
    a = a + numel(iRange);
end


trialRange = [5:12];
cumBases2 = [0; cumsum(DMdata.nBases(baseReind))];
cumTrialInds = DMdata.endTrialIndices(trialRange)-DMdata.endTrialIndices(trialRange(1));
iTrialRange = [DMdata.endTrialIndices(trialRange(1))+1:DMdata.endTrialIndices(trialRange(end)+1)];

%% Plot matrices:
f1 = figure;
ax1 = subplot(1,3,1);
imagesc(DM2(iTrialRange,:)'); colormap(ax1,flipud(gray)); clim([-.05 1]);
yline(cumBases2 + .5,'--k', LineWidth = 1)
xline(cumTrialInds+.5,'--k', LineWidth = 1)
yticks(cumBases2+.5)
xticks(cumTrialInds+.5)
set(ax1,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
box off; axis square;

ax2 = subplot(1,3,2);
imagesc(DMpca.score(iTrialRange,:)'); colormap(ax2,redblue); clim([-1 1])
xline(cumTrialInds+.5,'--k', LineWidth = 1)
xticks(cumTrialInds+.5)
set(ax2,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
box off; axis square;

ax3 = subplot(1,3,3);
imagesc(DMpca.score(iTrialRange,1:DMpca.nPCvar)'); colormap(ax3,redblue); clim([-1 1])
xline(cumTrialInds+.5,'--k', LineWidth = 1)
xticks(cumTrialInds+.5)
set(ax3,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
box off; axis square;

if savebool
    figname = 'DK169_lobv_predictorMats';
    figFold = 'C:\Users\Dimitar\Dropbox\DK papers\2022_Kostadinov_FastSlowLearning\Figures\SuppFig2';
    savefig(f1, fullfile(figFold,[figname,'.fig']));
    saveas(f1,fullfile(figFold,[figname,'.png']));
    print(f1,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');
end

%% Deviance explained scatter vs full mode

% GLM parameters:
glmParams.dataType      = 'spk'; % 'fluo' or 'spk' for zF or spike events
glmParams.distr         = 'gaussian'; % 'gaussian' (for fluo data), 'poisson' (for spk counts), or 'binomial' (for binarized spk trains)
glmParams.link          = []; % Use default link function for given distribution
glmParams.alpha         = .5; % Others use 0.5 for gaussian GLMs and 0.8-0.9 for poisson or binomial GLMs
glmParams.nLambda       = 100; % Good number, don't change from 100



devExpSummary = struct;

for iSesh = 1:numel(sessionsv4_n)
    % iSesh = 1;
    glmRoot = fullfile(dataFold,'01_Behav+imaging',sprintf('Version_%i',sessionsv4_n(iSesh).version),sessionsv4_n(iSesh).name,...
        sprintf('%s_%s',sessionsv4_n(iSesh).date,sessionsv4_n(iSesh).fov),'GLM','glmRfiles');
    glmFold = fullfile(glmRoot,sprintf('GLMfull_%s_%s_alpha%i',glmParams.dataType,glmParams.distr,glmParams.alpha*100));
    GLMoutput = load(fullfile(glmFold,'GLMsig.mat')); GLMoutput = GLMoutput.GLMoutput;
    devExpSummary(iSesh).devExp_full = cat(1,GLMoutput.devExp_full);
    devExpSummary(iSesh).devExp_best = cat(1,GLMoutput.devExp);
    devExpSummary(iSesh).devExp_shuf = cat(2,GLMoutput.devExp_shuf)';
    devExpSummary(iSesh).nPCs = cat(1,GLMoutput.GLMrankPCA);
    devExpSummary(iSesh).sigBool = cat(1,GLMoutput.sigBool);
end

%% Now plot
sigBool = logical(cat(1,devExpSummary.sigBool));
devExp_best = cat(1,devExpSummary.devExp_best);
devExp_full = cat(1,devExpSummary.devExp_full);
f2 = figure;
ax1 = subplot(1,2,1); title('All cells')
scatter(devExp_best,devExp_full,'MarkerFaceColor','k','MarkerEdgeColor','none'); axis square;  hold on;
refline(1,0); xline(0,'k'); yline(0,'k');
set(ax1,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');

ax2 = subplot(1,2,2); title('split by significance')
scatter(devExp_best(sigBool),devExp_full(sigBool),'MarkerFaceColor','r','MarkerEdgeColor','none'); axis square;  hold on;
scatter(devExp_best(~sigBool),devExp_full(~sigBool),'MarkerFaceColor','b','MarkerEdgeColor','none'); axis square;  hold on;
refline(1,0); xline(0,'k'); yline(0,'k');
set(ax2,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');

if savebool
    figname = 'devExpScatters';
    figFold = 'C:\Users\Dimitar\Dropbox\DK papers\2022_Kostadinov_FastSlowLearning\Figures\SuppFig2';
    savefig(f2, fullfile(figFold,[figname,'.fig']));
    saveas(f2,fullfile(figFold,[figname,'.png']));
    print(f2,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');
end


