function seshAll = fig1_v4learning(seshNaive,seshLearning,seshTrained, nLearn, dataFold,paperFold,savebool)

naiveNames = cellstr(cat(1,seshNaive.name));
percC_All = zeros(nLearn+2,numel(seshLearning)/nLearn);
percC_Good = percC_All;

% Histogram binning, etc.
nBinsCoarse = 28;
nBinsFine = 55;

histBins_coarse = linspace(-1,2,nBinsCoarse)';
histMids_coarse = histBins_coarse(1:end-1)+diff(histBins_coarse(1:2))/2;
fillBins_coarse = [histBins_coarse histBins_coarse]'; fillBins_coarse = fillBins_coarse(:); fillBins_coarse = [fillBins_coarse; fillBins_coarse(end)];
diffCoarse = histBins_coarse(2)-histBins_coarse(1);

histBins_fine = linspace(-1,2,nBinsFine)';
histMids_fine = histBins_fine(1:end-1)+diff(histBins_fine(1:2))/2;
fillBins_fine = [histBins_fine histBins_fine]'; fillBins_fine = fillBins_fine(:); fillBins_fine = [fillBins_fine; fillBins_fine(end)];
diffFine = histBins_fine(2)-histBins_fine(1);

%% Extract data
a = 0;
stairsV4_allCoarse = zeros(numel(histBins_coarse)-1,nLearn+2,numel(seshLearning)/nLearn);
stairsV4_allFine = zeros(numel(histBins_fine)-1,nLearn+2,numel(seshLearning)/nLearn);
stairsV1C_all = zeros(numel(histBins_coarse)-1,numel(seshLearning)/nLearn);
stairsV1F_all = zeros(numel(histBins_fine)-1,numel(seshLearning)/nLearn);
seshAll = [];
for iMouse = 1:nLearn:numel(seshLearning)
    a = a + 1;
    sesh = seshLearning(iMouse:iMouse+nLearn-1);
    sesh = cat(2,sesh,seshTrained((a-1)*2+1:a*2));
    
    % Calculate percent correct
    percC = zeros(numel(sesh),1); percCgood = percC;
    for iSesh = 1:numel(sesh)
        behavData = load(fullfile(dataFold,'04_Behav-processed',sesh(iSesh).name,sesh(iSesh).date,'Behaviour_4s',sprintf('%s_%s.mat',sesh(iSesh).name,sesh(iSesh).date)));
        inits = behavData.datastruct.inits; sets = behavData.datastruct.sets;
        behavData = behavData.datastruct.paqdata.behav;
        
        if strcmpi(sesh(iSesh).name,'dk194') & strcmp(sesh(iSesh).date,'2020-10-02')
            behavData.trials = behavData.trials(1:180); % Account for lack of motivation at end of this session
        end    
        % Remove non-responsive trials at end
        behavData.trials = behavData.trials(1:find(cat(1,behavData.trials.realtrial) > 0,1, 'last'));
        
        percC(iSesh) = numel(behavData.offLC)/sum(numel(behavData.offLC)+numel(behavData.offLINC));
        percCgood(iSesh) = sum(cat(1,behavData.trials.correct))/sum(cat(1,behavData.trials.realtrial));
        trialsReal = behavData.trials(logical(cat(1,behavData.trials.realtrial)));
        sesh(iSesh).endProbCoarse = 100*histcounts(-cat(1,trialsReal.endposreal),histBins_coarse,'Normalization','probability')';
        sesh(iSesh).endProbFine = 100*histcounts(-cat(1,trialsReal.endposreal),histBins_fine,'Normalization','probability')';        
        
        behavStab = load(fullfile(dataFold,'04_Behav-processed',sesh(iSesh).name,sesh(iSesh).date,'Behaviour_4s',sprintf('%s_%s_stability.mat',sesh(iSesh).name,sesh(iSesh).date)));
        behavStab = behavStab.behavstability;
        sesh(iSesh).iTrialsMov = logical(cat(1,behavData.trials.realtrial));
        sesh(iSesh).percMov = sum(sesh(iSesh).iTrialsMov)/numel(sesh(iSesh).iTrialsMov);
        sesh(iSesh).iTrialsC = logical(cat(1,trialsReal.correct));
        iEndsReal = cat(1,trialsReal.inds); iEndsReal = iEndsReal(:,2)/inits.paqfs;        
        sesh(iSesh).rewIntervals = diff(iEndsReal(sesh(iSesh).iTrialsC));
        sesh(iSesh).movv = behavStab.movvall(151:300,:);
        sesh(iSesh).movvC = mean(sesh(iSesh).movv(:,sesh(iSesh).iTrialsC),2);        
        sesh(iSesh).offx = behavStab.offxall(51:200,:);
        sesh(iSesh).offxC = mean(sesh(iSesh).offx(:,sesh(iSesh).iTrialsC),2);
    end
    
    sesh = fig1_v4learningMovRsq(sesh);
    seshAll = cat(2,seshAll,sesh);
    
    percC_All(:,a) = percC;
    percC_Good(:,a) = percCgood;
    stairsV4_allCoarse(:,:,a) = cat(2,sesh.endProbCoarse);
    stairsV4_allFine(:,:,a) = cat(2,sesh.endProbFine);
    
    % See if there are naive dataset for this mouse and add them if so
    iNaiveMatch = cellfun( @(x) strcmpi(x,sesh(1).name),naiveNames);
    seshNaiveLocal = seshNaive(iNaiveMatch);
    
    endProbCoarse_naive = zeros(numel(histBins_coarse)-1,numel(seshNaiveLocal));
    endProbFine_naive = zeros(numel(histBins_fine)-1,numel(seshNaiveLocal));
    for iNaive = 1:numel(seshNaiveLocal)
        endProbCoarse_naive(:,iNaive) = 100*histcounts(-seshNaiveLocal(iNaive).xStopAll,histBins_coarse,'Normalization','probability')';
        endProbFine_naive(:,iNaive) = 100*histcounts(-seshNaiveLocal(iNaive).xStopAll,histBins_fine,'Normalization','probability')';
    end
    endProbCoarse_naiveMu = mean(endProbCoarse_naive,2);
    endProbFine_naiveMu = mean(endProbFine_naive,2);
    stairsV1C_all(:,a) = endProbCoarse_naiveMu;
    stairsV1F_all(:,a) = endProbFine_naiveMu;
    
    
    % Plot endpoint histogram for each dataset in v1, v4_1, and v4_n
    stairsV1C = [endProbCoarse_naiveMu; endProbCoarse_naiveMu(end)];
    fillV1C = [stairsV1C stairsV1C]'; fillV1C = fillV1C(:); fillV1C = [fillV1C(1); fillV1C]; fillV1C([1,numel(fillV1C)]) = 0;     %#ok<AGROW>
    stairsV4_1C = mean([sesh(1).endProbCoarse sesh(2).endProbCoarse],2); stairsV4_1C = [stairsV4_1C; stairsV4_1C(end)]; %#ok<AGROW>
    fillV4_1C = [stairsV4_1C stairsV4_1C]'; fillV4_1C = fillV4_1C(:); fillV4_1C = [fillV4_1C(1); fillV4_1C]; fillV4_1C([1,numel(fillV4_1C)]) = 0;     %#ok<AGROW>
    stairsV4_nC = mean([sesh(end-1).endProbCoarse sesh(end).endProbCoarse],2); stairsV4_nC = [stairsV4_nC; stairsV4_nC(end)]; %#ok<AGROW>
    fillV4_nC = [stairsV4_nC stairsV4_nC]'; fillV4_nC = fillV4_nC(:); fillV4_nC = [fillV4_nC(1); fillV4_nC]; fillV4_nC([1,numel(fillV4_nC)]) = 0; %#ok<AGROW>    
    stairsV1F = [endProbFine_naiveMu; endProbFine_naiveMu(end)];
    fillV1F = [stairsV1F stairsV1F]'; fillV1F = fillV1F(:); fillV1F = [fillV1F(1); fillV1F]; fillV1F([1,numel(fillV1F)]) = 0;     %#ok<AGROW>
    stairsV4_1F = mean([sesh(1).endProbFine sesh(2).endProbFine],2); stairsV4_1F = [stairsV4_1F; stairsV4_1F(end)]; %#ok<AGROW>
    fillV4_1F = [stairsV4_1F stairsV4_1F]'; fillV4_1F = fillV4_1F(:); fillV4_1F = [fillV4_1F(1); fillV4_1F]; fillV4_1F([1,numel(fillV4_1F)]) = 0;     %#ok<AGROW>
    stairsV4_nF = mean([sesh(end-1).endProbFine sesh(end).endProbFine],2); stairsV4_nF = [stairsV4_nF; stairsV4_nF(end)]; %#ok<AGROW>
    fillV4_nF = [stairsV4_nF stairsV4_nF]'; fillV4_nF = fillV4_nF(:); fillV4_nF = [fillV4_nF(1); fillV4_nF]; fillV4_nF([1,numel(fillV4_nF)]) = 0; %#ok<AGROW>
    
    fHist = figure;
    subplot(2,2,1); pbaspect([1.5 1 1]); hold on;
    stairs(histBins_coarse,stairsV4_nC,'-k','LineWidth',2);
    fill(fillBins_coarse,fillV4_nC,'r','FaceAlpha',0.6,'EdgeColor','none');
    stairs(histBins_coarse,stairsV4_1C,'-k','LineWidth',2);
    fill(fillBins_coarse,fillV4_1C,'k','FaceAlpha',0.2,'EdgeColor','none');
    stairs(histBins_coarse,stairsV1C,'-k','LineWidth',2);
    fill(fillBins_coarse,fillV1C,'k','FaceAlpha',0.6,'EdgeColor','none');
    xline([-1/3 1/3],'--r');
    xlim([-1.05 1.05]); xticks([-1:1/3:1]); xtickangle(0);
    ylim([0 20]); yticks([0:5:20]);    
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out'); box off;
    ylabel('Probability (percent)');

    subplot(2,2,2); pbaspect([1.5 1 1]); hold on;
    stairs(histBins_fine,stairsV4_nF,'-k','LineWidth',2);
    fill(fillBins_fine,fillV4_nF,'r','FaceAlpha',0.6,'EdgeColor','none');
    stairs(histBins_fine,stairsV4_1F,'-k','LineWidth',2);
    fill(fillBins_fine,fillV4_1F,'k','FaceAlpha',0.2,'EdgeColor','none');
    stairs(histBins_fine,stairsV1F,'-k','LineWidth',2);
    fill(fillBins_fine,fillV1F,'k','FaceAlpha',0.6,'EdgeColor','none');
    xline([-1/3 1/3],'--r');
    xlim([-1.05 1.05]); xticks([-1:1/3:1]); xtickangle(0)
    ylim([0 12]); yticks([0:2:12]);
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out'); box off;
    
    subplot(2,2,3); pbaspect([1.5 1 1]); hold on;
    xPlot = [histMids_coarse(1)-diffCoarse; histMids_coarse; histMids_coarse(end)+diffCoarse];    
    fill(xPlot,[0; stairsV4_nC(1:end-1);0],'r','FaceAlpha',0.6,'EdgeColor','none');
    fill(xPlot,[0; stairsV4_1C(1:end-1);0],'k','FaceAlpha',0.2,'EdgeColor','none');
    fill(xPlot,[0; stairsV1C(1:end-1);0],'k','FaceAlpha',0.6,'EdgeColor','none');
    xline([-1/3 1/3],'--r');
    xlim([-1.05 1.05]); xticks([-1:1/3:1]); xtickangle(0);
    ylim([0 20]); yticks([0:5:20]);    
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out'); box off;
    ylabel('Probability (percent)');
    xlabel('Endpoint position');
    
    subplot(2,2,4); pbaspect([1.5 1 1]); hold on;
    xPlot = [histMids_fine(1)-diffFine; histMids_fine; histMids_fine(end)+diffFine];    
    fill(xPlot,[0; stairsV4_nF(1:end-1);0],'r','FaceAlpha',0.6,'EdgeColor','none');
    fill(xPlot,[0; stairsV4_1F(1:end-1);0],'k','FaceAlpha',0.2,'EdgeColor','none');
    fill(xPlot,[0; stairsV1F(1:end-1);0],'k','FaceAlpha',0.6,'EdgeColor','none');
    xline([-1/3 1/3],'--r');
    xlim([-1.05 1.05]); xticks([-1:1/3:1]); xtickangle(0);
    ylim([0 12]); yticks([0:2:12]);
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out'); box off;
    xlabel('Endpoint position');
    
    suptitle([sesh(1).name,' endpoints - naive (dark), learning (light), and expert (red)']);
    
    if savebool
        figFold = fullfile(paperFold,'Figures','Fig1','Sub1_v4Learning');
        if ~exist(figFold,'dir'); mkdir(figFold); end
        figname = [sesh(1).name,'_endProb'];
        savefig(fHist, fullfile(figFold,[figname,'.fig']));
        saveas(fHist,fullfile(figFold,[figname,'.png']));
        print(fHist,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');
    end
end


%% Plot it
% Make coarse average traces
stairsV4_1C_all = squeeze(mean(stairsV4_allCoarse(:,1:2,:),2));
stairsV4_nC_all = squeeze(mean(stairsV4_allCoarse(:,end-1:end,:),2));
stairsV1C_all = stairsV1C_all(:,~isnan(stairsV1C_all(1,:)));

stairsV1C_mu = mean(stairsV1C_all,2); stairsV1C_muPlot = [stairsV1C_mu; stairsV1C_mu(end)];
stairsV1C_sem = std(stairsV1C_all,[],2)/sqrt(size(stairsV1C_all,2));
stairsV4_1C_mu = mean(stairsV4_1C_all,2); stairsV4_1C_muPlot = [stairsV4_1C_mu; stairsV4_1C_mu(end)];
stairsV4_1C_sem = std(stairsV4_1C_all,[],2)/sqrt(size(stairsV4_1C_all,2));
stairsV4_nC_mu = mean(stairsV4_nC_all,2); stairsV4_nC_muPlot = [stairsV4_nC_mu; stairsV4_nC_mu(end)];
stairsV4_nC_sem = std(stairsV4_nC_all,[],2)/sqrt(size(stairsV4_nC_all,2));
fillV1C_mu = [stairsV1C_muPlot stairsV1C_muPlot]'; fillV1C_mu = fillV1C_mu(:); fillV1C_mu = [fillV1C_mu(1); fillV1C_mu]; fillV1C_mu([1,numel(fillV1C_mu)]) = 0; %#ok<*NASGU>
fillV4_1C_mu = [stairsV4_1C_muPlot stairsV4_1C_muPlot]'; fillV4_1C_mu = fillV4_1C_mu(:); fillV4_1C_mu = [fillV4_1C_mu(1); fillV4_1C_mu]; fillV4_1C_mu([1,numel(fillV4_1C_mu)]) = 0; %#ok<*NASGU>
fillV4_nC_mu = [stairsV4_nC_muPlot stairsV4_nC_muPlot]'; fillV4_nC_mu = fillV4_nC_mu(:); fillV4_nC_mu = [fillV4_nC_mu(1); fillV4_nC_mu]; fillV4_nC_mu([1,numel(fillV4_nC_mu)]) = 0; %#ok<*NASGU>

% Make fine average traces
stairsV4_1F_all = squeeze(mean(stairsV4_allFine(:,1:2,:),2));
stairsV4_nF_all = squeeze(mean(stairsV4_allFine(:,end-1:end,:),2));
stairsV1F_all = stairsV1F_all(:,~isnan(stairsV1F_all(1,:)));

stairsV1F_mu = mean(stairsV1F_all,2); stairsV1F_muPlot = [stairsV1F_mu; stairsV1F_mu(end)];
stairsV1F_sem = std(stairsV1F_all,[],2)/sqrt(size(stairsV1F_all,2));
stairsV4_1F_mu = mean(stairsV4_1F_all,2); stairsV4_1F_muPlot = [stairsV4_1F_mu; stairsV4_1F_mu(end)];
stairsV4_1F_sem = std(stairsV4_1F_all,[],2)/sqrt(size(stairsV4_1F_all,2));
stairsV4_nF_mu = mean(stairsV4_nF_all,2); stairsV4_nF_muPlot = [stairsV4_nF_mu; stairsV4_nF_mu(end)];
stairsV4_nF_sem = std(stairsV4_nF_all,[],2)/sqrt(size(stairsV4_nF_all,2));
fillV1F_mu = [stairsV1F_muPlot stairsV1F_muPlot]'; fillV1F_mu = fillV1F_mu(:); fillV1F_mu = [fillV1F_mu(1); fillV1F_mu]; fillV1F_mu([1,numel(fillV1F_mu)]) = 0; %#ok<*NASGU>
fillV4_1F_mu = [stairsV4_1F_muPlot stairsV4_1F_muPlot]'; fillV4_1F_mu = fillV4_1F_mu(:); fillV4_1F_mu = [fillV4_1F_mu(1); fillV4_1F_mu]; fillV4_1F_mu([1,numel(fillV4_1F_mu)]) = 0; %#ok<*NASGU>
fillV4_nF_mu = [stairsV4_nF_muPlot stairsV4_nF_muPlot]'; fillV4_nF_mu = fillV4_nF_mu(:); fillV4_nF_mu = [fillV4_nF_mu(1); fillV4_nF_mu]; fillV4_nF_mu([1,numel(fillV4_nF_mu)]) = 0; %#ok<*NASGU>



fHistAll = figure;
subplot(2,2,1); pbaspect([1.5 1 1]); hold on;
stairs(histBins_coarse,stairsV4_nC_muPlot,'-k','LineWidth',2);
fill(fillBins_coarse,fillV4_nC_mu,'r','FaceAlpha',0.6,'EdgeColor','none');
stairs(histBins_coarse,stairsV4_1C_muPlot,'-k','LineWidth',2);
fill(fillBins_coarse,fillV4_1C_mu,'k','FaceAlpha',0.2,'EdgeColor','none');
stairs(histBins_coarse,stairsV1C_muPlot,'-k','LineWidth',2);
fill(fillBins_coarse,fillV1C_mu,'k','FaceAlpha',0.6,'EdgeColor','none');
errorbar(histMids_coarse,stairsV1C_mu,stairsV1C_sem,'k','CapSize',0,'LineWidth',2,'LineStyle','none');
errorbar(histMids_coarse,stairsV4_1C_mu,stairsV4_1C_sem,'k','CapSize',0,'LineWidth',2,'LineStyle','none');
errorbar(histMids_coarse,stairsV4_nC_mu,stairsV4_nC_sem,'k','CapSize',0,'LineWidth',2,'LineStyle','none');
xline([-1/3 1/3],'--r');
xlim([-1.05 1.05]); xticks([-1:1/3:1]); xtickangle(0)
ylim([0 20]); yticks([0:5:20]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out'); box off;
title('Coarse bars');
ylabel('Probability (percent)');

subplot(2,2,2); pbaspect([1.5 1 1]); hold on;
stairs(histBins_fine,stairsV4_nF_muPlot,'-k','LineWidth',2);
fill(fillBins_fine,fillV4_nF_mu,'r','FaceAlpha',0.6,'EdgeColor','none');
stairs(histBins_fine,stairsV4_1F_muPlot,'-k','LineWidth',2);
fill(fillBins_fine,fillV4_1F_mu,'k','FaceAlpha',0.2,'EdgeColor','none');
stairs(histBins_fine,stairsV1F_muPlot,'-k','LineWidth',2);
fill(fillBins_fine,fillV1F_mu,'k','FaceAlpha',0.6,'EdgeColor','none');
errorbar(histMids_fine,stairsV1F_mu,stairsV1F_sem,'k','CapSize',0,'LineWidth',2,'LineStyle','none');
errorbar(histMids_fine,stairsV4_1F_mu,stairsV4_1F_sem,'k','CapSize',0,'LineWidth',2,'LineStyle','none');
errorbar(histMids_fine,stairsV4_nF_mu,stairsV4_nF_sem,'k','CapSize',0,'LineWidth',2,'LineStyle','none');
xline([-1/3 1/3],'--r');
xlim([-1.05 1.05]); xticks([-1:1/3:1]); xtickangle(0)
ylim([0 10]); yticks([0:2:10]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out'); box off;
title('Fine bars');

subplot(2,2,3); pbaspect([1.5 1 1]); hold on;
h0 = boundedline(histMids_coarse,stairsV4_nC_mu,stairsV4_nC_sem,'cmap',[.5 0 0],'alpha'); set(h0,'linewidth',1.5);
h1 = boundedline(histMids_coarse,stairsV4_1C_mu,stairsV4_1C_sem,'cmap',[.5 .5 .5],'alpha'); set(h1,'linewidth',1.5);
h2 = boundedline(histMids_coarse,stairsV1C_mu,stairsV1C_sem,'cmap',[0 0 0],'alpha'); set(h2,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
xline([-1/3 1/3],'--r');
xlim([-1.05 1.05]); xticks([-1:1/3:1]); xtickangle(0)
ylim([0 20]); yticks([0:5:20]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out'); box off;
title('Coarse lines');
ylabel('Probability (percent)');
xlabel('Endpoint position');

subplot(2,2,4); pbaspect([1.5 1 1]); hold on;
h0 = boundedline(histMids_fine,stairsV4_nF_mu,stairsV4_nF_sem,'cmap',[.5 0 0],'alpha'); set(h0,'linewidth',1.5);
h1 = boundedline(histMids_fine,stairsV4_1F_mu,stairsV4_1F_sem,'cmap',[.5 .5 .5],'alpha'); set(h1,'linewidth',1.5);
h2 = boundedline(histMids_fine,stairsV1F_mu,stairsV1F_sem,'cmap',[0 0 0],'alpha'); set(h2,'linewidth',1.5);
xline([-1/3 1/3],'--r');
xlim([-1.05 1.05]); xticks([-1:1/3:1]); xtickangle(0)
ylim([0 10]); yticks([0:2:10]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out'); box off;
title('Fine lines');
xlabel('Endpoint position');
suptitle(['Summary v4 endpoints - naive (dark), learning (light), and expert (red)']);

f1 = figure;
subplot(1,2,1); hold on;

iDays = (1:nLearn)';
h1 = boundedline(iDays, mean(percC_All(1:nLearn,:),2), std(percC_All(1:nLearn,:),[],2)/sqrt(numel(seshLearning)/nLearn), 'alpha','cmap',[0 0 0]);
set(h1,'linewidth',2);
scatter([1:nLearn], mean(percC_All(1:nLearn,:),2),60,'ok','filled');
errorbar(iDays, mean(percC_All(1:nLearn,:),2), std(percC_All(1:6,:),[],2)/sqrt(numel(seshLearning)/nLearn),'k','CapSize',0,'LineWidth',2);
v4n = percC_All(nLearn+1:end,:); v4n = v4n(:);
scatter(nLearn+1.5, mean(v4n),60,'ok','filled');
errorbar(nLearn+1.5,mean(v4n),std(v4n)/sqrt(numel(v4n)),'k','CapSize',0,'LineWidth',2);
axis([.5 nLearn+2.5 .25 .75]); axis square; box off
xticks([0:nLearn+3]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
title('v4 learning curves');
ylabel({'Task performance','(% correct)'});
xlabel({'Days training on','expert task'});


perC_all2 = reshape(percC_All',[],nLearn/2+1)';
subplot(1,2,2); hold on;
iDays2 = (1:nLearn/2)';
h1 = boundedline(iDays2, mean(perC_all2(1:nLearn/2,:),2), std(percC_All(1:nLearn/2,:),[],2)/sqrt(size(perC_all2,2)), 'alpha','cmap',[0 0 0]);
set(h1,'linewidth',2);
scatter([1:nLearn/2], mean(perC_all2(1:nLearn/2,:),2),60,'ok','filled');
errorbar(iDays2, mean(perC_all2(1:nLearn/2,:),2), std(percC_All(1:nLearn/2,:),[],2)/sqrt(size(perC_all2,2)),'k','CapSize',0,'LineWidth',2);
v4n2 = perC_all2(nLearn/2+1:end,:); v4n2 = v4n2(:);
scatter(nLearn/2+1, mean(v4n2),60,'ok','filled');
errorbar(nLearn/2+1,mean(v4n2),std(v4n2)/sqrt(numel(v4n2)),'k','CapSize',0,'LineWidth',2);
axis([.5 nLearn/2+1.5 .25 .75]); axis square; box off
xticks([0:nLearn/2+2]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
title('v4 learning curves - all trials');
xlabel({'Days training on','expert task'});

f2 = figure;
subplot(1,2,1); hold on;

iDays = (1:nLearn)';
h1 = boundedline(iDays, mean(percC_Good(1:nLearn,:),2), std(percC_Good(1:nLearn,:),[],2)/sqrt(numel(seshLearning)/nLearn), 'alpha','cmap',[0 0 0]);
set(h1,'linewidth',2);
scatter([1:nLearn], mean(percC_Good(1:nLearn,:),2),60,'ok','filled');
errorbar(iDays, mean(percC_Good(1:nLearn,:),2), std(percC_Good(1:6,:),[],2)/sqrt(numel(seshLearning)/nLearn),'k','CapSize',0,'LineWidth',2);
v4n = percC_Good(nLearn+1:end,:); v4n = v4n(:);
scatter(nLearn+1.5, mean(v4n),60,'ok','filled');
errorbar(nLearn+1.5,mean(v4n),std(v4n)/sqrt(numel(v4n)),'k','CapSize',0,'LineWidth',2);
axis([.5 nLearn+2.5 .25 .75]); axis square; box off
xticks([0:nLearn+3]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
title('v4 learning curves');
ylabel({'Task performance','(% correct)'});
xlabel({'Days training on','expert task'});


perC_allGood2 = reshape(percC_Good',[],nLearn/2+1)';
subplot(1,2,2); hold on;
iDays2 = (1:nLearn/2)';
h1 = boundedline(iDays2, mean(perC_allGood2(1:nLearn/2,:),2), std(perC_allGood2(1:nLearn/2,:),[],2)/sqrt(size(perC_allGood2,2)), 'alpha','cmap',[0 0 0]);
set(h1,'linewidth',2);
scatter([1:nLearn/2], mean(perC_allGood2(1:nLearn/2,:),2),60,'ok','filled');
errorbar(iDays2, mean(perC_allGood2(1:nLearn/2,:),2), std(perC_allGood2(1:nLearn/2,:),[],2)/sqrt(size(perC_allGood2,2)),'k','CapSize',0,'LineWidth',2);
v4n2 = perC_allGood2(nLearn/2+1:end,:); v4n2 = v4n2(:);
scatter(nLearn/2+1, mean(v4n2),60,'ok','filled');
errorbar(nLearn/2+1,mean(v4n2),std(v4n2)/sqrt(numel(v4n2)),'k','CapSize',0,'LineWidth',2);
axis([.5 nLearn/2+1.5 .4 .7]); axis square; box off
xticks([0:nLearn/2+2]); xticklabels({'','1-2','3-4','5-6','Expert'});
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
title('v4 learning curves - good trials');
xlabel({'Days training on','expert task'});

[percC_P, percC_Table, percC_Stats] = kruskalwallis(perC_allGood2);
percC_Stats_post = multcompare(percC_Stats, 'Ctype','bonferroni');

statsAll = v2struct(percC_P,percC_Table,percC_Stats, percC_Stats_post);

if savebool    
    figname = 'v4LearningCurves_AllTrials';
    savefig(f1, fullfile(figFold,[figname,'.fig']));
    saveas(f1,fullfile(figFold,[figname,'.png']));
    print(f1,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');
    figname = 'v4LearningCurves_GoodTrials';
    savefig(f2, fullfile(figFold,[figname,'.fig']));
    saveas(f2,fullfile(figFold,[figname,'.png']));
    print(f2,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');
    figname = 'v4LearningCurves_endProb';
    savefig(fHistAll, fullfile(figFold,[figname,'.fig']));
    saveas(fHistAll,fullfile(figFold,[figname,'.png']));
    print(fHistAll,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');
    
    save(fullfile(figFold,'v4LearningCurves.mat'),'percC_All','percC_Good','histBins_coarse','histBins_fine','stairsV4_allCoarse','stairsV4_allFine','stairsV1C_all','stairsV1F_all','seshAll','seshNaive', 'statsAll');    

end

end
