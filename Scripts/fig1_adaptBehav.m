function fig1_adaptBehav(sessionsv5,trialGrp,dataFold,paperFold,savebool)
%% Make histogram data bins, etc.

nBinsCoarse = 28;
nBinsFine = 55;

histBins_coarse = linspace(-1,2,nBinsCoarse)';
histMids_coarse = histBins_coarse(1:end-1)+diff(histBins_coarse(1:2))/2;
% fillBins_coarse = [histBins_coarse histBins_coarse]'; fillBins_coarse = fillBins_coarse(:); fillBins_coarse = [fillBins_coarse; fillBins_coarse(end)];

histBins_fine = linspace(-1,2,nBinsFine)';
histMids_fine = histBins_fine(1:end-1)+diff(histBins_fine(1:2))/2;
% fillBins_fine = [histBins_fine histBins_fine]'; fillBins_fine = fillBins_fine(:); fillBins_fine = [fillBins_fine; fillBins_fine(end)];

endProbAll_coarse = zeros(numel(histMids_coarse),5,numel(sessionsv5));
endProbAll_fine = zeros(numel(histMids_fine),5,numel(sessionsv5));

%% Loop through behavioural sessions and pull out trials at different points
for iSesh = 1:numel(sessionsv5)
    
    behavStruct = load(fullfile(dataFold,'04_Behav-processed',sessionsv5(iSesh).name,sessionsv5(iSesh).date,'Behaviour_4s',sprintf('%s_%s.mat',sessionsv5(iSesh).name,sessionsv5(iSesh).date)));
    behavData = behavStruct.datastruct.paqdata.behav;
    inits = behavStruct.datastruct.inits;
    sets = Rot2p_changetrials(behavStruct.datastruct.sets);
    if strcmpi(sessionsv5(iSesh).name,'dk105') & strcmp(sessionsv5(iSesh).date,'2018-03-07')
        sets.changebacktrial = sets.changetrial+10*(sets.changebackmult-1);
    end
    
    trialsGood = behavData.trials(cat(1,behavData.trials.realtrial) > 0);
    behavStability = Rot2p_behavstability(trialsGood,inits,sets,0); close    
    
    iTrialsGood = cat(1,trialsGood.trialnum);
    iTrialsBase = [sets.changetrial-trialGrp+1:sets.changetrial]';
    iTrialsEarly = [sets.changetrial+1:sets.changetrial+trialGrp]';
    if ~sets.changeback
        iTrialsWash = [];
        sets.changebacktrial = iTrialsGood(end);
    else
        iTrialsWash = [sets.changebacktrial+1:sets.changebacktrial+trialGrp]';
    end
    
    trialsAdapt = trialsGood(iTrialsGood > sets.changetrial & iTrialsGood <= sets.changebacktrial);
    iRange = [trialsAdapt(1).trialnum trialsAdapt(1).trialnum+trialGrp-1];
    iTrialsAdapt = cat(1,trialsAdapt.trialnum);
    endPosPerf = struct;
    a = 0;
    while iRange(2) <= sets.changebacktrial
        trialsAdapt_range = trialsAdapt(find(iTrialsAdapt >= iRange(1),1,'first'):find(iTrialsAdapt <= iRange(2),1,'last'));
        endPos = cat(1,trialsAdapt_range.endpos);
        endPosReal = cat(1,trialsAdapt_range.endposreal);
        a = a+1;
        endPosPerf(a).iRange = iRange;
        endPosPerf(a).nUCOgain = [sum(endPos > 1/3) sum(endPos <= 1/3 & endPos >= -1/3) sum(endPos < -1/3)];
        endPosPerf(a).nUCOreal = [sum(endPosReal > 1/3) sum(endPosReal <= 1/3 & endPosReal >= -1/3) sum(endPosReal < -1/3)];
        endPosPerf(a).nOutcomes = sum(cat(1,trialsAdapt_range.outcome));
        endPosPerf(a).ratioGainReal = endPosPerf(a).nOutcomes(2)/endPosPerf(a).nUCOreal(2);
        endPosPerf(a).endPos = endPos;
        endPosPerf(a).endPosReal = endPosReal;
        endPosPerf(a).nTrials = numel(endPosReal);
        iRange = iRange + 1; 
    end
    % Find indices of worst adaptation range - if multiple, take the once
    % that's closes to the middle of the range for the minimum and the
    % middle one (arbitrary) for the maximum one
    
    endPosPerf = endPosPerf(cat(1,endPosPerf.nTrials) > trialGrp*.6);
    minRatio = min(cat(1,endPosPerf.ratioGainReal));
    iMinRatio = find(cat(1,endPosPerf.ratioGainReal) == minRatio);
    [~,iMinCtr] = min(abs(numel(endPosPerf)/2-iMinRatio));
    iMinRatio = iMinRatio(iMinCtr);
    maxRatio = max(cat(1,endPosPerf(iMinRatio+1:end).ratioGainReal));
    iMaxRatio = find(cat(1,endPosPerf(iMinRatio+1:end).ratioGainReal) == maxRatio)+iMinRatio;
    iMaxRatio = iMaxRatio(ceil(end/2));
    iTrialsMid = [endPosPerf(iMinRatio).iRange(1):endPosPerf(iMinRatio).iRange(2)]';
    iTrialsLate = [endPosPerf(iMaxRatio).iRange(1):endPosPerf(iMaxRatio).iRange(2)]';
    
    clear a endPos endPosReal iMinCtr iRange iTrialsAdapt maxRatio minRatio trialsAdapt_range
    
    trialInds = v2struct(iTrialsGood,iTrialsBase,iTrialsEarly,iTrialsMid,iTrialsLate,iTrialsWash);
    
    [singleSesh,endProbCoarse,endProbFine] = fig1_adaptBehav_sortTrials(sessionsv5(iSesh),trialInds,trialsGood,behavStability,histBins_coarse,histBins_fine);
    seshOut(iSesh) = singleSesh; %#ok<AGROW>
    endProbAll_coarse(:,:,iSesh) = endProbCoarse;
    endProbAll_fine(:,:,iSesh) = endProbFine;

    
    % Figure for each FOV
    fHist = figure;
    subplot(2,2,1); pbaspect([1.5 1 1]); hold on;
    xPlot = [histMids_coarse(1); histMids_coarse; histMids_coarse(end)];    
    fill(xPlot,[0; endProbCoarse(:,1);0],'k','FaceAlpha',0.8,'EdgeColor','none');
    fill(xPlot,[0; endProbCoarse(:,2);0],'r','FaceAlpha',0.5,'EdgeColor','none');
    fill(xPlot,[0; endProbCoarse(:,3);0],'m','FaceAlpha',0.5,'EdgeColor','none');
    fill(xPlot,[0; endProbCoarse(:,4);0],'b','FaceAlpha',0.5,'EdgeColor','none');
    fill(xPlot,[0; endProbCoarse(:,5);0],'k','FaceAlpha',0.3,'EdgeColor','none');    
    xline([-1/3 1/3],'--r');
    xlim([-1.1 1.1]); xticks([-1:1/3:1]); xticklabels({});
    ylim([0 30]); yticks([0:5:30]);
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out'); box off;
    title('Coarse binning');
    ylabel('Probability (percent)');
%     xlabel('Endpoint position');
    
    subplot(2,2,3); pbaspect([1.5 1 1]); hold on;
    xPlot = [histMids_fine(1); histMids_fine; histMids_fine(end)];    
    fill(xPlot,[0; endProbFine(:,1);0],'k','FaceAlpha',0.8,'EdgeColor','none');
    fill(xPlot,[0; endProbFine(:,2);0],'r','FaceAlpha',0.5,'EdgeColor','none');
    fill(xPlot,[0; endProbFine(:,3);0],'m','FaceAlpha',0.5,'EdgeColor','none');
    fill(xPlot,[0; endProbFine(:,4);0],'b','FaceAlpha',0.5,'EdgeColor','none');
    fill(xPlot,[0; endProbFine(:,5);0],'k','FaceAlpha',0.3,'EdgeColor','none'); 
    xline([-1/3 1/3],'--r');
    xlim([-1.1 1.1]); xticks([-1:1/3:1]); xticklabels({'-1','-2/3','-1/3','1','1/3','2/3','1'}); xtickangle(0)
    ylim([0 15]); yticks([0:3:15]);
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out'); box off;
    title('Fine binning');    
    xlabel('Endpoint position');
    
    endTrace = seshOut(iSesh).endTrace.endTraceMu;
    endTraceSEM = seshOut(iSesh).endTrace.endTraceSEM;
    endTraceReal = seshOut(iSesh).endTrace.endTraceRealMu;
    endTraceRealSEM = seshOut(iSesh).endTrace.endTraceRealSEM;
    endVel = seshOut(iSesh).endTrace.movVelMu;
    endVelSEM = seshOut(iSesh).endTrace.movVelSEM;
    
    subplot(2,2,2); pbaspect([1.5 1 1]); hold on;
    tPlot = (1:size(endTrace,1))'/100-1/50;
    rectangle('Position',[0 -1/3 ceil(tPlot(end)) 2/3],'FaceColor',[0 .5 .5 .3],'EdgeColor','none');
    h1 = boundedline(tPlot,endTrace,endTraceSEM,'cmap',[0 0 0],'alpha');
    set(h1,'linewidth',1.5);
    h2 = boundedline(tPlot,endTraceReal,endTraceRealSEM,'cmap',[.5 0 0],'alpha');
    set(h2,'linewidth',1.5);
    xline([0:2:ceil(tPlot(end))-1],'--k');
    ylim([-1.05 1.05]); yticks([-1:1/3:1]); yticklabels({'-1','-2/3','-1/3','1','1/3','2/3','1'})
    xlim([0 ceil(tPlot(end))]); xticks([1:2:ceil(tPlot(end))]); xticklabels({});
    title('Wheel displacement traces');
    ylabel('Object/wheel endpoint');
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
    
    subplot(2,2,4); pbaspect([1.5 1 1]); hold on;    
%     rectangle('Position',[0 -1/3 ceil(tPlot(end)) 2/3],'FaceColor',[0 .5 .5 .3],'EdgeColor','none');
    h1 = boundedline(tPlot,endVel,endVelSEM,'cmap',[0 0 0],'alpha');
    set(h1,'linewidth',1.5);
    xline([0:2:ceil(tPlot(end))-1],'--k');
    ylim([-2.5 15]); yticks([-2.5:2.5:15]);
    xlim([0 ceil(tPlot(end))]); xticks([1:2:ceil(tPlot(end))]); xticklabels({'Baseline','Early Adapt.','Mid Adapt','Late Adapt.','Washout'});
    title('Wheel velocity traces');
    ylabel('Velocity (cm/s)');
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
    
    suptitle(sprintf('%s on %s: Adaptation endpoints',seshOut(iSesh).name,seshOut(iSesh).date));
    
    if savebool
        figFold = fullfile(paperFold,'Figures','Fig1','Sub1_v5Adaptation');
        if ~exist(figFold,'dir'); mkdir(figFold); end
        figname = sprintf('%s_%s_endProb',seshOut(iSesh).name,seshOut(iSesh).date);
        savefig(fHist, fullfile(figFold,[figname,'.fig']));
        saveas(fHist,fullfile(figFold,[figname,'.png']));
        print(fHist,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');
    end
end

%% Plot group data
% Percent correct and endpoint positions
endPos_grp = 1-cat(2,seshOut.endPosMu);
endPosReal_grp = 1-cat(2,seshOut.endPosRealMu);
percC_grp = cat(2,seshOut.percC);
endTrace_grp = zeros(numel(seshOut(1).endTrace.endTraceMu),numel(seshOut));
endTraceReal_grp = endTrace_grp;
for iSesh = 1:numel(seshOut)
    endTrace_grp(:,iSesh) = seshOut(iSesh).endTrace.endTraceMu;
    endTraceReal_grp(:,iSesh) = seshOut(iSesh).endTrace.endTraceRealMu;
end

[percC_P, percC_Table, percC_Stats] = kruskalwallis(percC_grp);
percC_Stats_post = multcompare(percC_Stats, 'Ctype','bonferroni');
statsAll = v2struct(percC_P,percC_Table,percC_Stats, percC_Stats_post);

fBar = figure;
subplot(2,2,1); pbaspect([1.5 1 1]); hold on;
rectangle('Position',[.5 2/3 5.5 2/3],'FaceColor',[0 .5 .5 .3],'EdgeColor','none');
iRow = 1;
bar(iRow,mean(endPos_grp(iRow,:),2,'omitnan'),0.5,'EdgeColor','none','FaceColor',[0 0 0],'FaceAlpha',0.8);
errorbar(iRow, mean(endPos_grp(iRow,:),2,'omitnan'), std(endPos_grp(iRow,:),[],2,'omitnan')/sqrt(numel(seshOut)),'k','CapSize',0,'LineWidth',2);
iRow = iRow+1;
bar(iRow,mean(endPos_grp(iRow,:),2,'omitnan'),0.5,'EdgeColor','none','FaceColor',[.5 0 0],'FaceAlpha',0.5);
errorbar(iRow, mean(endPos_grp(iRow,:),2,'omitnan'), std(endPos_grp(iRow,:),[],2,'omitnan')/sqrt(numel(seshOut)),'k','CapSize',0,'LineWidth',2);
iRow = iRow+1;
bar(iRow,mean(endPos_grp(iRow,:),2,'omitnan'),0.5,'EdgeColor','none','FaceColor',[.5 0 .5],'FaceAlpha',0.5);
errorbar(iRow, mean(endPos_grp(iRow,:),2,'omitnan'), std(endPos_grp(iRow,:),[],2,'omitnan')/sqrt(numel(seshOut)),'k','CapSize',0,'LineWidth',2);
iRow = iRow+1;
bar(iRow,mean(endPos_grp(iRow,:),2,'omitnan'),0.5,'EdgeColor','none','FaceColor',[0 0 .5],'FaceAlpha',0.5);
errorbar(iRow, mean(endPos_grp(iRow,:),2,'omitnan'), std(endPos_grp(iRow,:),[],2,'omitnan')/sqrt(numel(seshOut)),'k','CapSize',0,'LineWidth',2);
iRow = iRow+1;
bar(iRow,mean(endPos_grp(iRow,:),2,'omitnan'),0.5,'EdgeColor','none','FaceColor',[0 0 0],'FaceAlpha',0.3);
errorbar(iRow, mean(endPos_grp(iRow,:),2,'omitnan'), std(endPos_grp(iRow,:),[],2,'omitnan')/sqrt(numel(seshOut)),'k','CapSize',0,'LineWidth',2);
ylim([0 2]); yticks([0:1/3:2]);yticklabels({'0','1/3','2/3','1','4/3','5/3','2'})
xlim([0.5 5.5]); xticks([1:5]); xticklabels({}); %xticklabels({'Trained day 1','Trained day N','Adapt. baseline',sprintf('Adapt. early (%s)',num2str(adaptTrialGrp)),'Adapt. late','Adapt. washout'});
ylabel('Object endpoint');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
title('Object displacements');

subplot(2,2,3); pbaspect([1.5 1 1]); hold on;
rectangle('Position',[.5 2/3 5.5 2/3],'FaceColor',[0 .5 .5 .3],'EdgeColor','none');
iRow = 1;
bar(iRow,mean(endPosReal_grp(iRow,:),2,'omitnan'),0.5,'EdgeColor','none','FaceColor',[0 0 0],'FaceAlpha',0.8);
errorbar(iRow, mean(endPosReal_grp(iRow,:),2,'omitnan'), std(endPosReal_grp(iRow,:),[],2,'omitnan')/sqrt(numel(seshOut)),'k','CapSize',0,'LineWidth',2);
iRow = iRow+1;
bar(iRow,mean(endPosReal_grp(iRow,:),2,'omitnan'),0.5,'EdgeColor','none','FaceColor',[.5 0 0],'FaceAlpha',0.5);
errorbar(iRow, mean(endPosReal_grp(iRow,:),2,'omitnan'), std(endPosReal_grp(iRow,:),[],2,'omitnan')/sqrt(numel(seshOut)),'k','CapSize',0,'LineWidth',2);
iRow = iRow+1;
bar(iRow,mean(endPosReal_grp(iRow,:),2,'omitnan'),0.5,'EdgeColor','none','FaceColor',[.5 0 .5],'FaceAlpha',0.5);
errorbar(iRow, mean(endPosReal_grp(iRow,:),2,'omitnan'), std(efndPosReal_grp(iRow,:),[],2,'omitnan')/sqrt(numel(seshOut)),'k','CapSize',0,'LineWidth',2);
iRow = iRow+1;
bar(iRow,mean(endPosReal_grp(iRow,:),2,'omitnan'),0.5,'EdgeColor','none','FaceColor',[0 0 .5],'FaceAlpha',0.5);
errorbar(iRow, mean(endPosReal_grp(iRow,:),2,'omitnan'), std(endPosReal_grp(iRow,:),[],2,'omitnan')/sqrt(numel(seshOut)),'k','CapSize',0,'LineWidth',2);
iRow = iRow+1;
bar(iRow,mean(endPosReal_grp(iRow,:),2,'omitnan'),0.5,'EdgeColor','none','FaceColor',[0 0 0],'FaceAlpha',0.3);
errorbar(iRow, mean(endPosReal_grp(iRow,:),2,'omitnan'), std(endPosReal_grp(iRow,:),[],2,'omitnan')/sqrt(numel(seshOut)),'k','CapSize',0,'LineWidth',2);
ylim([0 2]); yticks([0:1/3:2]);yticklabels({'0','1/3','2/3','1','4/3','5/3','2'})
xlim([0.5 5.5]); xticks([1:5]); xticklabels({'Baseline','Early Adapt.','Mid Adapt','Late Adapt.','Washout'});
ylabel('Object endpoint');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
title('Wheel displacements');

subplot(2,2,2); pbaspect([1.5 1 1]); hold on;
iRow = 1;
bar(iRow,mean(percC_grp(iRow,:),2,'omitnan'),0.5,'EdgeColor','none','FaceColor',[0 0 0],'FaceAlpha',0.8);
errorbar(iRow, mean(percC_grp(iRow,:),2,'omitnan'), std(percC_grp(iRow,:),[],2,'omitnan')/sqrt(numel(seshOut)),'k','CapSize',0,'LineWidth',2);
iRow = iRow+1;
bar(iRow,mean(percC_grp(iRow,:),2,'omitnan'),0.5,'EdgeColor','none','FaceColor',[.5 0 0],'FaceAlpha',0.5);
errorbar(iRow, mean(percC_grp(iRow,:),2,'omitnan'), std(percC_grp(iRow,:),[],2,'omitnan')/sqrt(numel(seshOut)),'k','CapSize',0,'LineWidth',2);
iRow = iRow+1;
bar(iRow,mean(percC_grp(iRow,:),2,'omitnan'),0.5,'EdgeColor','none','FaceColor',[.5 0 .5],'FaceAlpha',0.5);
errorbar(iRow, mean(percC_grp(iRow,:),2,'omitnan'), std(percC_grp(iRow,:),[],2,'omitnan')/sqrt(numel(seshOut)),'k','CapSize',0,'LineWidth',2);
iRow = iRow+1;
bar(iRow,mean(percC_grp(iRow,:),2,'omitnan'),0.5,'EdgeColor','none','FaceColor',[0 0 .5],'FaceAlpha',0.5);
errorbar(iRow, mean(percC_grp(iRow,:),2,'omitnan'), std(percC_grp(iRow,:),[],2,'omitnan')/sqrt(numel(seshOut)),'k','CapSize',0,'LineWidth',2);
iRow = iRow+1;
bar(iRow,mean(percC_grp(iRow,:),2,'omitnan'),0.5,'EdgeColor','none','FaceColor',[0 0 0],'FaceAlpha',0.3);
errorbar(iRow, mean(percC_grp(iRow,:),2,'omitnan'), std(percC_grp(iRow,:),[],2,'omitnan')/sqrt(numel(seshOut)),'k','CapSize',0,'LineWidth',2);
ylim([15 65]); yticks([0:5:100]);
xlim([0.5 5.5]); xticks([1:5]); xticklabels({}); %xticklabels({'Baseline','Early Adapt.','Mid Adapt','Late Adapt.','Washout'});
ylabel('Percent correct');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');


subplot(2,2,4); pbaspect([1.5 1 1]); hold on;
tPlot = (1:size(endTrace_grp,1))'/100-1/50;
rectangle('Position',[0 -1/3 ceil(tPlot(end)) 2/3],'FaceColor',[0 .5 .5 .3],'EdgeColor','none');
h1 = boundedline(tPlot,mean(endTrace_grp,2,'omitnan'),std(endTrace_grp,[],2,'omitnan')/sqrt(numel(seshOut)),'cmap',[0 0 0],'alpha');
set(h1,'linewidth',1.5);
h2 = boundedline(tPlot,mean(endTraceReal_grp,2,'omitnan'),std(endTraceReal_grp,[],2,'omitnan')/sqrt(numel(seshOut)),'cmap',[.5 0 0],'alpha');
set(h2,'linewidth',1.5);
xline([0:2:ceil(tPlot(end))-1],'--k');
ylim([-1.05 1.05]); yticks([-1:1/3:1]);yticklabels({'-1','-2/3','-1/3','1','1/3','2/3','1'})
xlim([0 ceil(tPlot(end))]); xticks([1:2:ceil(tPlot(end))]); xticklabels({'Baseline','Early Adapt.','Mid Adapt','Late Adapt.','Washout'});
title('Wheel displacement traces');
ylabel('Object/wheel endpoint');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');



%% Endpoint histograms
fHistAll = figure; 
subplot(1,2,1); pbaspect([1.5 1 1]); hold on;
rectangle('Position',[-1/3 0 2/3 15],'FaceColor',[0 .5 .5 .3],'EdgeColor','none');
h1 = boundedline(histMids_coarse,mean(endProbAll_coarse(:,1,:),3,'omitnan'),std(endProbAll_coarse(:,1,:),[],3,'omitnan')/size(endProbAll_coarse,3),...
    'cmap',[0 0 0],'alpha'); set(h1,'linewidth',1.5);
h2 = boundedline(histMids_coarse,mean(endProbAll_coarse(:,2,:),3,'omitnan'),std(endProbAll_coarse(:,2,:),[],3,'omitnan')/size(endProbAll_coarse,3),...
    'cmap',[.5 0 0],'alpha'); set(h2,'linewidth',1.5);
h3 = boundedline(histMids_coarse,mean(endProbAll_coarse(:,3,:),3,'omitnan'),std(endProbAll_coarse(:,3,:),[],3,'omitnan')/size(endProbAll_coarse,3),...
    'cmap',[.5 0 .5],'alpha'); set(h3,'linewidth',1.5);
h4 = boundedline(histMids_coarse,mean(endProbAll_coarse(:,4,:),3,'omitnan'),std(endProbAll_coarse(:,4,:),[],3,'omitnan')/size(endProbAll_coarse,3),...
    'cmap',[0 0 .5],'alpha'); set(h4,'linewidth',1.5);
h5 = boundedline(histMids_coarse,mean(endProbAll_coarse(:,5,:),3,'omitnan'),std(endProbAll_coarse(:,5,:),[],3,'omitnan')/size(endProbAll_coarse,3),...
    'cmap',[.5 .5 .5],'alpha'); set(h5,'linewidth',1.5);
% xline([-1/3 1/3],'--r')
xlim([-1.1 1.1]); xticks([-1:1/3:1]); xtickangle(0)
ylim([0 15]); yticks([0:3:15]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out'); box off;
title('Coarse lines');
ylabel('Probability (percent)');
xlabel('Endpoint position');

subplot(1,2,2); pbaspect([1.5 1 1]); hold on;
rectangle('Position',[-1/3 0 2/3 15],'FaceColor',[0 .5 .5 .3],'EdgeColor','none');
h1 = boundedline(histMids_fine,mean(endProbAll_fine(:,1,:),3,'omitnan'),std(endProbAll_fine(:,1,:),[],3,'omitnan')/size(endProbAll_coarse,3),...
    'cmap',[0 0 0],'alpha'); set(h1,'linewidth',1.5);
h2 = boundedline(histMids_fine,mean(endProbAll_fine(:,2,:),3,'omitnan'),std(endProbAll_fine(:,2,:),[],3,'omitnan')/size(endProbAll_coarse,3),...
    'cmap',[.5 0 0],'alpha'); set(h2,'linewidth',1.5);
h3 = boundedline(histMids_fine,mean(endProbAll_fine(:,3,:),3,'omitnan'),std(endProbAll_fine(:,3,:),[],3,'omitnan')/size(endProbAll_coarse,3),...
    'cmap',[.5 0 .5],'alpha'); set(h3,'linewidth',1.5);
h4 = boundedline(histMids_fine,mean(endProbAll_fine(:,4,:),3,'omitnan'),std(endProbAll_fine(:,4,:),[],3,'omitnan')/size(endProbAll_coarse,3),...
    'cmap',[0 0 .5],'alpha'); set(h4,'linewidth',1.5);
h5 = boundedline(histMids_fine,mean(endProbAll_fine(:,5,:),3,'omitnan'),std(endProbAll_fine(:,5,:),[],3,'omitnan')/size(endProbAll_coarse,3),...
    'cmap',[.5 .5 .5],'alpha'); set(h5,'linewidth',1.5);
% xline([-1/3 1/3],'--r');
xlim([-1.1 1.1]); xticks([-1:1/3:1]); xtickangle(0)
ylim([0 15]); yticks([0:3:15]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out'); box off;
title('Fine lines');
ylabel('Probability (percent)');
xlabel('Endpoint position');


%% Save maybe

if savebool
    figFold = fullfile(paperFold,'Figures','Fig1','Sub1_v5Adaptation');
    if ~exist(figFold,'dir'); mkdir(figFold); end
    
%     figname = 'v5LearningCurves_AllTrials';
%     savefig(f1, fullfile(figFold,[figname,'.fig']));
%     saveas(f1,fullfile(figFold,[figname,'.png']));
%     print(f1,fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    figname = 'v5SummaryBars_GoodTrials';
    savefig(fBar, fullfile(figFold,[figname,'.fig']));
    saveas(fBar,fullfile(figFold,[figname,'.png']));
    print(fBar,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');
    figname = 'v5AdaptationCurves_endProb';
    savefig(fHistAll, fullfile(figFold,[figname,'.fig']));
    saveas(fHistAll,fullfile(figFold,[figname,'.png']));
    print(fHistAll,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');
    sessionsv5 = seshOut;
    save(fullfile(figFold,'vAdaptationCurves.mat'),'sessionsv5','endProbAll_coarse','endProbAll_fine','histBins_coarse','histBins_fine', 'statsAll');

end

end

%%

%{
%% Version 5 section
if ~exist(fullfile(figFold,'percMovCorrect_v5.mat'),'file')
% if ~exist(fullfile(figFold,['percMovCorrect_v5_',num2str(adaptTrialGrp),'trials.mat']),'file')
% If this has not already been done, iterate over session and extract 
% percent correct and fraction of trials with movement
for i = 1:numel(sessionsv5)
    dataDir = dir(fullfile(behavFold,sessionsv5(i).name,sessionsv5(i).date,'Behav*',[sessionsv5(i).name,'_',sessionsv5(i).date,'.mat']));
    datastruct = load(fullfile(dataDir.folder,dataDir.name)); datastruct = datastruct.datastruct;
    trials = datastruct.paqdata.behav.trials;
    inits = datastruct.inits;
    sets = Rot2p_changetrials(datastruct.sets);
    gaindiv = 37.5/360*4.5*pi/sets.rotgain;
    if sets.vrfs == 60
        smspan = inits.paqfs*.04; % 40 ms smoothing span
    elseif sets. vrfs == 30
        smspan = inits.paqfs*.06; % 60 ms smoothing span
    end
    if strcmpi(sessionsv5(i).name,'dk105') & strcmp(sessionsv5(i).date,'2018-03-07')
        sets.changebacktrial = sets.changetrial+10*(sets.changebackmult-1);
    end
    trials = trials(1:find(cat(1,trials.realtrial) > 0,1, 'last'));
    
    behavstability = Rot2p_behavstability(trials,inits,sets,0);
    
    sessionsv5(i).nTrials = numel(trials);
    sessionsv5(i).nTrialsReal = sum(cat(1,trials.realtrial));
    sessionsv5(i).nTrialsBase = sum(cat(1,trials(1:sets.changetrial).realtrial));
    sessionsv5(i).nTrialsEarly = sum(cat(1,trials(sets.changetrial+1:sets.changetrial+adaptTrialGrp).realtrial));
    sessionsv5(i).nTrialsLate = sum(cat(1,trials(sets.changetrial+121-adaptTrialGrp:sets.changetrial+120).realtrial));
    
    if ~sets.changeback
        sessionsv5(i).nTrialsWash = 0;
    else
        nWashTrials = numel(trials)-sets.changebacktrial;
        if nWashTrials > adaptTrialGrp; nWashTrials = adaptTrialGrp; end
        sessionsv5(i).nTrialsWash = sum(cat(1,trials(sets.changebacktrial+1:sets.changebacktrial+nWashTrials).realtrial));
    end
    
    sessionsv5(i).iTrialsReal = false(sessionsv5(i).nTrials,1);
    sessionsv5(i).movBool = false(sessionsv5(i).nTrials,1);
    sessionsv5(i).corrBool = false(sessionsv5(i).nTrials,1);

    % Determine if movement happened and if it resulted in correct trial
    for j = 1:numel(trials)
        x0 = trials(j).data(inits.paqfs*inits.mintime+1:end-inits.paqfs*inits.mintime,10);
        x1 = smooth(x0,smspan/numel(x0),'lowess');
        v1 = diff(x1)*inits.paqfs*-gaindiv; v1(end+1) = v1(end);
        if numel(v1) > 15*inits.paqfs; v1 = v1(1:15*inits.paqfs); end
        if numel(v1)/inits.paqfs <= sets.responsetime
            sessionsv5(i).iTrialsReal(j) = true;
        end
        if any(abs(v1) > 1.5) % movement threshold is 1.5 cm/s of wheel translation
%         if any(posTrace > sets.rotgain*360/45/pi()) % movement threshold is 1 mm of wheel translation
            sessionsv5(i).movBool(j) = true;
        end
        if sessionsv5(i).iTrialsReal(j) & trials(j).correct
            sessionsv5(i).corrBool(j) = true;
        end
    end
    sessionsv5(i).movPerc = sum(sessionsv5(i).movBool)/sessionsv5(i).nTrials;
    if sessionsv5(i).nTrialsReal > 0
        sessionsv5(i).corrPerc = sum(sessionsv5(i).corrBool)/sessionsv5(i).nTrialsReal;
        sessionsv5(i).corrPercBase = sum(sessionsv5(i).corrBool(1:sets.changetrial))/sessionsv5(i).nTrialsBase;
        sessionsv5(i).corrPercEarly = sum(sessionsv5(i).corrBool(sets.changetrial+1:sets.changetrial+adaptTrialGrp))/sessionsv5(i).nTrialsEarly;
        sessionsv5(i).corrPercLate = sum(sessionsv5(i).corrBool(sets.changetrial+121-adaptTrialGrp:sets.changetrial+120))/sessionsv5(i).nTrialsLate;
        if ~sets.changeback
            sessionsv5(i).corrPercWash = nan;
        else
            sessionsv5(i).corrPercWash = sum(sessionsv5(i).corrBool(sets.changebacktrial+1:sets.changebacktrial+nWashTrials))/sessionsv5(i).nTrialsWash;
        end
        xEndAll = nan(size(sessionsv5(i).iTrialsReal)); xEndAll(behavstability.params.trialinds) = behavstability.params.xendall;
        
        sessionsv5(i).xEndMu_Base = nanmean(xEndAll(1:sets.changetrial));
        sessionsv5(i).xEndMed_Base = nanmedian(xEndAll(1:sets.changetrial));
        sessionsv5(i).xEndMu_Early = nanmean(xEndAll(sets.changetrial+1:sets.changetrial+adaptTrialGrp));
        sessionsv5(i).xEndMed_Early = nanmedian(xEndAll(sets.changetrial+1:sets.changetrial+adaptTrialGrp));
        sessionsv5(i).xEndMu_Late = nanmean(xEndAll(sets.changetrial+121-adaptTrialGrp:sets.changetrial+120));
        sessionsv5(i).xEndMed_Late = nanmedian(xEndAll(sets.changetrial+121-adaptTrialGrp:sets.changetrial+120));
        if ~sets.changeback
            sessionsv5(i).xEndMu_Wash = nan;
            sessionsv5(i).xEndMed_Wash = nan;
        else
            sessionsv5(i).xEndMu_Wash = nanmean(xEndAll(sets.changebacktrial+1:sets.changebacktrial+nWashTrials));
            sessionsv5(i).xEndMed_Wash = nanmedian(xEndAll(sets.changebacktrial+1:sets.changebacktrial+nWashTrials));
        end
    else
        sessionsv5(i).corrPerc = 0;
    end
    sessionsv5(i).sets = sets;
    sessionsv5(i).params = rmfield(behavstability.params,{'trialnum','blocksize','rhoU','rhoC','rhoUmu','rhoUsd','rhoCmu','rhoCsd','rhoOmu','rhoOsd','numC'});
    sessionsv5(i).offXall = behavstability.offxall;
    sessionsv5(i).movVall = behavstability.movvall;
    
    switch lower(sessionsv5(i).name)
        case 'dk103'
            sessionsv5(i).matchColor = matchColors(1,:);
        case 'dk105'
            sessionsv5(i).matchColor = matchColors(2,:);
        case 'dk169'
            sessionsv5(i).matchColor = matchColors(3,:);
        case 'dk171'
            sessionsv5(i).matchColor = matchColors(4,:);
        case 'dk194'
            sessionsv5(i).matchColor = matchColors(5,:);
        case 'dk199'
            sessionsv5(i).matchColor = matchColors(6,:);
        case 'dk052'
            sessionsv5(i).matchColor = matchColors(7,:);
        case 'dk063'
            sessionsv5(i).matchColor = matchColors(8,:);
        case 'dk070'
            sessionsv5(i).matchColor = matchColors(9,:);
        case 'dk138'
            sessionsv5(i).matchColor = matchColors(10,:);
    end
end


if savebool
    if ~exist(figFold,'dir'); mkdir(figFold); end
    save(fullfile(figFold,'percMovCorrect_v5.mat'),'sessionsv5');
%     save(fullfile(figFold,['percMovCorrect_v5_',num2str(adaptTrialGrp),'trials.mat']),'sessionsv5');
    if ~exist(figFold2,'dir'); mkdir(figFold2); end
    save(fullfile(figFold2,'percMovCorrect_v5.mat'),'sessionsv5');
%     save(fullfile(figFold2,['percMovCorrect_v5_',num2str(adaptTrialGrp),'trials.mat']),'sessionsv5');
end
    
else
    % If this has already been done, just load the file to save time
    sessionsv5 = load(fullfile(figFold,'percMovCorrect_v5.mat'),'sessionsv5');
%     sessionsv5 = load(fullfile(figFold,['percMovCorrect_v5_',num2str(adaptTrialGrp),'trials.mat']));
    sessionsv5 = sessionsv5.sessionsv5;
    
    for i = 1:numel(sessionsv5)
        sets = sessionsv5(i).sets;
        
        sessionsv5(i).nTrialsBase = sum(sessionsv5(i).iTrialsReal(1:sets.changetrial));
        sessionsv5(i).nTrialsEarly = sum(sessionsv5(i).iTrialsReal(sets.changetrial+1:sets.changetrial+adaptTrialGrp));
        sessionsv5(i).nTrialsLate = sum(sessionsv5(i).iTrialsReal(sets.changetrial+121-adaptTrialGrp:sets.changetrial+120));
%         sessionsv5(i).nTrialsBase = sum(cat(1,trials(1:sets.changetrial).realtrial));
%         sessionsv5(i).nTrialsEarly = sum(cat(1,trials(sets.changetrial+1:sets.changetrial+adaptTrialGrp).realtrial));
%         sessionsv5(i).nTrialsLate = sum(cat(1,trials(sets.changetrial+121-adaptTrialGrp:sets.changetrial+120).realtrial));
        
        if ~sets.changeback
            sessionsv5(i).nTrialsWash = 0;
        else
            nWashTrials = sessionsv5(i).nTrials-sets.changebacktrial;
            if nWashTrials > adaptTrialGrp; nWashTrials = adaptTrialGrp; end
            sessionsv5(i).nTrialsWash = sum(sessionsv5(i).iTrialsReal(sets.changebacktrial+1:sets.changebacktrial+nWashTrials));
%             sessionsv5(i).nTrialsWash = sum(cat(1,trials(sets.changebacktrial+1:sets.changebacktrial+nWashTrials).realtrial));
        end
        
        if sessionsv5(i).nTrialsReal > 0
            sessionsv5(i).corrPerc = sum(sessionsv5(i).corrBool)/sessionsv5(i).nTrialsReal;
            sessionsv5(i).corrPercBase = sum(sessionsv5(i).corrBool(1:sets.changetrial))/sessionsv5(i).nTrialsBase;
            sessionsv5(i).corrPercEarly = sum(sessionsv5(i).corrBool(sets.changetrial+1:sets.changetrial+adaptTrialGrp))/sessionsv5(i).nTrialsEarly;
            sessionsv5(i).corrPercLate = sum(sessionsv5(i).corrBool(sets.changetrial+121-adaptTrialGrp:sets.changetrial+120))/sessionsv5(i).nTrialsLate;
            if ~sets.changeback
                sessionsv5(i).corrPercWash = nan;
            else
                sessionsv5(i).corrPercWash = sum(sessionsv5(i).corrBool(sets.changebacktrial+1:sets.changebacktrial+nWashTrials))/sessionsv5(i).nTrialsWash;
            end
            xEndAll = nan(size(sessionsv5(i).iTrialsReal)); xEndAll(sessionsv5(i).params.trialinds) = sessionsv5(i).params.xendall;
            
            sessionsv5(i).xEndMu_Base = nanmean(xEndAll(1:sets.changetrial));
            sessionsv5(i).xEndMed_Base = nanmedian(xEndAll(1:sets.changetrial));
            sessionsv5(i).xEndMu_Early = nanmean(xEndAll(sets.changetrial+1:sets.changetrial+adaptTrialGrp));
            sessionsv5(i).xEndMed_Early = nanmedian(xEndAll(sets.changetrial+1:sets.changetrial+adaptTrialGrp));
            sessionsv5(i).xEndMu_Late = nanmean(xEndAll(sets.changetrial+121-adaptTrialGrp:sets.changetrial+120));
            sessionsv5(i).xEndMed_Late = nanmedian(xEndAll(sets.changetrial+121-adaptTrialGrp:sets.changetrial+120));
            
            if ~sets.changeback
                sessionsv5(i).xEndMu_Wash = nan;
                sessionsv5(i).xEndMed_Wash = nan;
            else
                sessionsv5(i).xEndMu_Wash = nanmean(xEndAll(sets.changebacktrial+1:sets.changebacktrial+nWashTrials));
                sessionsv5(i).xEndMed_Wash = nanmedian(xEndAll(sets.changebacktrial+1:sets.changebacktrial+nWashTrials));
            end
            
        else
            sessionsv5(i).corrPerc = 0;
        end
    end
    
end

clear datastruct gaindiv i inits j nWashTrials sets smspan trials v1 x0 x1 xEndAll
%}
