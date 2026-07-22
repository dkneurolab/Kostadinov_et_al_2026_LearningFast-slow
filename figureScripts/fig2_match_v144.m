function matchPCs_v144 = fig2_match_v144(matchData,matchSpont,matchFOVs0,dataFold,paperFold,matchColors,savebool)
%% Set sampling stuff - assume second dataset of first recording is v4
fs.imfs = matchSpont(1).sets2.imfs;
fs.behavfs = matchSpont(1).sets2.fs/50;
fs.rewDelay = matchSpont(1).sets2.rewdelay;

%% Get sessions that have v145
a = 0;
for i = 1:numel(matchFOVs0)
    if ~isempty(matchFOVs0(i).v1date) && ~isempty(matchFOVs0(i).v4_1date)&& ~isempty(matchFOVs0(i).v4_Ndate)
        a = a+1;
        matchFOVs(a) = matchFOVs0(i);
    end
end
matchFOVs = rmfield(matchFOVs,{'iv1sesh','iv4_1sesh','iv4_Nsesh','iv5sesh','v5date'});

%% Get matching v145 data in one structure

for i = 1:numel(matchFOVs)
    for j = 1:numel(matchData)
        if strcmpi(matchFOVs(i).name,matchData(j).name)
            if strcmpi(matchFOVs(i).v1date,matchData(j).date1) && ...
                    strcmpi(matchFOVs(i).v4_Ndate,matchData(j).date2)
                matchFOVs(i).col = matchData(j).matchColor;
                matchFOVs(i).v14rois = [matchData(j).iROIs1 matchData(j).iROIs2];
                v14on = [matchData(j).on1 matchData(j).on2];
                v14mov = [matchData(j).mov1 matchData(j).mov2];
                v14rewT = [matchData(j).rewT1 matchData(j).rewT2];
                v14rewC = [matchData(j).rewC1 matchData(j).rewC2];
                v14rewR = [matchData(j).rewR1 matchData(j).rewR2];
                v14spkRate = {matchSpont(j).spkRate1; matchSpont(j).spkRate2};
                v14corrMat = {corrcoef(matchSpont(j).zF1All); corrcoef(matchSpont(j).zF2All)};
            elseif strcmpi(matchFOVs(i).v4_1date,matchData(j).date1) && ...
                    strcmpi(matchFOVs(i).v4_Ndate,matchData(j).date2)
                matchFOVs(i).v44rois = [matchData(j).iROIs1 matchData(j).iROIs2];
                v44on = [matchData(j).on1 matchData(j).on2];
                v44mov = [matchData(j).mov1 matchData(j).mov2];
                v44rewT = [matchData(j).rewT1 matchData(j).rewT2];
                v44rewC = [matchData(j).rewC1 matchData(j).rewC2];
                v44rewR = [matchData(j).rewR1 matchData(j).rewR2];
                v44spkRate = {matchSpont(j).spkRate1; matchSpont(j).spkRate2};
                v44corrMat = {corrcoef(matchSpont(j).zF1All); corrcoef(matchSpont(j).zF2All)};
            end
        end
    end
    
    a = ismember(matchFOVs(i).v14rois(:,2),matchFOVs(i).v44rois(:,2));
    a2 = ismember(matchFOVs(i).v44rois(:,2),matchFOVs(i).v14rois(:,2));
    
    
    v14rois_match0 = matchFOVs(i).v14rois(a,:);
    v44rois_match0 = matchFOVs(i).v44rois(a2,:);
    
    [~, v14order] = sort(v14rois_match0(:,2),'ascend');
    [~, v44order] = sort(v44rois_match0(:,2),'ascend');
    v14rois_match = v14rois_match0(v14order,:);
    v44rois_match = v44rois_match0(v44order,:);
    
    matchFOVs(i).v144rois = [v14rois_match(:,1) v44rois_match];
    matchFOVs(i).v144on = cat(2,v14on(1),v44on);
    matchFOVs(i).v144mov = cat(2,v14mov(1),v44mov);
    matchFOVs(i).v144rewT = cat(2,v14rewT(1),v44rewT);
    matchFOVs(i).v144rewC = cat(2,v14rewC(1),v44rewC);
    matchFOVs(i).v144rewR = cat(2,v14rewR(1),v44rewR);
    matchFOVs(i).v144spkRate = cat(1,v14spkRate(1),v44spkRate);
    matchFOVs(i).v144corrMat = cat(1,v14corrMat(1),v44corrMat);
    
    clear v14on v14mov v14rewT v14rewC v14rewR v14spkRate v14corrMat v44on v44mov v44rewT v44rewC v44rewR v44spkRate v44corrMat 
    clear v14rois_match v14rois_match0 v14order v45order v45rois_match v45rois_match0 a a2
end

clear matchData

%% Process data by location

[lobvfigs,lobvPCs] = fig1_match_v144plot(matchFOVs(cat(1,matchFOVs.ilobv)),'Lobule V',fs);
[sim2figs,sim2PCs] = fig1_match_v144plot(matchFOVs(cat(1,matchFOVs.isim2)),'Simplex',fs);

lobvPCs = fig1_match_v144plotSpont(lobvPCs, dataFold, paperFold, savebool);
sim2PCs = fig1_match_v144plotSpont(sim2PCs, dataFold, paperFold, savebool);

lobvFRcat = cat(1,lobvPCs.v144spkRate);
nLobv = size(lobvFRcat,1);
for i = 1:numel(lobvPCs)
    lobvFR(i,:) = mean(lobvPCs(i).v144spkRate);
    lobvCorr(i,:) = lobvPCs(i).corrPlot(1,:);
end

sim2FRcat = cat(1,sim2PCs.v144spkRate);
nSim2 = size(sim2FRcat,1);
for j = 1:numel(sim2PCs)
    sim2FR(j,:) = mean(sim2PCs(j).v144spkRate);
    sim2Corr(j,:) = sim2PCs(j).corrPlot(1,:);
end

figCorrFR = figure;
subplot(2,3,1); hold on;
bar([1 4 7],mean(lobvFRcat),0.3,'EdgeColor','none','FaceColor',[.6 .6 .6],'FaceAlpha',0.5);
errorbar([1 4 7], mean(lobvFRcat), std(lobvFRcat),'k','CapSize',0,'LineWidth',2,'LineStyle','none');
plot([2*ones(1,nLobv); 3*ones(1,nLobv)], [lobvFRcat(:,1) lobvFRcat(:,2)]','LineWidth',1,'Color',[0 0 0 0.2]);
plot([5*ones(1,nLobv); 6*ones(1,nLobv)], [lobvFRcat(:,2) lobvFRcat(:,3)]','LineWidth',1,'Color',[0 0 0 0.2]);
ylim([0 3]); yticks([0:0.5:3]);
xlim([0 8]); xticks([1,4,7]); xticklabels({'Naive','Trained','Adapt.'})
ylabel('LobV PC Firing rate (spks/s)');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
axis square;

subplot(2,3,2); hold on;
bar([1 4 7],mean(lobvFR),0.3,'EdgeColor','none','FaceColor',[.6 .6 .6],'FaceAlpha',0.5);
errorbar([1 4 7], mean(lobvFR), std(lobvFR),'k','CapSize',0,'LineWidth',2,'LineStyle','none');
plot([2*ones(1,i); 3*ones(1,i)], [lobvFR(:,1) lobvFR(:,2)]','LineWidth',1,'Color',[0 0 0 0.2]);
plot([5*ones(1,i); 6*ones(1,i)], [lobvFR(:,2) lobvFR(:,3)]','LineWidth',1,'Color',[0 0 0 0.2]);
ylim([0 3]); yticks([0:0.5:3]);
xlim([0 8]); xticks([1,4,7]); xticklabels({'Naive','Trained','Adapt.'})
ylabel('LobV FOV Firing rate (spks/s)');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
axis square;

subplot(2,3,4); hold on;
bar([1 4 7],mean(sim2FRcat),0.3,'EdgeColor','none','FaceColor',[.6 .6 .6],'FaceAlpha',0.5);
errorbar([1 4 7], mean(sim2FRcat), std(sim2FRcat),'k','CapSize',0,'LineWidth',2,'LineStyle','none');
plot([2*ones(1,nSim2); 3*ones(1,nSim2)], [sim2FRcat(:,1) sim2FRcat(:,2)]','LineWidth',1,'Color',[0 0 0 0.2]);
plot([5*ones(1,nSim2); 6*ones(1,nSim2)], [sim2FRcat(:,2) sim2FRcat(:,3)]','LineWidth',1,'Color',[0 0 0 0.2]);
ylim([0 3]); yticks([0:0.5:3]);
xlim([0 8]); xticks([1,4,7]); xticklabels({'Naive','Trained','Adapt.'})
ylabel('Simplex PC Firing rate (spks/s)');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
axis square;

subplot(2,3,5); hold on;
bar([1 4 7],mean(sim2FR),0.3,'EdgeColor','none','FaceColor',[.6 .6 .6],'FaceAlpha',0.5);
errorbar([1 4 7], mean(sim2FR), std(sim2FR),'k','CapSize',0,'LineWidth',2,'LineStyle','none');
plot([2*ones(1,j); 3*ones(1,j)], [sim2FR(:,1) sim2FR(:,2)]','LineWidth',1,'Color',[0 0 0 0.2]);
plot([5*ones(1,j); 6*ones(1,j)], [sim2FR(:,2) sim2FR(:,3)]','LineWidth',1,'Color',[0 0 0 0.2]);
ylim([0 3]); yticks([0:0.5:3]);
xlim([0 8]); xticks([1,4,7]); xticklabels({'Naive','Trained','Adapt.'})
ylabel('Simplex FOV Firing rate (spks/s)');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
axis square;

subplot(2,3,3); hold on;
bar([1 4],mean(lobvCorr),0.25,'EdgeColor','none','FaceColor',[.6 0 .6],'FaceAlpha',0.5);
errorbar([1 4], mean(lobvCorr), std(lobvCorr),'k','CapSize',0,'LineWidth',2,'LineStyle','none');
plot([2*ones(1,i); 3*ones(1,i)], [lobvCorr(:,1) lobvCorr(:,2)]','LineWidth',1,'Color',[0 0 0 0.2]);
bar(7+[1 4],mean(sim2Corr),0.25,'EdgeColor','none','FaceColor',[0 .6 .6],'FaceAlpha',0.5);
errorbar(7+[1 4], mean(sim2Corr), std(sim2Corr),'k','CapSize',0,'LineWidth',2,'LineStyle','none');
plot(7+[2*ones(1,j); 3*ones(1,j)], [sim2Corr(:,1) sim2Corr(:,2)]','LineWidth',1,'Color',[0 0 0 0.2]);
ylim([0 1]); yticks([0:0.2:1]);
xlim([0 12]); xticks([1,4,8,11]); xticklabels({' Naive\newlineTrained','Learning\newlineTrained',' Naive\newlineTrained','Learning\newlineTrained'})
ylabel('Correlation matrix similarity');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
xtickangle(0)
title(['\color[rgb]{.5 0 .5} LobV','\color{black} and','\color[rgb]{0 .5 .5} Simplex','\color{black} fields of view'])

matchPCs_v144.lobvPCs = lobvPCs;
matchPCs_v144.sim2PCs = sim2PCs;

if savebool
    figFold = fullfile(paperFold,'Figures','Fig4','Sub4_summary');
    if ~exist(figFold,'dir'); mkdir(figFold); end
    figFold2 = fullfile(dataFold,'paperFigs','Fig4','Sub4_summary');
    if ~exist(figFold2,'dir'); mkdir(figFold2); end
    
    save(fullfile(figFold,'matchPCs_v144.mat'),'lobvPCs','sim2PCs');
    save(fullfile(figFold2,'matchPCs_v144.mat'),'lobvPCs','sim2PCs');
    
    figname = 'spontCorrFR';
    savefig(figCorrFR, fullfile(figFold,[figname,'.fig']));
    saveas(figCorrFR,fullfile(figFold,[figname,'.png']));
    print(figCorrFR,fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    savefig(figCorrFR, fullfile(figFold2,[figname,'.fig']));
    saveas(figCorrFR,fullfile(figFold2,[figname,'.png']));
    print(figCorrFR,fullfile(figFold2,[figname,'.eps']), '-depsc', '-painters');
    
    
    figname = 'lobV_behav_heatmap';
    savefig(lobvfigs(1), fullfile(figFold,[figname,'.fig']));
    saveas(lobvfigs(1),fullfile(figFold,[figname,'.png']));
    print(lobvfigs(1),fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    savefig(lobvfigs(1), fullfile(figFold2,[figname,'.fig']));
    saveas(lobvfigs(1),fullfile(figFold2,[figname,'.png']));
    print(lobvfigs(1),fullfile(figFold2,[figname,'.eps']), '-depsc', '-painters');
    
    figname = 'lobV_reward_heatmap';
    savefig(lobvfigs(2), fullfile(figFold,[figname,'.fig']));
    saveas(lobvfigs(2),fullfile(figFold,[figname,'.png']));
    print(lobvfigs(2),fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    savefig(lobvfigs(2), fullfile(figFold2,[figname,'.fig']));
    saveas(lobvfigs(2),fullfile(figFold2,[figname,'.png']));
    print(lobvfigs(2),fullfile(figFold2,[figname,'.eps']), '-depsc', '-painters');
    
    figname = 'lobV_behav_summary';
    savefig(lobvfigs(3), fullfile(figFold,[figname,'.fig']));
    saveas(lobvfigs(3),fullfile(figFold,[figname,'.png']));
    print(lobvfigs(3),fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    savefig(lobvfigs(3), fullfile(figFold2,[figname,'.fig']));
    saveas(lobvfigs(3),fullfile(figFold2,[figname,'.png']));
    print(lobvfigs(3),fullfile(figFold2,[figname,'.eps']), '-depsc', '-painters');
    
    figname = 'lobV_reward_summary';
    savefig(lobvfigs(4), fullfile(figFold,[figname,'.fig']));
    saveas(lobvfigs(4),fullfile(figFold,[figname,'.png']));
    print(lobvfigs(4),fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    savefig(lobvfigs(4), fullfile(figFold2,[figname,'.fig']));
    saveas(lobvfigs(4),fullfile(figFold2,[figname,'.png']));
    print(lobvfigs(4),fullfile(figFold2,[figname,'.eps']), '-depsc', '-painters');
    
    figname = 'sim2_behav_heatmap';
    savefig(sim2figs(1), fullfile(figFold,[figname,'.fig']));
    saveas(sim2figs(1),fullfile(figFold,[figname,'.png']));
    print(sim2figs(1),fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    savefig(sim2figs(1), fullfile(figFold2,[figname,'.fig']));
    saveas(sim2figs(1),fullfile(figFold2,[figname,'.png']));
    print(sim2figs(1),fullfile(figFold2,[figname,'.eps']), '-depsc', '-painters');
    
    figname = 'sim2_reward_heatmap';
    savefig(sim2figs(2), fullfile(figFold,[figname,'.fig']));
    saveas(sim2figs(2),fullfile(figFold,[figname,'.png']));
    print(sim2figs(2),fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    savefig(sim2figs(2), fullfile(figFold2,[figname,'.fig']));
    saveas(sim2figs(2),fullfile(figFold2,[figname,'.png']));
    print(sim2figs(2),fullfile(figFold2,[figname,'.eps']), '-depsc', '-painters');
    
    figname = 'sim2_behav_summary';
    savefig(sim2figs(3), fullfile(figFold,[figname,'.fig']));
    saveas(sim2figs(3),fullfile(figFold,[figname,'.png']));
    print(sim2figs(3),fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    savefig(sim2figs(3), fullfile(figFold2,[figname,'.fig']));
    saveas(sim2figs(3),fullfile(figFold2,[figname,'.png']));
    print(sim2figs(3),fullfile(figFold2,[figname,'.eps']), '-depsc', '-painters');
    
    figname = 'sim2_reward_summary';
    savefig(sim2figs(4), fullfile(figFold,[figname,'.fig']));
    saveas(sim2figs(4),fullfile(figFold,[figname,'.png']));
    print(sim2figs(4),fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    savefig(sim2figs(4), fullfile(figFold2,[figname,'.fig']));
    saveas(sim2figs(4),fullfile(figFold2,[figname,'.png']));
    print(sim2figs(4),fullfile(figFold2,[figname,'.eps']), '-depsc', '-painters');
    
end

end
% load(...)
% v1rfs = RFs;
% load(...)
% v4fs = RFs;
% load(...);
% v5fs = RFs;
%
% load(...)
% v41rois = regi.rois.idcs;
% load(...)
% v45rois = regi.rois.idcs;
%
% [a,b] = ismember(v41rois(:,1),v45rois(:,1));
% [a2,b2] = ismember(v45rois(:,1),v41rois(:,1));
%
% v41rois_match0 = v41rois(a,:);
% v45rois_match0 = v45rois(a2,:);
%
% [~, v41order] = sort(v41rois_match0(:,1),'ascend');
% [~, v45order] = sort(v45rois_match0(:,1),'ascend');
% v41rois_match = v41rois_match0(v41order,:);
% v45rois_match = v45rois_match0(v45order,:);
%
% v145rois = [v41rois_match(:,2) v45rois_match];