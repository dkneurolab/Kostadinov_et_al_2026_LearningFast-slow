function matchStatsSummary = fig2_matchStats(matchData, matchPCs_v145, dataFold, paperFold, vSave, savebool)

figFold = fullfile(paperFold,'Figures','Fig4','Sub4_summary');
if ~exist(figFold,'dir'); mkdir(figFold); end
figFold2 = fullfile(dataFold,'paperFigs','Fig4','Sub4_summary');
if ~exist(figFold2,'dir'); mkdir(figFold2); end

matchStatsMat = struct;
a = 0;
for i = 1:numel(matchData)
    if matchData(i).v1 == 1 & matchData(i).v2 == 4
        a = a + 1;
        matchStatsMat(a).name = matchData(i).name;
        matchStatsMat(a).iv14 = true;
        matchStatsMat(a).iv45 = false;
        if strcmpi(matchData(i).fov,'lobv')
            matchStatsMat(a).ilobv = true;
            matchStatsMat(a).isim2 = false;
        elseif strcmpi(matchData(i).fov,'sim2')
            matchStatsMat(a).ilobv = false;
            matchStatsMat(a).isim2 = true;
        end
        matchStatsMat(a).nRois = numel(matchData(i).iROIs1);
    end
    if matchData(i).v1 == 4 & matchData(i).v2 == 5
        a = a + 1;
        matchStatsMat(a).name = matchData(i).name;
        matchStatsMat(a).iv14 = false;
        matchStatsMat(a).iv45 = true;
        if strcmpi(matchData(i).fov,'lobv')
            matchStatsMat(a).ilobv = true;
            matchStatsMat(a).isim2 = false;
        elseif strcmpi(matchData(i).fov,'sim2')
            matchStatsMat(a).ilobv = false;
            matchStatsMat(a).isim2 = true;
        end
        matchStatsMat(a).nRois = numel(matchData(i).iROIs1);
    end
end

nRois = cat(1,matchStatsMat.nRois);

matchStats.nLobv_14rois = sum(nRois(cat(1,matchStatsMat.iv14) & cat(1,matchStatsMat.ilobv)));
matchStats.nLobv_14fovs = sum(cat(1,matchStatsMat.iv14) & cat(1,matchStatsMat.ilobv));
matchStats.nLobv_45rois = sum(nRois(cat(1,matchStatsMat.iv45) & cat(1,matchStatsMat.ilobv)));
matchStats.nLobv_45fovs = sum(cat(1,matchStatsMat.iv45) & cat(1,matchStatsMat.ilobv));
matchStats.nLobv_145rois = size(cat(1,matchPCs_v145.lobvPCs.v145rois),1);
matchStats.nLobv_145fovs = numel(matchPCs_v145.lobvPCs);
matchStats.nSim2_14rois = sum(nRois(cat(1,matchStatsMat.iv14) & cat(1,matchStatsMat.isim2)));
matchStats.nSim2_14fovs = sum(cat(1,matchStatsMat.iv14) & cat(1,matchStatsMat.isim2));
matchStats.nSim2_45rois = sum(nRois(cat(1,matchStatsMat.iv45) & cat(1,matchStatsMat.isim2)));
matchStats.nSim2_45fovs = sum(cat(1,matchStatsMat.iv45) & cat(1,matchStatsMat.isim2));
matchStats.nLobv_14rois = sum(nRois(cat(1,matchStatsMat.iv14) & cat(1,matchStatsMat.ilobv)));
matchStats.nSim2_145rois = size(cat(1,matchPCs_v145.sim2PCs.v145rois),1);
matchStats.nSim2_145fovs = numel(matchPCs_v145.sim2PCs);

matchStatsSummary.matchStatsMat = matchStatsMat;
matchStatsSummary.matchStats = matchStats;


if savebool
    save(fullfile(figFold,sprintf('matchStats_%s.mat',vSave)),'matchStatsMat','matchStats');
    save(fullfile(figFold2,sprintf('matchStats_%s.mat',vSave)),'matchStatsMat','matchStats');
    
    
end

end