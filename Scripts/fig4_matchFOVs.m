function fig4_matchFOVs(sessionsv1,v1,sessionsv4_n,v2,sortFOV,dataFold,paperFold,savebool)

% Load in big matchdata structure that has all task-aligned activity for
% all matched FOV pairs
matchData = load(fullfile(dataFold,'paperFigs','Fig1','Sub3_behavFigs','matchDataBehav.mat')); 
matchData = matchData.matchData;
% sortFOV = 2; % This is the FOV that has the sortmap

%% Get task-aligned data for example heatmap
aa = 0;
matchV14 = struct;
for iSesh = 1:numel(sessionsv4_n)
    for jSesh = 1:numel(sessionsv1)
        if strcmpi(sessionsv4_n(iSesh).name,sessionsv1(jSesh).name) && ...
                strcmpi(sessionsv4_n(iSesh).fov,sessionsv1(jSesh).fov)
            
            % Get matched IDs
            v1Sesh = sessionsv1(jSesh); 
            v4Sesh = sessionsv4_n(iSesh);
            aa = aa + 1;      
            for iMatch = 1:numel(matchData)
                if strcmpi(matchData(iMatch).name, v1Sesh.name) && ...
                    ((strcmpi(matchData(iMatch).date1, v1Sesh.date) && strcmpi(matchData(iMatch).date2, v4Sesh.date)) || ...
                    (strcmpi(matchData(iMatch).date1, v4Sesh.date) && strcmpi(matchData(iMatch).date2, v1Sesh.date)))                                       
                    
                    matchV14(aa).name = v1Sesh.name;
                    if strcmpi(v1Sesh.fov,'lobv'); matchV14(aa).ilobv = true; else; matchV14(aa).ilobv = false; end
                    if strcmpi(v1Sesh.fov,'sim2'); matchV14(aa).isim2 = true; else; matchV14(aa).isim2 = false; end                    
                    matchV14(aa).('on') = fig1_roiMatch_spkStruct(matchData(iMatch),'on');
                    matchV14(aa).('mov') = fig1_roiMatch_spkStruct(matchData(iMatch),'mov');
                    matchV14(aa).('rewT') = fig1_roiMatch_spkStruct(matchData(iMatch),'rewT');
                    matchV14(aa).('rewC') = fig1_roiMatch_spkStruct(matchData(iMatch),'rewC');
                    matchV14(aa).('rewR') = fig1_roiMatch_spkStruct(matchData(iMatch),'rewR');
                end
            end
        end
    end
end

%%
matchV14lobv = matchV14(cat(1,matchV14.ilobv));

for ilobv = 1:numel(matchV14lobv)
    tempStruct = fig3_matchFOVs_PlotSort(matchV14lobv(ilobv).on,sortFOV,false);
    lobVon(ilobv) = tempStruct; %#ok<*AGROW>
    tempStruct = fig3_matchFOVs_PlotSort(matchV14lobv(ilobv).mov,sortFOV,true);
    lobVmov(ilobv) = tempStruct;
    tempStruct = fig3_matchFOVs_PlotSort(matchV14lobv(ilobv).rewT,sortFOV,false);
    lobVrewT(ilobv) = tempStruct;
    tempStruct = fig3_matchFOVs_PlotSort(matchV14lobv(ilobv).rewC,sortFOV,false);
    lobVrewC(ilobv) = tempStruct;
    tempStruct = fig3_matchFOVs_PlotSort(matchV14lobv(ilobv).rewR,sortFOV,false);
    lobVrewR(ilobv) = tempStruct;
    
end

figV_OnMovRew = fig3_matchFOVs_PlotFig(lobVon,lobVmov,lobVrewT,'Lobule V',{'Trial onset','Wheel movement','Trial Reward'},{'Naive','Expert'});
figV_RewRCT = fig3_matchFOVs_PlotFig(lobVrewR,lobVrewC,lobVrewT,'Lobule V',{'Random reward','Cued reward','Trial Reward'},{'Naive','Expert'});

%%
matchV14lsim2 = matchV14(cat(1,matchV14.isim2));

for isim2 = 1:numel(matchV14lsim2)
    
    tempStruct = fig3_matchFOVs_PlotSort(matchV14lsim2(isim2).on,sortFOV,false);
    sim2on(isim2) = tempStruct;
    tempStruct = fig3_matchFOVs_PlotSort(matchV14lsim2(isim2).mov,sortFOV,true);
    sim2mov(isim2) = tempStruct;
    tempStruct = fig3_matchFOVs_PlotSort(matchV14lsim2(isim2).rewT,sortFOV,false);
    sim2rewT(isim2) = tempStruct;
    tempStruct = fig3_matchFOVs_PlotSort(matchV14lsim2(isim2).rewC,sortFOV,false);
    sim2rewC(isim2) = tempStruct;
    tempStruct = fig3_matchFOVs_PlotSort(matchV14lsim2(isim2).rewR,sortFOV,false);
    sim2rewR(isim2) = tempStruct;
    
end

figSim_OnMovRew = fig3_matchFOVs_PlotFig(sim2on,sim2mov,sim2rewT,'Simplex',{'Trial onset','Wheel movement','Trial Reward'},{'Naive','Expert'});
figSim_RewRCT = fig3_matchFOVs_PlotFig(sim2rewR,sim2rewC,sim2rewT,'Simplex',{'Random reward','Cued reward','Trial Reward'},{'Naive','Expert'});


%% Save maybe
if savebool
    figFold = fullfile(paperFold,'Figures','Fig3',sprintf('Sub1_%s%s',v1,v2'));
    if ~exist(figFold,'dir'); mkdir(figFold); end
    save(fullfile(figFold,sprintf('match_%s_%s.mat',v1,v2)),'matchV14');
    
    figname = 'lobV_OnMovRew_heatmap';
    savefig(figV_OnMovRew, fullfile(figFold,[figname,'.fig']));
    saveas(figV_OnMovRew,fullfile(figFold,[figname,'.png']));
    print(figV_OnMovRew,fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    
    figname = 'lobV_RewRCT_heatmap';
    savefig(figV_RewRCT, fullfile(figFold,[figname,'.fig']));
    saveas(figV_RewRCT,fullfile(figFold,[figname,'.png']));
    print(figV_RewRCT,fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    
    figname = 'sim2_OnMovRew_heatmap';
    savefig(figSim_OnMovRew, fullfile(figFold,[figname,'.fig']));
    saveas(figSim_OnMovRew,fullfile(figFold,[figname,'.png']));
    print(figSim_OnMovRew,fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    
    figname = 'sim2_RewRCT_heatmap';
    savefig(figSim_RewRCT, fullfile(figFold,[figname,'.fig']));
    saveas(figSim_RewRCT,fullfile(figFold,[figname,'.png']));
    print(figSim_RewRCT,fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');

end


