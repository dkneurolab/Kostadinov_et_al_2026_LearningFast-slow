function roiStatsAll = fig2_Supp_roiCorrs_rev1(matchData, seshv4, nEx, paperFold,savebool)
% Wrapper function to plot sparse rois for each fov, correlate shapes to
% each other, compare to nearest neighbour

figFold = fullfile(paperFold,'Figures','Fig2','SuppFig3_roiCorrs_rev1');
if ~exist(figFold,'dir')
    mkdir(figFold);
elseif exist(fullfile(figFold,'roiStatsAll.mat'),'file')
    roiStatsAll = load(fullfile(figFold,'roiStatsAll.mat')); 
    
    roiStatsAll = roiStatsAll.roiStatsAll;
end


%% Loop through, process, and maybe save if data does not already exist

if ~exist('roiStatsAll','var')
tstart = tic;
roiStatsAll = struct;
for imatch = 1:numel(matchData)
    tstartlocal = tic;
    [roiStats, f1, f2] = fig2_Supp_roiCorrs_singleFOV_rev1(matchData(imatch), nEx);
    if imatch == 1
        roiStatsAll = roiStats;
    else
        roiStatsAll(imatch) = roiStats;
    end

    if savebool
        date1 = strrep(matchData(imatch).date1, '-', '');
        date2 = strrep(matchData(imatch).date2, '-', '');                
        figname = sprintf('%s_%s_%sv%i_%sv%i_fovimgs',matchData(imatch).name, matchData(imatch).fov, date1, matchData(imatch).v1, date2, matchData(imatch).v2);
        savefig(f1, fullfile(figFold,[figname,'.fig']));
        saveas(f1,fullfile(figFold,[figname,'.png']));
        print(f1,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');
        figname = sprintf('%s_%s_%sv%i_%sv%i_roiimgs',matchData(imatch).name,matchData(imatch).fov, date1, matchData(imatch).v1, date2, matchData(imatch).v2);
        savefig(f2, fullfile(figFold,[figname,'.fig']));
        saveas(f2,fullfile(figFold,[figname,'.png']));
        print(f2,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');
    end

fprintf('Finished %i of %i datasets.\n',imatch, numel(matchData))
fprintf('Dataset took: %.0f seconds.\n',toc(tstartlocal));
fprintf('Total time: %.0f seconds.\n',toc(tstart));
end
end

%% Match up v1-v4n-v5 FOVs and plot overlapping ROIs

for iv4 = 1:numel(seshv4)
    iv14 = [];
    iv45 = [];
    for jmatch = 1:numel(matchData)
        if strcmpi(seshv4(iv4).name, matchData(jmatch).name) && strcmpi(seshv4(iv4).date, matchData(jmatch).date2) && ...
                strcmpi(seshv4(iv4).fov, matchData(jmatch).fov) && matchData(jmatch).v1 == 1
            iv14 = jmatch;
        elseif strcmpi(seshv4(iv4).name, matchData(jmatch).name) && strcmpi(seshv4(iv4).date, matchData(jmatch).date1) && ...
                strcmpi(seshv4(iv4).fov, matchData(jmatch).fov) && matchData(jmatch).v2 == 5
            iv45 = jmatch;
        end        
    end
    if ~isempty(iv14) && ~isempty(iv45)
        [f3, f4] = fig2_Supp_roiCorrs_v145FOV_rev1(matchData(iv14), matchData(iv45),nEx);
        if savebool
            figname = sprintf('%s_%s_v145_mImgOverlay',matchData(iv14).name, matchData(iv14).fov);
            savefig(f3, fullfile(figFold,[figname,'.fig']));
            saveas(f3,fullfile(figFold,[figname,'.png']));
            print(f3,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');
            figname = sprintf('%s_%s_v145_roiOverlay',matchData(iv14).name, matchData(iv14).fov);
            savefig(f4, fullfile(figFold,[figname,'.fig']));
            saveas(f4,fullfile(figFold,[figname,'.png']));
            print(f4,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');
        end

    end
end

%% Calculate summary stats

[roiStatsAll, figSummary] = fig2_Supp_roiCorrs_summary_rev1(roiStatsAll);


%% save maybe
if savebool
    figFold = fullfile(paperFold,'Figures','Fig2','SuppFig3_roiCorrs_rev1');
    if ~exist(figFold,'dir'); mkdir(figFold); end
    figname = 'roiCrrs_summary';
    savefig(figSummary, fullfile(figFold,[figname,'.fig']));
    saveas(figSummary,fullfile(figFold,[figname,'.png']));
    print(figSummary,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');
    save(fullfile(figFold,'roiStatsAll.mat'),'roiStatsAll');
end

end