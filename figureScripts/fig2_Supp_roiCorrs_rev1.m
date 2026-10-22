function roiStatsAll = fig2_Supp_roiCorrs_rev1(matchData, sessionsv4_n, nEx, paperFold,savebool)
% Wrapper function to plot sparse rois for each fov, correlate shapes to
% each other, compare to nearest neighbour

%% Loop through, process, and maybe save
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
        figFold = fullfile(paperFold,'Figures','Fig2','SuppFig3_roiCorrs_rev1');
        if ~exist(figFold,'dir'); mkdir(figFold); end
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

%% Match up v1-v4n-v5 FOVs and plot overlapping ROIs


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
    save(fullfile(figFold,'roiStatsAll.mat','roiStatsAll'));
end

end