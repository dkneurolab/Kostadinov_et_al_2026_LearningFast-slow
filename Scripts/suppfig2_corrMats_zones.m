function suppfig2_corrMats_zones(seshs,dataFold,paperFold,savebool)
imfs = 30;
% Mov period: (computed as mean activity over -300 to 0 ms from movement onset)
movRange = [-imfs*0.3:0];
% Rew period: (computed as mean activity over 0 to +200 ms from reward delivery)
rewRange = [1:imfs*0.2];
% max value of log-based colormap
logplotmax = 1;

%% Loop through FOVs and extract corr matrix data, compute corr matrix, and extract mov and rew RFs

% parentPath = fullfile(dataFold, 'v4FOVs');
parentPath = fullfile(dataFold, '01_Behav+imaging', 'Version_4');

fovStruct = struct;

for ifov = 1:numel(seshs)
    foldName = sprintf('%s/%s_%s',seshs(ifov).name,seshs(ifov).date,seshs(ifov).fov);
    % % Loading in all the stuff:
    % Load in fov params as datdk
    datdk = load(fullfile(parentPath,foldName,'fovparams.mat'));
    datdk = datdk.fovparams;
    % Load in sets and spont map
    try
        zonesets = load(fullfile(parentPath,foldName,'zonesets.mat'));
        spontData = zonesets.zonesets.Fspont;
        sets = zonesets.zonesets;
        clear zonesets
    catch
        setsonly = load(fullfile(parentPath,foldName,'setsonly.mat'));
        spontData = setsonly.sets.Fspont;
        sets = setsonly.sets;
        clear setsonly
    end    
    % Load data for RF map
    RFs = load(fullfile(parentPath,foldName,'RFs','RFs.mat'));
    
    % % Compute pixel data for FOV:
    roiStats = datdk.final;
    xpix = cell(numel(roiStats),1);
    ypix = xpix;
    lam = xpix;
    for i = 1:numel(roiStats)
        xpix{i,1} = roiStats(i).xpix;
        ypix{i,1} = roiStats(i).ypix;
        lam{i,1} = roiStats(i).lam;
    end    
    xRange = (datdk.ops.xrange(1):datdk.ops.xrange(2));
    yRange = (datdk.ops.yrange(1):datdk.ops.yrange(2));
    
    try % This may have different names depending on whether it was new or old file
        dims = size(datdk.ops.meanImg);
    catch
        dims = size(datdk.ops.mimg);
    end    
    
    % % Fill in data struct with fov features
    fovStruct(ifov).name = seshs(ifov).name;
    fovStruct(ifov).date = seshs(ifov).date;
    fovStruct(ifov).fov = seshs(ifov).fov;
    fovStruct(ifov).spontData = spontData;
    spontCorr = corrcoef(spontData);
    spontCorr(1:size(spontCorr,1)+1:end) = 0;
    fovStruct(ifov).corrMat = spontCorr;
    fovStruct(ifov).spkMov = RFs.RFs.mover_trial.SPKmean;
    fovStruct(ifov).spkRew = RFs.RFs.trialrew.SPKmean;

    % % Compute mov and rew vectors:
    % 1. calculate baseline (first second) and divide through by it
    spkMov_base = mean(fovStruct(ifov).spkMov(1:imfs,:));
    spkRew_base = mean(fovStruct(ifov).spkRew(1:imfs,:));

    spkMov_mat = fovStruct(ifov).spkMov./spkMov_base;
    spkRew_mat = fovStruct(ifov).spkRew./spkRew_base;
    % 2. Average over relevant period
    spkMov_sig0 = mean(spkMov_mat(end/2+movRange,:));
    spkRew_sig0 = mean(spkRew_mat(end/2+rewRange,:));

    % 3. Take the log
    fovStruct(ifov).spkMov_sig = log(spkMov_sig0);
    fovStruct(ifov).spkRew_sig = log(spkRew_sig0);
    spkMov_plot0 = log(spkMov_sig0)';
    spkMov_plot0(spkMov_plot0 > logplotmax) = logplotmax;
    spkMov_plot0(spkMov_plot0 < -logplotmax) = -logplotmax;
    spkMov_plot0(isnan(spkMov_plot0)) = 0;
    spkRew_plot0 = log(spkRew_sig0)';
    spkRew_plot0(spkRew_plot0 > logplotmax) = logplotmax;
    spkRew_plot0(spkRew_plot0 < -logplotmax) = -logplotmax;
    spkRew_plot0(isnan(spkRew_plot0)) = 0;

    % 4. Convert log values to colormap cool hsv color    
    spkMov_plot = round(100*(spkMov_plot0+logplotmax) / (2*logplotmax))+1;
    spkRew_plot = round(100*(spkRew_plot0+logplotmax) / (2*logplotmax))+1;
    % spkMov_plot = repmat((spkMov_plot0+logplotmax) / (2*logplotmax),1,3);
    % spkRew_plot = (spkRew_plot0+logplotmax) / (2*logplotmax);
    cmap101 = rgb2hsv(redblue(101)); % 100 levels of rgb 'cool'
    spkMov_plot_col = cmap101(spkMov_plot,:);
    spkRew_plot_col = cmap101(spkRew_plot,:);

    % % Test plotting
    [re_xpix,re_ypix] = s2p_reindex_rois(xpix,ypix,xRange,yRange);
    spkMov_roi_im = s2p_rois2image(re_xpix,re_ypix,lam,dims,spkMov_plot_col);        
    spkRew_roi_im = s2p_rois2image(re_xpix,re_ypix,lam,dims,spkRew_plot_col);     
    % figure; 
    % subplot(1,2,1);
    % imagesc(spkMov_roi_im); axis off; box off; axis square;
    % subplot(1,2,2);
    % imagesc(spkRew_roi_im); axis off; box off; axis square;

    fovStruct(ifov).spkMov_roi_im = spkMov_roi_im;
    fovStruct(ifov).spkRew_roi_im = spkRew_roi_im;

end

a = 4;

%% For each FOV - make 1x3 figure of corr matrix, move mapping and rew mapping

for ifov = 1:numel(fovStruct)
    corrMat = fovStruct(ifov).corrMat;
    foldName = sprintf('%s: %s, %s',fovStruct(ifov).name,fovStruct(ifov).date,fovStruct(ifov).fov);
    fPlot = figure;
    ax1 = subplot(1,3,1);
    imagesc(fovStruct(ifov).corrMat, [0, 0.4])
    axis square; colormap(ax1, 'viridis'); box off;
    xticks([0.5,50:50:size(corrMat,1)]); yticks([0.5,50:50:size(corrMat,1)]);
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    title(foldName)
    colorbar
    
    ax2 = subplot(1,3,2);
    imagesc(fovStruct(ifov).spkMov_roi_im, [-4, 4])
    axis square; colormap(ax2,'redblue'); box off;
    xticks([0.5,64:64:size(fovStruct(ifov).spkMov_roi_im,1)]); yticks([0.5,64:64:size(fovStruct(ifov).spkMov_roi_im,2)]);
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    colorbar
    ax3 = subplot(1,3,3);
    imagesc(fovStruct(ifov).spkRew_roi_im, [-4, 4])
    axis square; colormap(ax3, 'redblue'); box off;
    xticks([0.5,64:64:size(fovStruct(ifov).spkRew_roi_im,1)]); yticks([0.5,64:64:size(fovStruct(ifov).spkRew_roi_im,2)]);
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    colorbar

    % Save individual files
    if savebool
        figname = sprintf('%s_%s_%s_corr+MovRewHeat',fovStruct(ifov).name,fovStruct(ifov).date,fovStruct(ifov).fov);    
        savFold = fullfile(paperFold,'Figures/SuppFig2_new');
        if ~isfolder(savFold); mkdir(savFold); end
        savefig(fPlot, fullfile(savFold,[figname,'.fig']));
        saveas(fPlot,fullfile(savFold,[figname,'.png']));
        print(fPlot,fullfile(savFold,[figname,'.eps']), '-depsc', '-vector');    
    end

end

%% Save data structure at end
if savebool
    savFold = fullfile(paperFold,'Figures/SuppFig2_new');
    if ~isfolder(savFold); mkdir(savFold); end
    save(fullfile(savFold,'corr+MovRewHeat.mat'), 'fovStruct')
end

end