function [roiStats, f1, f2] = fig2_Supp_roiCorrs_singleFOV_rev1(matchData, nEx)
% Plot sparse rois for each fov - correlate shapes to each other, compare
% to nearest neighbour

um2px = 670/512;

%% Split up naming info based on which image moved
if matchData.tForm.tform_direction(1) == 2
    titleRef = sprintf('%s, %s on %s: Version %i', upper(matchData.name), matchData.fov, matchData.date2, matchData.v2);
    titleWarp = sprintf('%s, %s on %s: Version %i', upper(matchData.name), matchData.fov, matchData.date1, matchData.v1);
    iRoiRef = matchData.iROIs2;
    iRoiWarp = matchData.iROIs1;
    if numel(iRoiRef) > nEx
        skipper = floor(numel(iRoiRef) / nEx);
        iRoiRef_sub = iRoiRef(floor(skipper/2):skipper:end); iRoiRef_sub = iRoiRef_sub(1:nEx);
        iRoiWarp_sub = iRoiWarp(floor(skipper/2):skipper:end); iRoiWarp_sub = iRoiWarp_sub(1:nEx);
    end
    figlabel1 = titleWarp;
    figlabel2 = titleRef;
else
    titleRef = sprintf('%s, %s on %s: Version %i', upper(matchData.name), matchData.fov, matchData.date1, matchData.v1);
    titleWarp = sprintf('%s, %s on %s: Version %i', upper(matchData.name), matchData.fov, matchData.date2, matchData.v2);
    iRoiRef = matchData.iROIs1;
    iRoiWarp = matchData.iROIs2;
    if numel(iRoiRef) > nEx
        skipper = floor(numel(iRoiRef) / nEx);
        iRoiRef_sub = iRoiRef(floor(skipper/2):skipper:end); iRoiRef_sub = iRoiRef_sub(1:nEx);
        iRoiWarp_sub = iRoiWarp(floor(skipper/2):skipper:end); iRoiWarp_sub = iRoiWarp_sub(1:nEx);
    end
    figlabel1 = titleRef;
    figlabel2 = titleWarp;
end

%% Load in settings and assemble cell arrays
rng(2349);
roiColors = rand(1e5,1);
tform = matchData.tForm.tform;

% Static image first
imDataRef = matchData.imData(matchData.tForm.tform_direction(1));
if ~isfield(imDataRef.ops,'meanImg')
    imDataRef.ops.meanImg = imDataRef.ops.mimg;
    dimsRef = size(imDataRef.ops.meanImg);
    xrangeRef  = imDataRef.ops.xrange;
    yrangeRef  = imDataRef.ops.yrange;
else
    dimsRef = size(imDataRef.ops.meanImg);
    xrangeRef  = imDataRef.ops.xrange;
    yrangeRef  = imDataRef.ops.yrange;
end

nRoisRef = size(imDataRef.final,2);
xpixRef = cell(nRoisRef,1); ypixRef = xpixRef; lamRef = xpixRef;
for j = 1:nRoisRef
        ypixRef{j,1} = imDataRef.final(j).ypix;
        xpixRef{j,1} = imDataRef.final(j).xpix;        
        lamRef{j,1} = imDataRef.final(j).lam;
end
[re_xpixRef,re_ypixRef] = s2p_reindex_rois(xpixRef,ypixRef,xrangeRef,yrangeRef);
roiImgRef_matched = s2p_rois2image(re_xpixRef(iRoiRef_sub), re_ypixRef(iRoiRef_sub), lamRef(iRoiRef_sub), dimsRef, roiColors(iRoiRef_sub));
roiImgRef_gray = s2p_rois2image(re_xpixRef(iRoiRef_sub), re_ypixRef(iRoiRef_sub), lamRef(iRoiRef_sub), dimsRef, roiColors(iRoiRef_sub), zeros(size(iRoiRef_sub)));
meanImgRef = imDataRef.ops.meanImg;
mu = median(meanImgRef(:));
sd1 = mean(abs(meanImgRef(meanImgRef<mu) - mu));
sd2 = mean(abs(meanImgRef(meanImgRef>mu) - mu));
meanImgRef_scaled = meanImgRef;
meanImgRef_scaled(meanImgRef_scaled(:) < (mu-3*sd1)) = (mu-3*sd1);
meanImgRef_scaled(meanImgRef_scaled(:) > (mu+3*sd2)) = (mu+3*sd2);
meanImgRef_scaled = meanImgRef_scaled-min(meanImgRef_scaled(:));
meanImgRef_scaled = meanImgRef_scaled/max(meanImgRef_scaled(:));

% Warped image next
imDataWarp = matchData.imData(matchData.tForm.tform_direction(2));
dimsWarp = size(imDataWarp.ops.meanImg);
xrangeWarp  = imDataWarp.ops.xrange;
yrangeWarp  = imDataWarp.ops.yrange;
nRoisWarp = size(imDataWarp.final,2);
xpixWarp = cell(nRoisWarp,1); ypixWarp = xpixWarp; lamWarp = xpixWarp;
for j = 1:size(imDataWarp.final,2)
        ypixWarp{j,1} = imDataWarp.final(j).ypix;
        xpixWarp{j,1} = imDataWarp.final(j).xpix;        
        lamWarp{j,1} = imDataWarp.final(j).lam;
        
end
[re_xpixWarp,re_ypixWarp] = s2p_reindex_rois(xpixWarp,ypixWarp,xrangeWarp,yrangeWarp);
roiImgWarp_matched = s2p_rois2image(re_xpixWarp(iRoiWarp_sub),re_ypixWarp(iRoiWarp_sub),lamWarp(iRoiWarp_sub),dimsWarp,roiColors(iRoiRef_sub));
roiImgWarp_gray = s2p_rois2image(re_xpixWarp(iRoiWarp_sub), re_ypixWarp(iRoiWarp_sub), lamWarp(iRoiWarp_sub), dimsWarp, roiColors(iRoiRef_sub), zeros(size(iRoiRef_sub)));
meanImgWarp = imDataWarp.ops.meanImg;
mu = median(meanImgWarp(:));
sd1 = mean(abs(meanImgWarp(meanImgWarp<mu) - mu));
sd2 = mean(abs(meanImgWarp(meanImgWarp>mu) - mu));
meanImgWarp_scaled = meanImgWarp;
meanImgWarp_scaled(meanImgWarp_scaled(:) < (mu-3*sd1)) = (mu-3*sd1);
meanImgWarp_scaled(meanImgWarp_scaled(:) > (mu+3*sd2)) = (mu+3*sd2);
meanImgWarp_scaled = meanImgWarp_scaled-min(meanImgWarp_scaled(:));
meanImgWarp_scaled = meanImgWarp_scaled/max(meanImgWarp_scaled(:));

meanImgWarp_tform = imwarp(meanImgWarp_scaled,tform,'OutputView',imref2d(dimsWarp));
roiImgWarp_matched_tform = imwarp(roiImgWarp_matched,tform,'OutputView',imref2d(dimsWarp));
roiImgWarp_gray_tform = imwarp(roiImgWarp_gray,tform,'OutputView',imref2d(dimsWarp));

%%
f1 = figure;
sgtitle(sprintf('%s: %s matched v%i-v%i, fov alignment', matchData.name, matchData.fov, matchData.v1, matchData.v2));
% Mean images
subplot(2,4,1);
imagesc(meanImgRef_scaled, [0 1]);
colormap(gray);
axis image; axis off; box off
title(titleRef);
subplot(2,4,2);
imagesc(meanImgWarp_scaled, [0 1]);
colormap(gray);
axis image; axis off; box off
title(titleWarp);
subplot(2,4,3);
imagesc(meanImgWarp_tform,[0 1]);
colormap(gray);
axis image; axis off; box off
title([titleWarp,' warped']);
subplot(2,4,4);
meanImgMerge = cat(3,meanImgWarp_tform, meanImgRef_scaled, meanImgRef_scaled);
imagesc(meanImgMerge);
axis image; axis off; box off;
title('Mean image overlay');

% Matched ROIs
subplot(2,4,5);
imagesc(roiImgRef_matched);
axis image; axis off; box off;
subplot(2,4,6);
imagesc(roiImgWarp_matched_tform);
axis image; axis off; box off;
%roiImg overlay
subplot(2,4,7);
roiImgMerge = cat(3,roiImgWarp_gray_tform(:,:,1), roiImgRef_gray(:,:,1), roiImgRef_gray(:,:,1));
imagesc(roiImgMerge);
axis image; axis off; box off;
title('ROI image overlay');

clear meanImgMerge meanImgRef meanImgRef_scaled meanImgWarp meanImgWarp_scaled meanImgWarp_tform ...
    roiImgMerge roiImgRef roiImgRef_gray roiImgRef_matched roiImgWarp roiImgWarp_gray roiImgWarp_gray_tform roiImgWarp_matched roiImgWarp_matched_tform ...
    mu sd1 sd2 j

%% Loop through individual ROIs, make 
f2 = figure;
sgtitle(sprintf('%s: %s matched v%i-v%i, example rois', matchData.name, matchData.fov, matchData.v1, matchData.v2));
irow = 0;
r_real = nan(numel(iRoiRef),1); r_next = r_real;
percFootprint_real = r_real; percFootprint_next = r_real;
ctrDist_real = r_real; ctrDist_next = r_real;

for iroi = 1:numel(iRoiRef)
    roi_ref = iRoiRef(iroi);
    roiImgRef = s2p_rois2image(re_xpixRef(roi_ref), re_ypixRef(roi_ref), lamRef(roi_ref), dimsRef, roiColors(roi_ref), zeros(size(roi_ref)));
    roi_warp = iRoiWarp(iroi);
    roi_next = iRoiWarp(iroi)+1;
    if roi_next > numel(re_xpixWarp); roi_next = roi_next - 2; end
    roiImgWarp0 = s2p_rois2image(re_xpixWarp(roi_warp), re_ypixWarp(roi_warp), lamWarp(roi_warp), dimsWarp, roiColors(roi_warp), zeros(roi_warp));
    roiImgWarp = imwarp(roiImgWarp0,tform,'OutputView',imref2d(dimsWarp));
    roiImgNext0 = s2p_rois2image(re_xpixWarp(roi_next), re_ypixWarp(roi_next), lamWarp(roi_next), dimsWarp, roiColors(roi_next), zeros(roi_next));
    roiImgNext = imwarp(roiImgNext0,tform,'OutputView',imref2d(dimsWarp));

    % Inter-centroid dist:
    xyRef0 = double([re_xpixRef{roi_ref} re_ypixRef{roi_ref}]);
    xyRef = transformPointsInverse(tform, xyRef0);
    xyWarp = double([re_xpixWarp{roi_warp} re_ypixWarp{roi_warp}]);
    xyNext = double([re_xpixWarp{roi_next} re_ypixWarp{roi_next}]);
    xDist_real = median(xyRef(:,1)) - median(xyWarp(:,1));
    xDist_next = median(xyRef(:,1)) - median(xyNext(:,1));
    yDist_real = median(xyRef(:,2)) - median(xyWarp(:,2));
    yDist_next = median(xyRef(:,2)) - median(xyNext(:,2));
    ctrDist_real(iroi) = sqrt(double(xDist_real^2 + yDist_real^2))*um2px;
    ctrDist_next(iroi) = sqrt(double(xDist_next^2 + yDist_next^2))*um2px;

    
    xminRef = find(mean(roiImgRef(:,:,1),2) > 0,1); xmaxRef = find(mean(roiImgRef(:,:,1),2) > 0,1,'last');
    yminRef = find(mean(roiImgRef(:,:,1),1) > 0,1); ymaxRef = find(mean(roiImgRef(:,:,1),1) > 0,1,'last');
    xminWarp = find(mean(roiImgWarp(:,:,1),2) > 0,1); xmaxWarp = find(mean(roiImgWarp(:,:,1),2) > 0,1,'last');
    yminWarp = find(mean(roiImgWarp(:,:,1),1) > 0,1); ymaxWarp = find(mean(roiImgWarp(:,:,1),1) > 0,1,'last');
    xminNext = find(mean(roiImgNext(:,:,1),2) > 0,1); xmaxNext = find(mean(roiImgNext(:,:,1),2) > 0,1,'last');
    yminNext = find(mean(roiImgNext(:,:,1),1) > 0,1); ymaxNext = find(mean(roiImgNext(:,:,1),1) > 0,1,'last');
    xminAll = min([xminRef, xminWarp, xminNext]); yminAll = min([yminRef, yminWarp, yminNext]);
    xmaxAll = max([xmaxRef, xmaxWarp, xmaxNext]); ymaxAll = max([ymaxRef, ymaxWarp, ymaxNext]);
    a = 1;
    while xmaxAll-xminAll < 20
        if a > 0
            xmaxAll = xmaxAll + 1;
        else
            xminAll = xminAll - 1;
        end
        a = -a;
    end
    if xminAll < 1
        xmaxAll = xMaxAll + (abs(xminAll))+1;
        xminAll = 1;
    end
    if xmaxAll > 512
        xminAll = xminAll - (xmaxAll-512);
        xmaxAll = 512;
    end
    a = 1;
    while ymaxAll-yminAll < 20
        if a > 0
            ymaxAll = ymaxAll + 1;
        else
            yminAll = yminAll - 1;
        end
        a = -a;
    end
    if yminAll < 1
        ymaxAll = ymaxAll + (abs(yminAll))+1;
        yminAll = 1;
    end
    if ymaxAll > 512
        yminAll = yminAll - (ymaxAll-512);
        ymaxAll = 512;
    end
    
    % Crop to common rectangle
    roiImgRef_roi = roiImgRef(xminAll:xmaxAll,yminAll:ymaxAll);
    roiImgWarp_roi = roiImgWarp(xminAll:xmaxAll,yminAll:ymaxAll);
    roiImgNext_roi = roiImgNext(xminAll:xmaxAll,yminAll:ymaxAll);
    % Shift Next image so that it's center is in the same place as the Warp
    xDist_warp_next = median(xyWarp(:,1)) - median(xyNext(:,1));
    yDist_warp_next = median(xyWarp(:,2)) - median(xyNext(:,2));
    roiImgNext_roi = circshift(roiImgNext_roi,round(xDist_warp_next),2);
    roiImgNext_roi = circshift(roiImgNext_roi,round(yDist_warp_next),1);
    
    % % Crop images and pad if they are too small
    % roiImgRef_roi = roiImgRef(xminRef:xmaxRef,yminRef:ymaxRef,1);
    % npad = [0,0];
    % if size(roiImgRef_roi,1) < 20; npad(1) = ceil((20-size(roiImgRef_roi,1))/2); end
    % if size(roiImgRef_roi,2) < 20; npad(2) = ceil((20-size(roiImgRef_roi,2))/2); end
    % roiImgRef_roi = padarray(roiImgRef_roi,npad);    
    % 
    % roiImgWarp_roi = roiImgWarp(xminWarp:xmaxWarp,yminWarp:ymaxWarp,1);
    % npad = [0,0];
    % if size(roiImgWarp_roi,1) < 20; npad(1) = ceil((20-size(roiImgWarp_roi,1))/2); end
    % if size(roiImgWarp_roi,2) < 20; npad(2) = ceil((20-size(roiImgWarp_roi,2))/2); end
    % roiImgWarp_roi = padarray(roiImgWarp_roi,npad);
    % 
    % roiImgNext_roi = roiImgNext(xminNext:xmaxNext,yminNext:ymaxNext,1);
    % npad = [0,0];
    % if size(roiImgNext_roi,1) < 20; npad(1) = ceil((20-size(roiImgNext_roi,1))/2); end
    % if size(roiImgNext_roi,2) < 20; npad(2) = ceil((20-size(roiImgNext_roi,2))/2); end
    % roiImgNext_roi = padarray(roiImgNext_roi,npad);

    % refit roi matching for individual rois
    [optimizer, metric] = imregconfig('monomodal');   % 'multimodal' if intensity scales differ
    optimizer.MaximumIterations = 500;
    optimizer.MaximumStepLength = 0.01;
    optimizer.MinimumStepLength = 1e-6;
    optimizer.RelaxationFactor  = 0.7;
    % Translation first, then rigid seeded with it (more stable)
    % tform_Warp0 = imregtform(roiImgWarp_roi, roiImgRef_roi, 'translation', optimizer, metric);
    tform_Warp  = imregtform(roiImgWarp_roi, roiImgRef_roi, 'rigid', optimizer, metric); %, 'InitialTransformation', tform_Warp0);
    % tform_Next0 = imregtform(roiImgNext_roi, roiImgRef_roi, 'translation', optimizer, metric);
    tform_Next  = imregtform(roiImgNext_roi, roiImgRef_roi, 'rigid', optimizer, metric); %, 'InitialTransformation', tform_Next0);

    Rfixed    = imref2d(size(roiImgRef_roi));
    roiImgWarp_roi_tform = imwarp(roiImgWarp_roi, tform_Warp, 'OutputView', Rfixed);
    roiImgNext_roi_tform = imwarp(roiImgNext_roi, tform_Next, 'OutputView', Rfixed);
    r_real(iroi) = corr(roiImgRef_roi(:), roiImgWarp_roi_tform(:));
    r_next(iroi) = corr(roiImgRef_roi(:), roiImgNext_roi_tform(:));

    % Footprint overlap
    roiImgRef_roi_bool = logical(roiImgRef_roi); % npixRef = sum(roiImgRef_roi_bool);
    roiImgWarp_roi_tform_bool = logical(roiImgWarp_roi_tform); % npixWarp = sum(roiImgWarp_roi_tform_bool);
    roiImgNext_roi_tform_bool = logical(roiImgNext_roi_tform); % npixNext = sum(roiImgNext_roi_tform_bool);
    percFootprint_real(iroi) = sum(roiImgRef_roi_bool(:)+roiImgWarp_roi_tform_bool(:) > 1) / min([sum(roiImgRef_roi_bool(:)),sum(roiImgWarp_roi_tform_bool(:))]);
    percFootprint_next(iroi) = sum(roiImgRef_roi_bool(:)+roiImgNext_roi_tform_bool(:) > 1) / min([sum(roiImgRef_roi_bool(:)),sum(roiImgNext_roi_tform_bool(:))]);
    
    % Plot real and next (comment out when doing at scale)
    if ismember(roi_ref, iRoiRef_sub)
        figure(f2);
        subplot(3,nEx,irow+1);
        imagesc(roiImgRef_roi); axis image;
        title(sprintf('iroi: %i, roi %i',iroi, roi_ref))
        subplot(3,nEx,nEx+irow+1);
        imagesc(roiImgWarp_roi_tform); axis image;
        title(sprintf('iroi: %i, roi %i',iroi, roi_warp))
        subplot(3,nEx,2*nEx+irow+1);
        imagesc(roiImgNext_roi_tform); axis image;    
        title(sprintf('iroi: %i, roi %i',iroi, roi_next))
        irow = irow+1;
        colormap viridis
    end

    if percFootprint_real(iroi) < 0.2
        figure;
        subplot(2,3,1);
        imagesc(roiImgRef_roi); axis image;
        title(sprintf('iroi: %i, roi %i',iroi, roi_ref))
        subplot(2,3,2);
        imagesc(roiImgWarp_roi_tform); axis image;
        title(sprintf('iroi: %i, roi %i',iroi, roi_warp))
        subplot(2,3,3);
        imagesc(roiImgNext_roi_tform); axis image;    
        title(sprintf('iroi: %i, roi %i',iroi, roi_next))
        subplot(2,3,4);
        imagesc(roiImgRef_roi); axis image;
        subplot(2,3,5);
        imagesc(roiImgWarp_roi); axis image;
        subplot(2,3,6);
        imagesc(roiImgNext_roi); axis image;
        colormap viridis        
    end

end

%% Assemble struct to pass back
roiStats.name = matchData.name;
roiStats.fov = matchData.fov;
roiStats.dateInj = matchData.dateInj;
roiStats.date1 = matchData.date1;
roiStats.date2 = matchData.date2;
roiStats.v1 = matchData.v1;
roiStats.v2 = matchData.v2;
roiStats.r_real = r_real;
roiStats.r_next = r_next;
roiStats.percFoot_real = percFootprint_real;
roiStats.percFoot_next = percFootprint_next;
roiStats.ctrDist_real = ctrDist_real;
roiStats.ctrDist_next = ctrDist_next;

fprintf('\n---Summary stats - %s: %s v%i-v%i---\n', matchData.name, matchData.fov, matchData.v1, matchData.v2);
fprintf('roi correlations (r):\n');
fprintf('  mean real: %.2f, mean next: %.2f\n',mean(r_real,'omitnan'), mean(r_next,'omitnan'));
fprintf('  %i of %i (%.2f%%) higher for real\n', sum(r_real > r_next), numel(r_real), 100*sum(r_real > r_next) / numel(r_real));
fprintf('footprint overlap (%%):\n');
fprintf('  mean real: %.2f, mean next: %.2f\n',mean(percFootprint_real,'omitnan'), mean(percFootprint_next,'omitnan'));
fprintf('  %i of %i (%.2f%%) higher for real\n', sum(percFootprint_real > percFootprint_next), numel(percFootprint_real), 100*sum(percFootprint_real > percFootprint_next) / numel(percFootprint_real));
fprintf('Centroid distances (um):\n');
fprintf('  mean real: %.2f, mean next: %.2f\n',mean(ctrDist_real,'omitnan'), mean(ctrDist_next,'omitnan'));
fprintf('  %i of %i (%.2f%%) closer for real\n', sum(ctrDist_real < ctrDist_next), numel(ctrDist_real), 100*sum(ctrDist_real < ctrDist_next) / numel(ctrDist_real));


end



