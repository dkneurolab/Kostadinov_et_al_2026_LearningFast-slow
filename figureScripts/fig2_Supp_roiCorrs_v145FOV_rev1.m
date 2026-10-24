function [f3, f4] = fig2_Supp_roiCorrs_v145FOV_rev1(matchData14,matchData45, nEx)
% Plot sparse rois for each fov

um2px = 670/512;

%% Split up naming info based on which image moved
if matchData14.tForm.tform_direction(1) == 2
    titleRef = sprintf('%s, %s on %s: Version %i', upper(matchData14.name), matchData14.fov, matchData14.date2, matchData14.v2);
    titleWarp14 = sprintf('%s, %s on %s: Version %i', upper(matchData14.name), matchData14.fov, matchData14.date1, matchData14.v1);
    iRoiRef14 = matchData14.iROIs2;
    iRoiWarp14 = matchData14.iROIs1;
else
    titleRef = sprintf('%s, %s on %s: Version %i', upper(matchData14.name), matchData14.fov, matchData14.date1, matchData14.v1);
    titleWarp14 = sprintf('%s, %s on %s: Version %i', upper(matchData14.name), matchData14.fov, matchData14.date2, matchData14.v2);
    iRoiRef14 = matchData14.iROIs1;
    iRoiWarp14 = matchData14.iROIs2;
end

if matchData45.tForm.tform_direction(1) == 2
    titleRef = sprintf('%s, %s on %s: Version %i', upper(matchData45.name), matchData45.fov, matchData45.date2, matchData45.v2);
    titleWarp45 = sprintf('%s, %s on %s: Version %i', upper(matchData45.name), matchData45.fov, matchData45.date1, matchData45.v1);
    iRoiRef45 = matchData45.iROIs2;
    iRoiWarp45 = matchData45.iROIs1;
else
    titleRef = sprintf('%s, %s on %s: Version %i', upper(matchData45.name), matchData45.fov, matchData45.date1, matchData45.v1);
    titleWarp45 = sprintf('%s, %s on %s: Version %i', upper(matchData45.name), matchData45.fov, matchData45.date2, matchData45.v2);
    iRoiRef45 = matchData45.iROIs1;
    iRoiWarp45 = matchData45.iROIs2;
end

[a,b] = ismember(iRoiRef14, iRoiRef45);
iRoiRef_sub0 = iRoiRef14(a);
iRoiWarp14_sub0 = iRoiWarp14(a);
% iRoiRef145_sub = iRoiRef45(b(a)); % Produces same index + order as iRoiRef14(a);
iRoiWarp45_sub0 = iRoiWarp45(b(a));

if numel(iRoiRef_sub0) > nEx
    skipper = floor(numel(iRoiRef_sub0) / nEx);
    iRoiRef_sub = iRoiRef_sub0(floor(skipper/2):skipper:end); iRoiRef_sub = iRoiRef_sub(1:nEx);
    iRoiWarp14_sub = iRoiWarp14_sub0(floor(skipper/2):skipper:end); iRoiWarp14_sub = iRoiWarp14_sub(1:nEx);
    iRoiWarp45_sub = iRoiWarp45_sub0(floor(skipper/2):skipper:end); iRoiWarp45_sub = iRoiWarp45_sub(1:nEx);
end

%% Load in settings and assemble cell arrays
rng(2349);
roiColors = rand(1e5,1);
tform14 = matchData14.tForm.tform;
tform45 = matchData45.tForm.tform;

% Static image first
imDataRef145 = matchData14.imData(matchData14.tForm.tform_direction(1));
if ~isfield(imDataRef145.ops,'meanImg')
    imDataRef145.ops.meanImg = imDataRef145.ops.mimg;
    dimsRef145 = size(imDataRef145.ops.meanImg);
    xrangeRef145  = imDataRef145.ops.xrange;
    yrangeRef145  = imDataRef145.ops.yrange;
else
    dimsRef145 = size(imDataRef145.ops.meanImg);
    xrangeRef145  = imDataRef145.ops.xrange;
    yrangeRef145  = imDataRef145.ops.yrange;
end

nRoisRef145 = size(imDataRef145.final,2);
xpixRef145 = cell(nRoisRef145,1); ypixRef145 = xpixRef145; lamRef145 = xpixRef145;
for j = 1:nRoisRef145
        ypixRef145{j,1} = imDataRef145.final(j).ypix;
        xpixRef145{j,1} = imDataRef145.final(j).xpix;        
        lamRef145{j,1} = imDataRef145.final(j).lam;
end
[re_xpixRef,re_ypixRef] = s2p_reindex_rois(xpixRef145,ypixRef145,xrangeRef145,yrangeRef145);
roiImgRef_matched = s2p_rois2image(re_xpixRef(iRoiRef_sub), re_ypixRef(iRoiRef_sub), lamRef145(iRoiRef_sub), dimsRef145, roiColors(1:numel(iRoiRef_sub)));
roiImgRef_gray = s2p_rois2image(re_xpixRef(iRoiRef_sub), re_ypixRef(iRoiRef_sub), lamRef145(iRoiRef_sub), dimsRef145, roiColors(1:numel(iRoiRef_sub)), zeros(size(iRoiRef_sub)));
roiImgRef0_matched = s2p_rois2image(re_xpixRef(iRoiRef_sub0), re_ypixRef(iRoiRef_sub0), lamRef145(iRoiRef_sub0), dimsRef145, roiColors(1:numel(iRoiRef_sub0)));
roiImgRef0_gray = s2p_rois2image(re_xpixRef(iRoiRef_sub0), re_ypixRef(iRoiRef_sub0), lamRef145(iRoiRef_sub0), dimsRef145, roiColors(1:numel(iRoiRef_sub0)), zeros(size(iRoiRef_sub0)));
meanImgRef = imDataRef145.ops.meanImg;
mu = median(meanImgRef(:));
sd1 = mean(abs(meanImgRef(meanImgRef<mu) - mu));
sd2 = mean(abs(meanImgRef(meanImgRef>mu) - mu));
meanImgRef_scaled = meanImgRef;
meanImgRef_scaled(meanImgRef_scaled(:) < (mu-3*sd1)) = (mu-3*sd1);
meanImgRef_scaled(meanImgRef_scaled(:) > (mu+3*sd2)) = (mu+3*sd2);
meanImgRef_scaled = meanImgRef_scaled-min(meanImgRef_scaled(:));
meanImgRef_scaled = meanImgRef_scaled/max(meanImgRef_scaled(:));

% Warped image next 14
imDataWarp14 = matchData14.imData(matchData14.tForm.tform_direction(2));
dimsWarp14 = size(imDataWarp14.ops.meanImg);
xrangeWarp14  = imDataWarp14.ops.xrange;
yrangeWarp14  = imDataWarp14.ops.yrange;
nRoisWarp = size(imDataWarp14.final,2);
xpixWarp14 = cell(nRoisWarp,1); ypixWarp14 = xpixWarp14; lamWarp14 = xpixWarp14;
for j = 1:size(imDataWarp14.final,2)
        ypixWarp14{j,1} = imDataWarp14.final(j).ypix;
        xpixWarp14{j,1} = imDataWarp14.final(j).xpix;        
        lamWarp14{j,1} = imDataWarp14.final(j).lam;
end
[re_xpixWarp14,re_ypixWarp14] = s2p_reindex_rois(xpixWarp14,ypixWarp14,xrangeWarp14,yrangeWarp14);
roiImgWarp14_matched = s2p_rois2image(re_xpixWarp14(iRoiWarp14_sub),re_ypixWarp14(iRoiWarp14_sub),lamWarp14(iRoiWarp14_sub),dimsWarp14,roiColors(1:numel(iRoiWarp14_sub)));
roiImgWarp14_gray = s2p_rois2image(re_xpixWarp14(iRoiWarp14_sub), re_ypixWarp14(iRoiWarp14_sub), lamWarp14(iRoiWarp14_sub), dimsWarp14, roiColors(1:numel(iRoiWarp14_sub)), zeros(size(iRoiWarp14_sub)));
roiImgWarp14_0_matched = s2p_rois2image(re_xpixWarp14(iRoiWarp14_sub0), re_ypixWarp14(iRoiWarp14_sub0), lamWarp14(iRoiWarp14_sub0), dimsWarp14, roiColors(1:numel(iRoiWarp14_sub0)));
roiImgWarp14_0_gray = s2p_rois2image(re_xpixWarp14(iRoiWarp14_sub0), re_ypixWarp14(iRoiWarp14_sub0), lamWarp14(iRoiWarp14_sub0), dimsWarp14, roiColors(1:numel(iRoiWarp14_sub0)), zeros(size(iRoiWarp14_sub0)));
meanImgWarp14 = imDataWarp14.ops.meanImg;
mu = median(meanImgWarp14(:));
sd1 = mean(abs(meanImgWarp14(meanImgWarp14<mu) - mu));
sd2 = mean(abs(meanImgWarp14(meanImgWarp14>mu) - mu));
meanImgWarp14_scaled = meanImgWarp14;
meanImgWarp14_scaled(meanImgWarp14_scaled(:) < (mu-3*sd1)) = (mu-3*sd1);
meanImgWarp14_scaled(meanImgWarp14_scaled(:) > (mu+3*sd2)) = (mu+3*sd2);
meanImgWarp14_scaled = meanImgWarp14_scaled-min(meanImgWarp14_scaled(:));
meanImgWarp14_scaled = meanImgWarp14_scaled/max(meanImgWarp14_scaled(:));

meanImgWarp14_tform = imwarp(meanImgWarp14_scaled,tform14,'OutputView',imref2d(dimsWarp14));
roiImgWarp14_matched_tform = imwarp(roiImgWarp14_matched,tform14,'OutputView',imref2d(dimsWarp14));
roiImgWarp14_gray_tform = imwarp(roiImgWarp14_gray,tform14,'OutputView',imref2d(dimsWarp14));
roiImgWarp14_0_matched_tform = imwarp(roiImgWarp14_0_matched,tform14,'OutputView',imref2d(dimsWarp14));
roiImgWarp14_0_gray_tform = imwarp(roiImgWarp14_0_gray,tform14,'OutputView',imref2d(dimsWarp14));

% Warped image next 45
imDataWarp45 = matchData45.imData(matchData45.tForm.tform_direction(2));
dimsWarp45 = size(imDataWarp45.ops.meanImg);
xrangeWarp45  = imDataWarp45.ops.xrange;
yrangeWarp45  = imDataWarp45.ops.yrange;
nRoisWarp = size(imDataWarp45.final,2);
xpixWarp45 = cell(nRoisWarp,1); ypixWarp45 = xpixWarp45; lamWarp45 = xpixWarp45;
for j = 1:size(imDataWarp45.final,2)
        ypixWarp45{j,1} = imDataWarp45.final(j).ypix;
        xpixWarp45{j,1} = imDataWarp45.final(j).xpix;        
        lamWarp45{j,1} = imDataWarp45.final(j).lam;
end
[re_xpixWarp45,re_ypixWarp45] = s2p_reindex_rois(xpixWarp45,ypixWarp45,xrangeWarp45,yrangeWarp45);
roiImgWarp45_matched = s2p_rois2image(re_xpixWarp45(iRoiWarp45_sub),re_ypixWarp45(iRoiWarp45_sub),lamWarp45(iRoiWarp45_sub),dimsWarp45,roiColors(1:numel(iRoiWarp45_sub)));
roiImgWarp45_gray = s2p_rois2image(re_xpixWarp45(iRoiWarp45_sub), re_ypixWarp45(iRoiWarp45_sub), lamWarp45(iRoiWarp45_sub), dimsWarp45, roiColors(1:numel(iRoiWarp45_sub)), zeros(size(iRoiWarp45_sub)));
roiImgWarp45_0_matched = s2p_rois2image(re_xpixWarp45(iRoiWarp45_sub0), re_ypixWarp45(iRoiWarp45_sub0), lamWarp45(iRoiWarp45_sub0), dimsWarp45, roiColors(1:numel(iRoiWarp45_sub0)));
roiImgWarp45_0_gray = s2p_rois2image(re_xpixWarp45(iRoiWarp45_sub0), re_ypixWarp45(iRoiWarp45_sub0), lamWarp45(iRoiWarp45_sub0), dimsWarp45, roiColors(1:numel(iRoiWarp45_sub0)), zeros(size(iRoiWarp45_sub0)));
meanImgWarp45 = imDataWarp45.ops.meanImg;
mu = median(meanImgWarp45(:));
sd1 = mean(abs(meanImgWarp45(meanImgWarp45<mu) - mu));
sd2 = mean(abs(meanImgWarp45(meanImgWarp45>mu) - mu));
meanImgWarp45_scaled = meanImgWarp45;
meanImgWarp45_scaled(meanImgWarp45_scaled(:) < (mu-3*sd1)) = (mu-3*sd1);
meanImgWarp45_scaled(meanImgWarp45_scaled(:) > (mu+3*sd2)) = (mu+3*sd2);
meanImgWarp45_scaled = meanImgWarp45_scaled-min(meanImgWarp45_scaled(:));
meanImgWarp45_scaled = meanImgWarp45_scaled/max(meanImgWarp45_scaled(:));

meanImgWarp45_tform = imwarp(meanImgWarp45_scaled,tform45,'OutputView',imref2d(dimsWarp45));
roiImgWarp45_matched_tform = imwarp(roiImgWarp45_matched,tform45,'OutputView',imref2d(dimsWarp45));
roiImgWarp45_gray_tform = imwarp(roiImgWarp45_gray,tform45,'OutputView',imref2d(dimsWarp45));
roiImgWarp45_0_matched_tform = imwarp(roiImgWarp45_0_matched,tform45,'OutputView',imref2d(dimsWarp45));
roiImgWarp45_0_gray_tform = imwarp(roiImgWarp45_0_gray,tform45,'OutputView',imref2d(dimsWarp45));

%%
f3 = figure;
sgtitle(sprintf('%s: %s matched v145, fov alignment', matchData14.name, matchData14.fov));
% sgtitle(sprintf('%s: %s matched v%i-v%i, fov alignment', matchData.name, matchData.fov, matchData.v1, matchData.v2));
% Mean images
subplot(2,4,1);
imagesc(meanImgRef_scaled, [0 1]);
colormap(gray);
axis image; axis off; box off
title(titleRef);
subplot(2,4,2);
imagesc(meanImgWarp14_scaled, [0 1]);
colormap(gray);
axis image; axis off; box off
title(titleWarp14);
subplot(2,4,3);
imagesc(meanImgWarp14_tform,[0 1]);
colormap(gray);
axis image; axis off; box off
title([titleWarp14,' warped']);
subplot(2,4,4);
meanImgMerge14 = cat(3,meanImgWarp14_tform, meanImgRef_scaled, meanImgRef_scaled);
imagesc(meanImgMerge14);
axis image; axis off; box off;
title('Mean image overlay - 14');
subplot(2,4,6);
imagesc(meanImgWarp45_scaled, [0 1]);
colormap(gray);
axis image; axis off; box off
title(titleWarp45);
subplot(2,4,7);
imagesc(meanImgWarp45_tform,[0 1]);
colormap(gray);
axis image; axis off; box off
title([titleWarp45,' warped']);
subplot(2,4,8);
meanImgMerge45 = cat(3,meanImgWarp45_tform, meanImgRef_scaled, meanImgRef_scaled);
imagesc(meanImgMerge45);
axis image; axis off; box off;
title('Mean image overlay - 45');
subplot(2,4,5);
meanImgMerge145 = cat(3,meanImgWarp14_tform, meanImgWarp45_tform, meanImgRef_scaled);
imagesc(meanImgMerge145);
axis image; axis off; box off;
title('Mean image overlay - 145');


f4 = figure;
sgtitle(sprintf('%s: %s matched v145, roi alignment', matchData14.name, matchData14.fov));
% Matched ROIs
subplot(2,4,1);
imagesc(roiImgWarp14_matched_tform);
axis image; axis off; box off;
title('Sparse ROIs - v1')
subplot(2,4,2);
imagesc(roiImgRef_matched);
axis image; axis off; box off;
title('Sparse ROIs - v4')
subplot(2,4,3);
imagesc(roiImgWarp45_matched_tform);
axis image; axis off; box off;
title('Sparse ROIs - v5')
subplot(2,4,4);
roiImgMerge = cat(3,roiImgWarp14_gray_tform(:,:,1), roiImgWarp45_gray_tform(:,:,1), roiImgRef_gray(:,:,1));
imagesc(roiImgMerge);
axis image; axis off; box off;
title('Sparse merge - v145');

% Matched ROIs
subplot(2,4,5);
imagesc(roiImgWarp14_0_matched_tform);
axis image; axis off; box off;
title('All ROIs - v1')
subplot(2,4,6);
imagesc(roiImgRef0_matched);
axis image; axis off; box off;
title('All ROIs - v4')
subplot(2,4,7);
imagesc(roiImgWarp45_0_matched_tform);
axis image; axis off; box off;
title('All ROIs - v5')
subplot(2,4,8);
roiImgMerge = cat(3,roiImgWarp14_0_gray_tform(:,:,1), roiImgWarp45_0_gray_tform(:,:,1), roiImgRef0_gray(:,:,1));
imagesc(roiImgMerge);
axis image; axis off; box off;
title('All merge - v145');

clear meanImgMerge45 meanImgRef meanImgRef_scaled meanImgWarp14 meanImgWarp14_scaled meanImgWarp14_tform ...
    roiImgMerge roiImgRef roiImgRef_gray roiImgRef_matched roiImgWarp roiImgWarp14_gray roiImgWarp14_gray_tform roiImgWarp14_matched roiImgWarp14_matched_tform ...
    mu sd1 sd2 j

