function figs = fig1_roiMatch_plot(matchData)
%% Split up naming info based on which image moved
if matchData.tForm.tform_direction(1) == 2
    titleRef = sprintf('%s, %s on %s: Version %i', upper(matchData.name), matchData.fov, matchData.date2, matchData.v2);
    titleWarp = sprintf('%s, %s on %s: Version %i', upper(matchData.name), matchData.fov, matchData.date1, matchData.v1);
    iRoiRef = matchData.iROIs2;
    iRoiWarp = matchData.iROIs1;
    figlabel1 = titleWarp;
    figlabel2 = titleRef;
else
    titleRef = sprintf('%s, %s on %s: Version %i', upper(matchData.name), matchData.fov, matchData.date1, matchData.v1);
    titleWarp = sprintf('%s, %s on %s: Version %i', upper(matchData.name), matchData.fov, matchData.date2, matchData.v2);
    iRoiRef = matchData.iROIs1;
    iRoiWarp = matchData.iROIs2;
    figlabel1 = titleRef;
    figlabel2 = titleWarp;
end
%% Load in settings and assemble cell arrays
roiColors = rand(1e5,1);
tform = matchData.tForm.tform;

% Static image first
imDataRef = matchData.imData(matchData.tForm.tform_direction(1));
if ~isfield(imDataRef.ops,'meanImg')
    imDataRef.ops.meanImg = imDataRef.ops.mimg;
    dims = size(imDataRef.ops.meanImg);
    xrange  = imDataRef.ops.xrange;
    yrange  = imDataRef.ops.yrange;
else
    dims = size(imDataRef.ops.meanImg);
    xrange  = imDataRef.ops.xrange;
    yrange  = imDataRef.ops.yrange;
end

nRois = size(imDataRef.final,2);
xpix = cell(nRois,1); ypix = xpix; lam = xpix;
for j = 1:size(imDataRef.final,2)
        ypix{j,1} = imDataRef.final(j).ypix;
        xpix{j,1} = imDataRef.final(j).xpix;        
        lam{j,1} = imDataRef.final(j).lam;
end
[re_xpix,re_ypix] = s2p_reindex_rois(xpix,ypix,xrange,yrange);
roiImgRef = s2p_rois2image(re_xpix,re_ypix,lam,dims,roiColors(1:nRois));
roiImgRef_matched = s2p_rois2image(re_xpix(iRoiRef),re_ypix(iRoiRef),lam(iRoiRef),dims,roiColors(iRoiRef));
roiImgRef_gray = s2p_rois2image(re_xpix(iRoiRef),re_ypix(iRoiRef),lam(iRoiRef),dims,...
    roiColors(iRoiRef),zeros(size(iRoiRef)));
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
dims = size(imDataWarp.ops.meanImg);
xrange  = imDataWarp.ops.xrange;
yrange  = imDataWarp.ops.yrange;
nRois = size(imDataWarp.final,2);
xpix = cell(nRois,1); ypix = xpix; lam = xpix;
for j = 1:size(imDataWarp.final,2)
        ypix{j,1} = imDataWarp.final(j).ypix;
        xpix{j,1} = imDataWarp.final(j).xpix;        
        lam{j,1} = imDataWarp.final(j).lam;
end
[re_xpix,re_ypix] = s2p_reindex_rois(xpix,ypix,xrange,yrange);
roiImgWarp = s2p_rois2image(re_xpix,re_ypix,lam,dims,roiColors(1:nRois));
roiImgWarp_matched = s2p_rois2image(re_xpix(iRoiWarp),re_ypix(iRoiWarp),lam(iRoiWarp),dims,roiColors(iRoiRef));
roiImgWarp_gray = s2p_rois2image(re_xpix(iRoiWarp),re_ypix(iRoiWarp),lam(iRoiWarp),dims,...
    roiColors(iRoiRef),zeros(size(iRoiRef)));
meanImgWarp = imDataWarp.ops.meanImg;
mu = median(meanImgWarp(:));
sd1 = mean(abs(meanImgWarp(meanImgWarp<mu) - mu));
sd2 = mean(abs(meanImgWarp(meanImgWarp>mu) - mu));
meanImgWarp_scaled = meanImgWarp;
meanImgWarp_scaled(meanImgWarp_scaled(:) < (mu-3*sd1)) = (mu-3*sd1);
meanImgWarp_scaled(meanImgWarp_scaled(:) > (mu+3*sd2)) = (mu+3*sd2);
meanImgWarp_scaled = meanImgWarp_scaled-min(meanImgWarp_scaled(:));
meanImgWarp_scaled = meanImgWarp_scaled/max(meanImgWarp_scaled(:));

meanImgWarp_tform = imwarp(meanImgWarp_scaled,tform,'OutputView',imref2d(dims));
roiImgWarp_matched_tform = imwarp(roiImgWarp_matched,tform,'OutputView',imref2d(dims));
roiImgWarp_gray_tform = imwarp(roiImgWarp_gray,tform,'OutputView',imref2d(dims));

clear dims j lam mu nRois re_xpix re_ypix sd1 sd2 tform xpix xrange ypix yrange

%% Fig 1 - Original FOV images

figs(1) = figure;
% Mean images
subplot(2,3,1);
imagesc(meanImgRef_scaled,[0 1]);
colormap(gray);
axis image; axis off; box off
title(titleRef);
subplot(2,3,4);
imagesc(meanImgWarp_scaled,[0 1]);
axis image; axis off; box off
title(titleWarp);

% All ROIs
subplot(2,3,2);
imagesc(roiImgRef);
axis image; axis off; box off;
subplot(2,3,5);
imagesc(roiImgWarp);
axis image; axis off; box off;

% Matched ROIs
subplot(2,3,3);
imagesc(roiImgRef_matched);
axis image; axis off; box off;
subplot(2,3,6);
imagesc(roiImgWarp_matched);
axis image; axis off; box off;

%% Fig 2 - Warped FOV images

figs(2) = figure;
% Mean images
subplot(2,3,1);
imagesc(meanImgRef_scaled, [0 1]);
colormap(gray);
axis image; axis off; box off
title(titleRef);
subplot(2,3,4);
imagesc(meanImgWarp_tform,[0 1]);
colormap(gray);
axis image; axis off; box off
title(titleWarp);

% Matched ROIs
subplot(2,3,2);
imagesc(roiImgRef_matched);
axis image; axis off; box off;
subplot(2,3,5);
imagesc(roiImgWarp_matched_tform);
axis image; axis off; box off;

% meanImg overlay
subplot(2,3,3);
meanImgMerge = cat(3,meanImgWarp_tform, meanImgRef_scaled, meanImgRef_scaled);
imagesc(meanImgMerge);
axis image; axis off; box off;
title('Mean image overlay');

%roiImg overlay
subplot(2,3,6);
roiImgMerge = cat(3,roiImgWarp_gray_tform(:,:,1), roiImgRef_gray(:,:,1), roiImgRef_gray(:,:,1));
imagesc(roiImgMerge);
axis image; axis off; box off;
title('ROI image overlay');

clear imDataRef imDataWarp meanImgMerge meanImgRef meanImgRef_scaled meanImgWarp meanImgWarp_scaled ...
    meanImgWarp_tform roiImgMerge roiImgRef roiImgRef_gray roiImgRef_matched roiImgWarp roiImgWarp_gray ...
    roiImgWarp_gray_tform roiImgWarp_matched roiImgWarp_matched_tform iRoiRef iRoiWarp

%% Fig 3 - On, mov, rew comparison
on = fig1_roiMatch_spkStruct(matchData,'on');
mov = fig1_roiMatch_spkStruct(matchData,'mov');
rewT = fig1_roiMatch_spkStruct(matchData,'rewT');
nRois = numel(matchData.iROIs1);
t2p = ((0.5:1:120)/30)-2;
tvr = ((0.5:1:400)/100)-2;
yscale = [0 10];

figs(3) = figure;
subplot(2,6,1);
imagesc(on(1).spkMap(end/2-29:end/2+30,:)',yscale);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([0.5:15:60.5]); xticklabels({});
yticks(0.5:20:nRois); yticklabels([1,20:20:nRois]);
box off;
ylabel(figlabel1);
title('Trial on');
subplot(2,6,2);
imagesc(on(2).spkMap(end/2-29:end/2+30,:)',yscale);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([0.5:15:60.5]); xticklabels({});
yticks(0.5:20:nRois); yticklabels({});
ylabel(figlabel2);
box off;
colormap viridis;

subplot(6,6,19);
h = boundedline(t2p,on(1).spkMu, on(1).spkSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 0 12]);
xticks([-1:0.5:1]); xticklabels({});
yticks([0:3:12]);
ylabel('Pop. response (spk/s)');

subplot(6,6,20);
h = boundedline(t2p,on(2).spkMu, on(2).spkSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 0 12]);
xticks([-1:0.5:1]); xticklabels({});
yticks([0:3:12]); yticklabels({});

subplot(6,6,25);
h = boundedline(tvr,on(1).lickMu,on(1).lickSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 0 1]);
xticks([-1:0.5:1]); xticklabels({});
yticks([0:.25:1]);
ylabel('Licking');

subplot(6,6,26);
h = boundedline(tvr,on(2).lickMu,on(2).lickSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 0 1]);
xticks([-1:0.5:1]); xticklabels({});
yticks([0:.25:1]); yticklabels({});

subplot(6,6,31);
h = boundedline(tvr,on(1).movMu,on(1).movSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 -5 20]);
xticks([-1:0.5:1]);
xlabel('Time (s)');
yticks([-5:5:20]);
ylabel('Movement (cm/s)');

subplot(6,6,32);
h = boundedline(tvr,on(2).movMu,on(2).movSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 -5 20]);
xticks([-1:0.5:1]); xticklabels({});
yticks([-5:5:20]); yticklabels({});


subplot(2,6,3);
imagesc(mov(1).spkMap(end/2-29:end/2+30,:)',yscale);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([0.5:15:60.5]); xticklabels({});
yticks(0.5:20:nRois); yticklabels([1,20:20:nRois]);
box off;
title('Movement on');
subplot(2,6,4);
imagesc(mov(2).spkMap(end/2-29:end/2+30,:)',yscale);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([0.5:15:60.5]); xticklabels({});
yticks(0.5:20:nRois); yticklabels({});
box off;
colormap viridis;

subplot(6,6,21);
h = boundedline(t2p,mov(1).spkMu,mov(1).spkSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 0 12]);
xticks([-1:0.5:1]); xticklabels({});
yticks([0:3:12]);

subplot(6,6,22);
h = boundedline(t2p,mov(2).spkMu,mov(2).spkSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 0 12]);
xticks([-1:0.5:1]); xticklabels({});
yticks([0:3:12]); yticklabels({});

subplot(6,6,27);
h = boundedline(tvr,mov(1).lickMu,mov(1).lickSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 0 1]);
xticks([-1:0.5:1]); xticklabels({});
yticks([0:.25:1]);

subplot(6,6,28);
h = boundedline(tvr,mov(2).lickMu,mov(2).lickSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 0 1]);
xticks([-1:0.5:1]); xticklabels({});
yticks([0:.25:1]); yticklabels({});

subplot(6,6,33);
h = boundedline(tvr,mov(1).movMu,mov(1).movSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 -5 20]);
xticks([-1:0.5:1]); xticklabels({});
yticks([-5:5:20]);

subplot(6,6,34);
h = boundedline(tvr,mov(2).movMu,mov(2).movSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 -5 20]);
xticks([-1:0.5:1]); xticklabels({});
yticks([-5:5:20]); yticklabels({});

subplot(2,6,5);
imagesc(rewT(1).spkMap(end/2-29:end/2+30,:)',yscale);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([0.5:15:60.5]); xticklabels({});
yticks(0.5:20:nRois); yticklabels([1,20:20:nRois]);
box off;
title('Trial reward');
subplot(2,6,6);
imagesc(rewT(2).spkMap(end/2-29:end/2+30,:)',yscale);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([0.5:15:60.5]); xticklabels({});
yticks(0.5:20:nRois); yticklabels({});
box off;
colormap viridis;

subplot(6,6,23);
h = boundedline(t2p,rewT(1).spkMu,rewT(1).spkSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 0 12]);
xticks([-1:0.5:1]); xticklabels({});
yticks([0:3:12]);
subplot(6,6,24);
h = boundedline(t2p,rewT(2).spkMu,rewT(2).spkSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 0 12]);
xticks([-1:0.5:1]); xticklabels({});
yticks([0:3:12]); yticklabels({});

subplot(6,6,29);
h = boundedline(tvr,rewT(1).lickMu,rewT(1).lickSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 0 1]);
xticks([-1:0.5:1]); xticklabels({});
yticks([0:.25:1]);

subplot(6,6,30);
h = boundedline(tvr,rewT(2).lickMu,rewT(2).lickSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 0 1]);
xticks([-1:0.5:1]); xticklabels({});
yticks([0:.25:1]); yticklabels({});

subplot(6,6,35);
h = boundedline(tvr,rewT(1).movMu,rewT(1).movSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 -5 20]);
xticks([-1:0.5:1]); xticklabels({});
yticks([-5:5:20]);

subplot(6,6,36);
h = boundedline(tvr,rewT(2).movMu,rewT(2).movSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 -5 20]);
xticks([-1:0.5:1]); xticklabels({});
yticks([-5:5:20]); yticklabels({});

%% Fig 4 - Mean responses - on, mov, rewT

figs(4) = figure;
subplot(3,3,1); hold on;
h = boundedline(t2p,on(1).spkMu, on(1).spkSD,'-b','alpha');
set(h,'linewidth',1.5);
h = boundedline(t2p,on(2).spkMu, on(2).spkSD,'-r','alpha');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 -3 12]);
xticks([-1:0.5:1]); xticklabels({});
yticks([-3:3:12]);
ylabel('Pop. response (spk/s)');
title('Trial onset');

subplot(3,3,2); hold on;
h = boundedline(t2p,mov(1).spkMu, mov(1).spkSD,'-b','alpha');
set(h,'linewidth',1.5);
h = boundedline(t2p,mov(2).spkMu, mov(2).spkSD,'-r','alpha');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 -3 12]);
xticks([-1:0.5:1]); xticklabels({});
yticks([-3:3:12]);
title('Movement onset');

subplot(3,3,3); hold on;
h = boundedline(t2p,rewT(1).spkMu, rewT(1).spkSD,'-b','alpha');
set(h,'linewidth',1.5);
h = boundedline(t2p,rewT(2).spkMu, rewT(2).spkSD,'-r','alpha');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 -3 12]);
xticks([-1:0.5:1]); xticklabels({});
yticks([-3:3:12]);
title('Trial reward');

subplot(3,3,4); hold on;
h = boundedline(tvr,on(1).lickMu,on(1).lickSD,'-b','alpha');
set(h,'linewidth',1.5);
h = boundedline(tvr,on(2).lickMu,on(2).lickSD,'-r','alpha');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 0 1]);
xticks([-1:0.5:1]); xticklabels({});
yticks([0:.25:1]);
ylabel('Licking');

subplot(3,3,5); hold on;
h = boundedline(tvr,mov(1).lickMu,mov(1).lickSD,'-b','alpha');
set(h,'linewidth',1.5);
h = boundedline(tvr,mov(2).lickMu,mov(2).lickSD,'-r','alpha');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 0 1]);
xticks([-1:0.5:1]); xticklabels({});
yticks([0:.25:1]);

subplot(3,3,6); hold on;
h = boundedline(tvr,rewT(1).lickMu,rewT(1).lickSD,'-b','alpha');
set(h,'linewidth',1.5);
h = boundedline(tvr,rewT(2).lickMu,rewT(2).lickSD,'-r','alpha');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 0 1]);
xticks([-1:0.5:1]); xticklabels({});
yticks([0:.25:1]);

subplot(3,3,7); hold on;
h = boundedline(tvr,on(1).movMu,on(1).movSD,'-b','alpha');
set(h,'linewidth',1.5);
h = boundedline(tvr,on(2).movMu,on(2).movSD,'-r','alpha');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 -5 20]);
xticks([-1:0.5:1]);
yticks([-5:5:20]);
ylabel('Movement (cm/s)');

subplot(3,3,8); hold on;
h = boundedline(tvr,mov(1).movMu,mov(1).movSD,'-b','alpha');
set(h,'linewidth',1.5);
h = boundedline(tvr,mov(2).movMu,mov(2).movSD,'-r','alpha');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 -5 20]);
xticks([-1:0.5:1]);
xlabel('Time (s)');
yticks([-5:5:20]);

subplot(3,3,9); hold on;
h = boundedline(tvr,rewT(1).movMu,rewT(1).movSD,'-b','alpha');
set(h,'linewidth',1.5);
h = boundedline(tvr,rewT(2).movMu,rewT(2).movSD,'-r','alpha');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 -5 20]);
xticks([-1:0.5:1]);
yticks([-5:5:20]);

suptitle({[figlabel1,' (blue)'];[figlabel2,' (red)']});

%% Fig 5 - Rand, trial, cued rew comparison

rewR = fig1_roiMatch_spkStruct(matchData,'rewR');
rewC = fig1_roiMatch_spkStruct(matchData,'rewC');
rewT = fig1_roiMatch_spkStruct(matchData,'rewT');

figs(5) = figure;
subplot(2,6,1);
imagesc(rewR(1).spkMap(end/2-29:end/2+30,:)',yscale);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([0.5:15:60.5]); xticklabels({});
yticks(0.5:20:nRois); yticklabels([1,20:20:nRois]);
box off;
ylabel(figlabel1);
title('Random reward');
subplot(2,6,2);
imagesc(rewR(2).spkMap(end/2-29:end/2+30,:)',yscale);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([0.5:15:60.5]); xticklabels({});
yticks(0.5:20:nRois); yticklabels({});
ylabel(figlabel2);
box off;
colormap viridis;

subplot(6,6,19);
h = boundedline(t2p,rewR(1).spkMu, rewR(1).spkSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 0 12]);
xticks([-1:0.5:1]); xticklabels({});
yticks([0:3:12]);
ylabel('Pop. response (spk/s)');

subplot(6,6,20);
h = boundedline(t2p,rewR(2).spkMu, rewR(2).spkSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 0 12]);
xticks([-1:0.5:1]); xticklabels({});
yticks([0:3:12]); yticklabels({});

subplot(6,6,25);
h = boundedline(tvr,rewR(1).lickMu,rewR(1).lickSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 0 1]);
xticks([-1:0.5:1]); xticklabels({});
yticks([0:.25:1]);
ylabel('Licking');

subplot(6,6,26);
h = boundedline(tvr,rewR(2).lickMu,rewR(2).lickSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 0 1]);
xticks([-1:0.5:1]); xticklabels({});
yticks([0:.25:1]); yticklabels({});

subplot(6,6,31);
h = boundedline(tvr,rewR(1).movMu,rewR(1).movSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 -5 20]);
xticks([-1:0.5:1]);
xlabel('Time (s)');
yticks([-5:5:20]);
ylabel('Movement (cm/s)');

subplot(6,6,32);
h = boundedline(tvr,rewR(2).movMu,rewR(2).movSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 -5 20]);
xticks([-1:0.5:1]); xticklabels({});
yticks([-5:5:20]); yticklabels({});


subplot(2,6,3);
imagesc(rewT(1).spkMap(end/2-29:end/2+30,:)',yscale);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([0.5:15:60.5]); xticklabels({});
yticks(0.5:20:nRois); yticklabels([1,20:20:nRois]);
box off;
title('Trial reward');
subplot(2,6,4);
imagesc(rewT(2).spkMap(end/2-29:end/2+30,:)',yscale);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([0.5:15:60.5]); xticklabels({});
yticks(0.5:20:nRois); yticklabels({});
box off;
colormap viridis;

subplot(6,6,21);
h = boundedline(t2p,rewT(1).spkMu,rewT(1).spkSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 0 12]);
xticks([-1:0.5:1]); xticklabels({});
yticks([0:3:12]);

subplot(6,6,22);
h = boundedline(t2p,rewT(2).spkMu,rewT(2).spkSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 0 12]);
xticks([-1:0.5:1]); xticklabels({});
yticks([0:3:12]); yticklabels({});

subplot(6,6,27);
h = boundedline(tvr,rewT(1).lickMu,rewT(1).lickSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 0 1]);
xticks([-1:0.5:1]); xticklabels({});
yticks([0:.25:1]);

subplot(6,6,28);
h = boundedline(tvr,rewT(2).lickMu,rewT(2).lickSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 0 1]);
xticks([-1:0.5:1]); xticklabels({});
yticks([0:.25:1]); yticklabels({});

subplot(6,6,33);
h = boundedline(tvr,rewT(1).movMu,rewT(1).movSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 -5 20]);
xticks([-1:0.5:1]); xticklabels({});
yticks([-5:5:20]);

subplot(6,6,34);
h = boundedline(tvr,rewT(2).movMu,rewT(2).movSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 -5 20]);
xticks([-1:0.5:1]); xticklabels({});
yticks([-5:5:20]); yticklabels({});

subplot(2,6,5);
imagesc(rewC(1).spkMap(end/2-29:end/2+30,:)',yscale);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([0.5:15:60.5]); xticklabels({});
yticks(0.5:20:nRois); yticklabels([1,20:20:nRois]);
box off;
title('Cued reward');
subplot(2,6,6);
imagesc(rewC(2).spkMap(end/2-29:end/2+30,:)',yscale);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([0.5:15:60.5]); xticklabels({});
yticks(0.5:20:nRois); yticklabels({});
box off;
colormap viridis;

subplot(6,6,23);
h = boundedline(t2p,rewC(1).spkMu,rewC(1).spkSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 0 12]);
xticks([-1:0.5:1]); xticklabels({});
yticks([0:3:12]);
subplot(6,6,24);
h = boundedline(t2p,rewC(2).spkMu,rewC(2).spkSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 0 12]);
xticks([-1:0.5:1]); xticklabels({});
yticks([0:3:12]); yticklabels({});

subplot(6,6,29);
h = boundedline(tvr,rewC(1).lickMu,rewC(1).lickSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 0 1]);
xticks([-1:0.5:1]); xticklabels({});
yticks([0:.25:1]);

subplot(6,6,30);
h = boundedline(tvr,rewC(2).lickMu,rewC(2).lickSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 0 1]);
xticks([-1:0.5:1]); xticklabels({});
yticks([0:.25:1]); yticklabels({});

subplot(6,6,35);
h = boundedline(tvr,rewC(1).movMu,rewC(1).movSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 -5 20]);
xticks([-1:0.5:1]); xticklabels({});
yticks([-5:5:20]);

subplot(6,6,36);
h = boundedline(tvr,rewC(2).movMu,rewC(2).movSD,'-k');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 -5 20]);
xticks([-1:0.5:1]); xticklabels({});
yticks([-5:5:20]); yticklabels({});

%% Fig 6 - - Mean responses - rewR, rewT, rewC

figs(6) = figure;
subplot(3,3,1); hold on;
h = boundedline(t2p,rewR(1).spkMu, rewR(1).spkSD,'-b','alpha');
set(h,'linewidth',1.5);
h = boundedline(t2p,rewR(2).spkMu, rewR(2).spkSD,'-r','alpha');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 -3 12]);
xticks([-1:0.5:1]); xticklabels({});
yticks([-3:3:12]);
ylabel('Pop. response (spk/s)');
title('Random reward');

subplot(3,3,2); hold on;
h = boundedline(t2p,rewT(1).spkMu, rewT(1).spkSD,'-b','alpha');
set(h,'linewidth',1.5);
h = boundedline(t2p,rewT(2).spkMu, rewT(2).spkSD,'-r','alpha');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 -3 12]);
xticks([-1:0.5:1]); xticklabels({});
yticks([-3:3:12]);
title('Trial reward');

subplot(3,3,3); hold on;
h = boundedline(t2p,rewC(1).spkMu, rewC(1).spkSD,'-b','alpha');
set(h,'linewidth',1.5);
h = boundedline(t2p,rewC(2).spkMu, rewC(2).spkSD,'-r','alpha');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 -3 12]);
xticks([-1:0.5:1]); xticklabels({});
yticks([-3:3:12]);
title('Cued reward');

subplot(3,3,4); hold on;
h = boundedline(tvr,rewR(1).lickMu,rewR(1).lickSD,'-b','alpha');
set(h,'linewidth',1.5);
h = boundedline(tvr,rewR(2).lickMu,rewR(2).lickSD,'-r','alpha');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 0 1]);
xticks([-1:0.5:1]); xticklabels({});
yticks([0:.25:1]);
ylabel('Licking');

subplot(3,3,5); hold on;
h = boundedline(tvr,rewT(1).lickMu,rewT(1).lickSD,'-b','alpha');
set(h,'linewidth',1.5);
h = boundedline(tvr,rewT(2).lickMu,rewT(2).lickSD,'-r','alpha');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 0 1]);
xticks([-1:0.5:1]); xticklabels({});
yticks([0:.25:1]);

subplot(3,3,6); hold on;
h = boundedline(tvr,rewC(1).lickMu,rewC(1).lickSD,'-b','alpha');
set(h,'linewidth',1.5);
h = boundedline(tvr,rewC(2).lickMu,rewC(2).lickSD,'-r','alpha');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 0 1]);
xticks([-1:0.5:1]); xticklabels({});
yticks([0:.25:1]);

subplot(3,3,7); hold on;
h = boundedline(tvr,rewR(1).movMu,rewR(1).movSD,'-b','alpha');
set(h,'linewidth',1.5);
h = boundedline(tvr,rewR(2).movMu,rewR(2).movSD,'-r','alpha');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 -5 20]);
xticks([-1:0.5:1]);
yticks([-5:5:20]);
ylabel('Movement (cm/s)');

subplot(3,3,8); hold on;
h = boundedline(tvr,rewT(1).movMu,rewT(1).movSD,'-b','alpha');
set(h,'linewidth',1.5);
h = boundedline(tvr,rewT(2).movMu,rewT(2).movSD,'-r','alpha');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 -5 20]);
xticks([-1:0.5:1]);
xlabel('Time (s)');
yticks([-5:5:20]);

subplot(3,3,9); hold on;
h = boundedline(tvr,rewC(1).movMu,rewC(1).movSD,'-b','alpha');
set(h,'linewidth',1.5);
h = boundedline(tvr,rewC(2).movMu,rewC(2).movSD,'-r','alpha');
set(h,'linewidth',1.5);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([-1 1 -5 20]);
xticks([-1:0.5:1]);
yticks([-5:5:20]);

suptitle({[figlabel1,' (blue)'];[figlabel2,' (red)']});

%%
figs(7) = figure;
subplot(1,2,1); hold on;
subplot(1,2,2); hold on;
if numel(matchData.iROIs1) >=5
    subplot(1,2,1);
    plot(matchData.imData(1).zF(:,end-4),'b');
    tSpks = find(matchData.imData(1).spks(:,end-4))';
    plot([tSpks;tSpks],7.5+[zeros(size(tSpks)); 2.5*ones(size(tSpks))],'-k','LineWidth',2)
    subplot(1,2,2);
    plot(matchData.imData(2).zF(:,end-4),'r');
    tSpks = find(matchData.imData(2).spks(:,end-4))';
    plot([tSpks;tSpks],7.5+[zeros(size(tSpks)); 2.5*ones(size(tSpks))],'-k','LineWidth',2)    
    yMax = 15;
end

if numel(matchData.iROIs1) >=25
    subplot(1,2,1);
    plot(matchData.imData(1).zF(:,end-24)+20,'b');
    tSpks = find(matchData.imData(1).spks(:,end-24))';
    plot([tSpks;tSpks],27.5+[zeros(size(tSpks)); 2.5*ones(size(tSpks))],'-k','LineWidth',2)
    subplot(1,2,2);
    plot(matchData.imData(2).zF(:,end-24)+20,'r');
    tSpks = find(matchData.imData(2).spks(:,end-24))';
    plot([tSpks;tSpks],27.5+[zeros(size(tSpks)); 2.5*ones(size(tSpks))],'-k','LineWidth',2)    
    yMax = 35;
end

if numel(matchData.iROIs1) >=45
    subplot(1,2,1);
    plot(matchData.imData(1).zF(:,end-44)+40,'b');
    tSpks = find(matchData.imData(1).spks(:,end-44))';
    plot([tSpks;tSpks],47.5+[zeros(size(tSpks)); 2.5*ones(size(tSpks))],'-k','LineWidth',2)
    subplot(1,2,2);
    plot(matchData.imData(2).zF(:,end-44)+40,'r');
    tSpks = find(matchData.imData(2).spks(:,end-44))';
    plot([tSpks;tSpks],47.5+[zeros(size(tSpks)); 2.5*ones(size(tSpks))],'-k','LineWidth',2)    
    yMax = 55;
end

if numel(matchData.iROIs1) >=65
    subplot(1,2,1);
    plot(matchData.imData(1).zF(:,end-64)+60,'b');
    tSpks = find(matchData.imData(1).spks(:,end-64))';
    plot([tSpks;tSpks],67.5+[zeros(size(tSpks)); 2.5*ones(size(tSpks))],'-k','LineWidth',2)
    subplot(1,2,2);
    plot(matchData.imData(2).zF(:,end-64)+60,'r');
    tSpks = find(matchData.imData(2).spks(:,end-64))';
    plot([tSpks;tSpks],67.5+[zeros(size(tSpks)); 2.5*ones(size(tSpks))],'-k','LineWidth',2)    
    yMax = 75;
end

if numel(matchData.iROIs1) >=85
    subplot(1,2,1);
    plot(matchData.imData(1).zF(:,end-84)+80,'b');
    tSpks = find(matchData.imData(1).spks(:,end-84))';
    plot([tSpks;tSpks],87.5+[zeros(size(tSpks)); 2.5*ones(size(tSpks))],'-k','LineWidth',2)
    subplot(1,2,2);
    plot(matchData.imData(2).zF(:,end-84)+80,'r');
    tSpks = find(matchData.imData(2).spks(:,end-84))';
    plot([tSpks;tSpks],87.5+[zeros(size(tSpks)); 2.5*ones(size(tSpks))],'-k','LineWidth',2)    
    yMax = 95;
end

if numel(matchData.iROIs1) >=105
    subplot(1,2,1);
    plot(matchData.imData(1).zF(:,end-104)+100,'b');
    tSpks = find(matchData.imData(1).spks(:,end-104))';
    plot([tSpks;tSpks],107.5+[zeros(size(tSpks)); 2.5*ones(size(tSpks))],'-k','LineWidth',2)
    subplot(1,2,2);
    plot(matchData.imData(2).zF(:,end-104)+100,'r');
    tSpks = find(matchData.imData(2).spks(:,end-104))';
    plot([tSpks;tSpks],107.5+[zeros(size(tSpks)); 2.5*ones(size(tSpks))],'-k','LineWidth',2)    
    yMax = 115;
end


subplot(1,2,1);
ylim([-5 yMax]); yticks([0:10:yMax]);
xlim([0 size(matchData.imData(1).zF,1)]); xticks([0:300:size(matchData.imData(1).zF,1)]);
% title({[figlabel1,' (blue)'];'Spikes (black)'});
subplot(1,2,2);
ylim([-5 yMax]); yticks([0:10:yMax]);
xlim([0 size(matchData.imData(1).zF,1)]); xticks([0:300:size(matchData.imData(1).zF,1)]);
% title({[figlabel2,' (red)'];'Spikes (black)'});

suptitle({[figlabel1,' (blue)'];[figlabel2,' (red)']});

%%
% figs(7) = figure;
% hold on;
% iYticks = [];
% if numel(matchData.iROIs1) >=5    
%     plot(matchData.imData(1).zF(:,end-4),'b');
%     plot(matchData.imData(2).zF(:,end-4)+5,'r');
%     iYticks = [iYticks; 0; 5];
% end
% if numel(matchData.iROIs1) >=15    
%     plot(matchData.imData(1).zF(:,end-14)+15,'b');
%     plot(matchData.imData(2).zF(:,end-14)+20,'r');
%     iYticks = [iYticks; 15; 20];
% end
% if numel(matchData.iROIs1) >=25    
%     plot(matchData.imData(1).zF(:,end-24)+30,'b');
%     plot(matchData.imData(2).zF(:,end-24)+35,'r');
%     iYticks = [iYticks; 30; 35];
% end
% if numel(matchData.iROIs1) >=35    
%     plot(matchData.imData(1).zF(:,end-34)+45,'b');
%     plot(matchData.imData(2).zF(:,end-34)+50,'r');
%     iYticks = [iYticks; 45; 50];
% end
% if numel(matchData.iROIs1) >=45    
%     plot(matchData.imData(1).zF(:,end-44)+60,'b');
%     plot(matchData.imData(2).zF(:,end-44)+65,'r');
%     iYticks = [iYticks; 60; 65];
% end
% if numel(matchData.iROIs1) >=55    
%     plot(matchData.imData(1).zF(:,end-54)+75,'b');
%     plot(matchData.imData(2).zF(:,end-54)+80,'r');
%     iYticks = [iYticks; 75; 80];
% end
% if numel(matchData.iROIs1) >=65    
%     plot(matchData.imData(1).zF(:,end-64)+90,'b');
%     plot(matchData.imData(2).zF(:,end-64)+95,'r');
%     iYticks = [iYticks; 90; 95];
% end
% ylim([-5 105]); yticks(iYticks);
% xticks(0:300:3000)

end