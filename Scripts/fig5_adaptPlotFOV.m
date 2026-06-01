function [figOut,dataOut] = fig5_adaptPlotFOV(dataIn,sets)
% Extract settings
trialGrp = sets.trialGrp;
imfs = sets.imfs;
vrfs = sets.vrfs;

trStarts = [dataIn.gainUp-trialGrp+1:round(trialGrp/2):dataIn.gainUp+120-trialGrp+1,dataIn.gainDown+1];
trEnds = [dataIn.gainUp:round(trialGrp/2):dataIn.gainUp+120,dataIn.gainDown+trialGrp];

%% Plot
spkMov = []; 
spkMovMap = [];
xMov = [];
vMov = [];
spkRew = [];
spkRewMap = [];
xRew = [];
vRew = [];
spkUnder = [];
spkUnderMap = [];
spkUnder_mov = [];
spkUnderMap_mov = [];
xUnder = [];
vUnder = [];
spkOver = [];
spkOverMap = [];
spkOver_mov = [];
spkOverMap_mov = [];
xOver = [];
vOver = [];

for iTr = 1:numel(trStarts)
    iTrGoodMov = find(ismember(dataIn.iTrialsMov,trStarts(iTr):trEnds(iTr)));
    iTrGoodRew = find(ismember(dataIn.iTrialsRew,trStarts(iTr):trEnds(iTr)));
    iTrGoodUnder = find(ismember(dataIn.iTrialsUnder,trStarts(iTr):trEnds(iTr)));
    iTrGoodOver = find(ismember(dataIn.iTrialsOver,trStarts(iTr):trEnds(iTr)));
    %spkMov
    spkMovMapLocal = mean(dataIn.spkMov(end/2-imfs+1:end/2+imfs,dataIn.sigBoolMov,iTrGoodMov),3,'omitnan');
    spkMovLocal = [mean(spkMovMapLocal,2,'omitnan') std(mean(dataIn.spkMov(end/2-imfs+1:end/2+imfs,dataIn.sigBoolMov,iTrGoodMov),2,'omitnan'),[],3,'omitnan')/sqrt(numel(iTrGoodMov))];
    spkMovMap = [spkMovMap; spkMovMapLocal]; %#ok<*AGROW>
    spkMov = [spkMov; spkMovLocal];
    %spkRew
    spkRewMapLocal = mean(dataIn.spkCor(end/2-imfs+1:end/2+imfs,dataIn.sigBoolRew,iTrGoodRew),3,'omitnan');
    spkRewLocal = [mean(spkRewMapLocal,2,'omitnan') std(mean(dataIn.spkCor(end/2-imfs+1:end/2+imfs,dataIn.sigBoolRew,iTrGoodRew),2,'omitnan'),[],3,'omitnan')/sqrt(numel(iTrGoodRew))];
    spkRewMap = [spkRewMap; spkRewMapLocal];
    spkRew = [spkRew; spkRewLocal];
    %spkUnder - reward cells
    spkUnderMapLocal = mean(dataIn.spkUnder(end/2-imfs+1:end/2+imfs,dataIn.sigBoolRew,iTrGoodUnder),3,'omitnan');
    spkUnderLocal = [mean(spkUnderMapLocal,2,'omitnan') std(mean(dataIn.spkUnder(end/2-imfs+1:end/2+imfs,dataIn.sigBoolRew,iTrGoodUnder),2,'omitnan'),[],3,'omitnan')/sqrt(numel(iTrGoodUnder))];
    spkUnderMap = [spkUnderMap; spkUnderMapLocal];
    spkUnder = [spkUnder; spkUnderLocal];
    %spkOver - reward cells
    spkOverMapLocal = mean(dataIn.spkOver(end/2-imfs+1:end/2+imfs,dataIn.sigBoolRew,iTrGoodOver),3,'omitnan');
    spkOverLocal = [mean(spkOverMapLocal,2,'omitnan') std(mean(dataIn.spkOver(end/2-imfs+1:end/2+imfs,dataIn.sigBoolRew,iTrGoodOver),2,'omitnan'),[],3,'omitnan')/sqrt(numel(iTrGoodOver))];
    spkOverMap = [spkOverMap; spkOverMapLocal];
    spkOver = [spkOver; spkOverLocal];
    %spkUnder - mov cells
    spkUnderMapLocal_mov = mean(dataIn.spkUnder(end/2-imfs+1:end/2+imfs,dataIn.sigBoolMov,iTrGoodUnder),3,'omitnan');
    spkUnderLocal_mov = [mean(spkUnderMapLocal_mov,2,'omitnan') std(mean(dataIn.spkUnder(end/2-imfs+1:end/2+imfs,dataIn.sigBoolMov,iTrGoodUnder),2,'omitnan'),[],3,'omitnan')/sqrt(numel(iTrGoodUnder))];
    spkUnderMap_mov = [spkUnderMap_mov; spkUnderMapLocal_mov];
    spkUnder_mov = [spkUnder_mov; spkUnderLocal_mov];
    %spkOver - mov cells
    spkOverMapLocal_mov = mean(dataIn.spkOver(end/2-imfs+1:end/2+imfs,dataIn.sigBoolMov,iTrGoodOver),3,'omitnan');
    spkOverLocal_mov = [mean(spkOverMapLocal_mov,2,'omitnan') std(mean(dataIn.spkOver(end/2-imfs+1:end/2+imfs,dataIn.sigBoolMov,iTrGoodOver),2,'omitnan'),[],3,'omitnan')/sqrt(numel(iTrGoodOver))];
    spkOverMap_mov = [spkOverMap_mov; spkOverMapLocal_mov];
    spkOver_mov = [spkOver_mov; spkOverLocal_mov];

    %xMov
    xMovLocal = [mean(dataIn.posMov(end/2-vrfs+1:end/2+vrfs,iTrGoodMov),2,'omitnan') ...
        std(dataIn.posMov(end/2-vrfs+1:end/2+vrfs,iTrGoodMov),[],2,'omitnan')];
    xMov = [xMov; xMovLocal];
    %vMov
    vMovLocal = [mean(dataIn.velMov(end/2-vrfs+1:end/2+vrfs,iTrGoodMov),2,'omitnan') ...
        std(dataIn.velMov(end/2-vrfs+1:end/2+vrfs,iTrGoodMov),[],2,'omitnan')];
    vMov = [vMov; vMovLocal];
    %xRew
    xRewLocal = [mean(dataIn.posCor(end/2-vrfs+1:end/2+vrfs,iTrGoodRew),2,'omitnan') ...
        std(dataIn.posCor(end/2-vrfs+1:end/2+vrfs,iTrGoodRew),[],2,'omitnan')];
    xRew = [xRew; xRewLocal];
    %vRew
    vRewLocal = [mean(dataIn.velCor(end/2-vrfs+1:end/2+vrfs,iTrGoodRew),2,'omitnan') ...
        std(dataIn.velCor(end/2-vrfs+1:end/2+vrfs,iTrGoodRew),[],2,'omitnan')];
    vRew = [vRew; vRewLocal];
    %xUnder
    xUnderLocal = [mean(dataIn.posUnder(end/2-vrfs+1:end/2+vrfs,iTrGoodUnder),2,'omitnan') ...
        std(dataIn.posUnder(end/2-vrfs+1:end/2+vrfs,iTrGoodUnder),[],2,'omitnan')];
    xUnder = [xUnder; xUnderLocal];
    %vUnder
    vUnderLocal = [mean(dataIn.velUnder(end/2-vrfs+1:end/2+vrfs,iTrGoodUnder),2,'omitnan') ...
        std(dataIn.velUnder(end/2-vrfs+1:end/2+vrfs,iTrGoodUnder),[],2,'omitnan')];
    vUnder = [vUnder; vUnderLocal];
    %xOver
    xOverLocal = [mean(dataIn.posOver(end/2-vrfs+1:end/2+vrfs,iTrGoodOver),2,'omitnan') ...
        std(dataIn.posOver(end/2-vrfs+1:end/2+vrfs,iTrGoodOver),[],2,'omitnan')];
    xOver = [xOver; xOverLocal];
    %vOver
    vOverLocal = [mean(dataIn.velOver(end/2-vrfs+1:end/2+vrfs,iTrGoodOver),2,'omitnan') ...
        std(dataIn.velOver(end/2-vrfs+1:end/2+vrfs,iTrGoodOver),[],2,'omitnan')];
    vOver = [vOver; vOverLocal];
    
end

spkMovMap = spkMovMap*imfs;
coeffs = pca(spkMovMap);
[~,iSort] = sort(coeffs(:,1),'descend');
spkMovMap_sorted = spkMovMap(:,iSort);

spkRewMap = spkRewMap*imfs;
coeffs = pca(spkRewMap);
[~,iSort] = sort(coeffs(:,1),'descend');
spkRewMap_sorted = spkRewMap(:,iSort);

spkUnderMap = spkUnderMap*imfs;
coeffs = pca(spkUnderMap);
[~,iSort] = sort(coeffs(:,1),'descend');
spkUnderMap_sorted = spkUnderMap(:,iSort);

spkOverMap = spkOverMap*imfs;
coeffs = pca(spkOverMap);
[~,iSort] = sort(coeffs(:,1),'descend');
spkOverMap_sorted = spkOverMap(:,iSort);

spkUnderMap_mov = spkUnderMap_mov*imfs;
coeffs = pca(spkUnderMap_mov);
[~,iSort] = sort(coeffs(:,1),'descend');
spkUnderMap_mov_sorted = spkUnderMap_mov(:,iSort);

spkOverMap_mov = spkOverMap_mov*imfs;
coeffs = pca(spkOverMap_mov);
[~,iSort] = sort(coeffs(:,1),'descend');
spkOverMap_mov_sorted = spkOverMap_mov(:,iSort);

spkMov = spkMov*imfs;
spkRew = spkRew*imfs;
spkUnder = spkUnder*imfs;
spkOver = spkOver*imfs;
spkUnder_mov = spkUnder_mov*imfs;
spkOver_mov = spkOver_mov*imfs;
clear iTr iTrGood spkMovLocal spkRewLocal spkMovMapLocal spkRewMapLocal vMovLocal vRewLocal xMovLocal xRewLocal spkMovMapLocal spkRewMapLocal

%% Plot movement signals:
t2p = (1:size(spkMov,1))'/imfs-1/(2*imfs);
tvr = (1:size(xMov,1))'/vrfs-1/(2*vrfs);

figOut(1) = figure;
subplot(4,1,1);
imagesc(t2p, 1:size(spkMovMap,2), spkMovMap');  colormap(viridis); clim([0 10]); box off;
xlim([0 20]); xticks(0:20);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xline(2:2:18,'-w',LineWidth= 1);
xline((2:2:20)-1,'--w',LineWidth= 0.5);

subplot(4,1,2);
imagesc(t2p, 1:size(spkMovMap_sorted,2), spkMovMap_sorted');  colormap(viridis); clim([0 10]); box off;
xlim([0 20]); xticks(0:20);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xline(2:2:18,'-w',LineWidth= 1);
xline((2:2:20)-1,'--w',LineWidth= 0.5);

subplot(4,1,3);
h1 = boundedline(t2p,spkMov(:,1),spkMov(:,2),'cmap',[0 0 0],'alpha');
ylim([0 12]); yticks(0:2:12);
xlim([0 20]); xticks(0:20);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xline(2:2:18,'-b',LineWidth= 1);
xline((2:2:20)-1,'--b',LineWidth= 0.5);

% subplot(4,1,3);
% h2 = boundedline(tvr,1-xMov(:,1),xMov(:,2),'cmap',[0 0 0],'alpha');
% ylim([-.1 2.1]); yticks([0:1/3:2]);
% xlim([0 20]); xticks(0:20);
% set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
% xline([2:2:18],'-b',LineWidth= 1);
% xline([2:2:20]-1,'--b',LineWidth= 0.5);

subplot(4,1,4);
h3 = boundedline(tvr,vMov(:,1),vMov(:,2),'cmap',[0 0 0],'alpha');
ylim([-2.5 17.5]); yticks([0:2.5:15]);
xlim([0 20]); xticks(0:20);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xline([2:2:18],'-b',LineWidth= 1);
xline([2:2:20]-1,'--b',LineWidth= 0.5);
suptitle(sprintf('Mouse: %s, date: %s, fov: %s, Movement',sets.name,sets.date,sets.fov));

%% Plot reward
figOut(2) = figure;
subplot(4,1,1);
imagesc(t2p, 1:size(spkRewMap,2), spkRewMap');  colormap(viridis); clim([0 20]); box off;
xlim([0 20]); xticks(0:20);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xline([2:2:18],'-w',LineWidth= 1);
xline([2:2:20]-1,'--w',LineWidth= 0.5);

subplot(4,1,2);
imagesc(t2p, 1:size(spkRewMap_sorted,2), spkRewMap_sorted');  colormap(viridis); clim([0 20]); box off;
xlim([0 20]); xticks(0:20);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xline([2:2:18],'-w',LineWidth= 1);
xline([2:2:20]-1,'--w',LineWidth= 0.5);

subplot(4,1,3);
h1 = boundedline(t2p,spkRew(:,1),spkRew(:,2),'cmap',[0 0 0],'alpha');
ylim([0 20]); yticks([0:2:20]);
xlim([0 20]); xticks(0:20);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xline([2:2:18],'-b',LineWidth= 1);
xline([2:2:20]-1,'--b',LineWidth= 0.5);

% subplot(4,1,3);
% h2 = boundedline(tvr,1-xRew(:,1),xRew(:,2),'cmap',[0 0 0],'alpha');
% ylim([-.1 2.1]); yticks([0:1/3:2]);
% xlim([0 20]); xticks(0:20);
% set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
% xline([2:2:18],'-b',LineWidth= 1);
% xline([2:2:20]-1,'--b',LineWidth= 0.5);

subplot(4,1,4);
h3 = boundedline(tvr,vRew(:,1),vRew(:,2),'cmap',[0 0 0],'alpha');
ylim([-2.5 17.5]); yticks([0:2.5:15]);
xlim([0 20]); xticks(0:20);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xline([2:2:18],'-b',LineWidth= 1);
xline([2:2:20]-1,'--b',LineWidth= 0.5);
suptitle(sprintf('Mouse: %s, date: %s, fov: %s, Reward',sets.name,sets.date,sets.fov));

%% Plot undershoots - reward PCs
figOut(3) = figure;
subplot(4,1,1);
imagesc(t2p, 1:size(spkUnderMap,2), spkUnderMap');  colormap(viridis); clim([0 20]); box off;
xlim([0 20]); xticks(0:20);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xline([2:2:18],'-w',LineWidth= 1);
xline([2:2:20]-1,'--w',LineWidth= 0.5);

subplot(4,1,2);
imagesc(t2p, 1:size(spkUnderMap_sorted,2), spkUnderMap_sorted');  colormap(viridis); clim([0 20]); box off;
xlim([0 20]); xticks(0:20);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xline([2:2:18],'-w',LineWidth= 1);
xline([2:2:20]-1,'--w',LineWidth= 0.5);

subplot(4,1,3);
h1 = boundedline(t2p,spkUnder(:,1),spkUnder(:,2),'cmap',[0 0 0],'alpha');
ylim([0 20]); yticks([0:2:20]);
xlim([0 20]); xticks(0:20);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xline([2:2:18],'-b',LineWidth= 1);
xline([2:2:20]-1,'--b',LineWidth= 0.5);

subplot(4,1,4);
h3 = boundedline(tvr,vUnder(:,1),vUnder(:,2),'cmap',[0 0 0],'alpha');
ylim([-2.5 17.5]); yticks([0:2.5:15]);
xlim([0 20]); xticks(0:20);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xline([2:2:18],'-b',LineWidth= 1);
xline([2:2:20]-1,'--b',LineWidth= 0.5);
suptitle(sprintf('Mouse: %s, date: %s, fov: %s, Undershoots - reward PCs',sets.name,sets.date,sets.fov));

%% Plot overshoots - reward PCs
figOut(4) = figure;
subplot(4,1,1);
imagesc(t2p, 1:size(spkOverMap,2), spkOverMap');  colormap(viridis); clim([0 20]); box off;
xlim([0 20]); xticks(0:20);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xline([2:2:18],'-w',LineWidth= 1);
xline([2:2:20]-1,'--w',LineWidth= 0.5);

subplot(4,1,2);
imagesc(t2p, 1:size(spkOverMap_sorted,2), spkOverMap_sorted');  colormap(viridis); clim([0 20]); box off;
xlim([0 20]); xticks(0:20);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xline([2:2:18],'-w',LineWidth= 1);
xline([2:2:20]-1,'--w',LineWidth= 0.5);

subplot(4,1,3);
h1 = boundedline(t2p,spkOver(:,1),spkOver(:,2),'cmap',[0 0 0],'alpha');
ylim([0 20]); yticks([0:2:20]);
xlim([0 20]); xticks(0:20);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xline([2:2:18],'-b',LineWidth= 1);
xline([2:2:20]-1,'--b',LineWidth= 0.5);

subplot(4,1,4);
h3 = boundedline(tvr,vOver(:,1),vOver(:,2),'cmap',[0 0 0],'alpha');
ylim([-2.5 17.5]); yticks([0:2.5:15]);
xlim([0 20]); xticks(0:20);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xline([2:2:18],'-b',LineWidth= 1);
xline([2:2:20]-1,'--b',LineWidth= 0.5);
suptitle(sprintf('Mouse: %s, date: %s, fov: %s, Overshoots - reward PCs',sets.name,sets.date,sets.fov));

%% Plot undershoots - mov PCs
figOut(5) = figure;
subplot(4,1,1);
imagesc(t2p, 1:size(spkUnderMap_mov,2), spkUnderMap_mov');  colormap(viridis); clim([0 20]); box off;
xlim([0 20]); xticks(0:20);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xline([2:2:18],'-w',LineWidth= 1);
xline([2:2:20]-1,'--w',LineWidth= 0.5);

subplot(4,1,2);
imagesc(t2p, 1:size(spkUnderMap_mov_sorted,2), spkUnderMap_mov_sorted');  colormap(viridis); clim([0 20]); box off;
xlim([0 20]); xticks(0:20);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xline([2:2:18],'-w',LineWidth= 1);
xline([2:2:20]-1,'--w',LineWidth= 0.5);

subplot(4,1,3);
h1 = boundedline(t2p,spkUnder_mov(:,1),spkUnder_mov(:,2),'cmap',[0 0 0],'alpha');
ylim([0 20]); yticks([0:2:20]);
xlim([0 20]); xticks(0:20);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xline([2:2:18],'-b',LineWidth= 1);
xline([2:2:20]-1,'--b',LineWidth= 0.5);

subplot(4,1,4);
h3 = boundedline(tvr,vUnder(:,1),vUnder(:,2),'cmap',[0 0 0],'alpha');
ylim([-2.5 17.5]); yticks([0:2.5:15]);
xlim([0 20]); xticks(0:20);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xline([2:2:18],'-b',LineWidth= 1);
xline([2:2:20]-1,'--b',LineWidth= 0.5);
suptitle(sprintf('Mouse: %s, date: %s, fov: %s, Undershoots - movement PCs',sets.name,sets.date,sets.fov));

%% Plot overshoots - mov PCs
figOut(6) = figure;
subplot(4,1,1);
imagesc(t2p, 1:size(spkOverMap_mov,2), spkOverMap_mov');  colormap(viridis); clim([0 20]); box off;
xlim([0 20]); xticks(0:20);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xline([2:2:18],'-w',LineWidth= 1);
xline([2:2:20]-1,'--w',LineWidth= 0.5);

subplot(4,1,2);
imagesc(t2p, 1:size(spkOverMap_mov_sorted,2), spkOverMap_mov_sorted');  colormap(viridis); clim([0 20]); box off;
xlim([0 20]); xticks(0:20);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xline([2:2:18],'-w',LineWidth= 1);
xline([2:2:20]-1,'--w',LineWidth= 0.5);

subplot(4,1,3);
h1 = boundedline(t2p,spkOver_mov(:,1),spkOver_mov(:,2),'cmap',[0 0 0],'alpha');
ylim([0 20]); yticks([0:2:20]);
xlim([0 20]); xticks(0:20);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xline([2:2:18],'-b',LineWidth= 1);
xline([2:2:20]-1,'--b',LineWidth= 0.5);

subplot(4,1,4);
h3 = boundedline(tvr,vOver(:,1),vOver(:,2),'cmap',[0 0 0],'alpha');
ylim([-2.5 17.5]); yticks([0:2.5:15]);
xlim([0 20]); xticks(0:20);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xline([2:2:18],'-b',LineWidth= 1);
xline([2:2:20]-1,'--b',LineWidth= 0.5);
suptitle(sprintf('Mouse: %s, date: %s, fov: %s, Overshoots - movement PCs',sets.name,sets.date,sets.fov));

%% Wrap structure and send back
dataOut = v2struct(spkMov,spkMovMap,xMov,vMov,spkRew,spkRewMap,xRew,vRew,spkUnder,spkUnderMap,spkUnder_mov,spkUnderMap_mov,xUnder,vUnder,spkOver,spkOverMap,spkOver_mov,spkOverMap_mov,xOver,vOver);


end