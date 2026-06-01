function [figout,matchPCs] = fig1_match_v144plot(matchFOVs,fovname, fs)
%% Loop through and assemble struct with matching rois and x-validated v4 trials

matchPCs = rmfield(matchFOVs,{'v1date','v4_1date','v4_Ndate','ilobv','isim2','v14rois','v44rois'});

for i = 1:numel(matchFOVs)
    matchPCs(i).col = repmat(matchPCs(i).col,size(matchPCs(i).v144rois,1),1);
    
    matchPCs(i).v144on(3).spkMap0 = mean(matchPCs(i).v144on(3).spk(:,matchPCs(i).v144rois(:,3),1:2:end),3);
    matchPCs(i).v144on(3).spkMap = mean(matchPCs(i).v144on(3).spk(:,matchPCs(i).v144rois(:,3),2:2:end),3);
    matchPCs(i).v144on(3).spkMu = mean(matchPCs(i).v144on(3).spkMap,2);
    matchPCs(i).v144on(3).lickMu = mean(matchPCs(i).v144on(3).lick(:,2:2:end),2);
    matchPCs(i).v144on(3).movMu = mean(matchPCs(i).v144on(3).mov(:,2:2:end),2);
    
    matchPCs(i).v144mov(3).spkMap0 = mean(matchPCs(i).v144mov(3).spk(:,matchPCs(i).v144rois(:,3),1:2:end),3);
    matchPCs(i).v144mov(3).spkMap = mean(matchPCs(i).v144mov(3).spk(:,matchPCs(i).v144rois(:,3),2:2:end),3);
    matchPCs(i).v144mov(3).spkMu = mean(matchPCs(i).v144mov(3).spkMap,2);
    matchPCs(i).v144mov(3).lickMu = mean(matchPCs(i).v144mov(3).lick(:,2:2:end),2);
    matchPCs(i).v144mov(3).movMu = mean(matchPCs(i).v144mov(3).mov(:,2:2:end),2);
    
    matchPCs(i).v144rewT(3).spkMap0 = mean(matchPCs(i).v144rewT(3).spk(:,matchPCs(i).v144rois(:,3),1:2:end),3);
    matchPCs(i).v144rewT(3).spkMap = mean(matchPCs(i).v144rewT(3).spk(:,matchPCs(i).v144rois(:,3),2:2:end),3);
    matchPCs(i).v144rewT(3).spkMu = mean(matchPCs(i).v144rewT(3).spkMap,2);
    matchPCs(i).v144rewT(3).lickMu = mean(matchPCs(i).v144rewT(3).lick(:,2:2:end),2);
    matchPCs(i).v144rewT(3).movMu = mean(matchPCs(i).v144rewT(3).mov(:,2:2:end),2);
    
    matchPCs(i).v144rewC(3).spkMap0 = mean(matchPCs(i).v144rewC(3).spk(:,matchPCs(i).v144rois(:,3),1:2:end),3);
    matchPCs(i).v144rewC(3).spkMap = mean(matchPCs(i).v144rewC(3).spk(:,matchPCs(i).v144rois(:,3),2:2:end),3);
    matchPCs(i).v144rewC(3).spkMu = mean(matchPCs(i).v144rewC(3).spkMap,2);
    matchPCs(i).v144rewC(3).lickMu = mean(matchPCs(i).v144rewC(3).lick(:,2:2:end),2);
    matchPCs(i).v144rewC(3).movMu = mean(matchPCs(i).v144rewC(3).mov(:,2:2:end),2);
    
    matchPCs(i).v144rewR(3).spkMap0 = mean(matchPCs(i).v144rewR(3).spk(:,matchPCs(i).v144rois(:,3),1:2:end),3);
    matchPCs(i).v144rewR(3).spkMap = mean(matchPCs(i).v144rewR(3).spk(:,matchPCs(i).v144rois(:,3),2:2:end),3);
    matchPCs(i).v144rewR(3).spkMu = mean(matchPCs(i).v144rewR(3).spkMap,2);
    matchPCs(i).v144rewR(3).lickMu = mean(matchPCs(i).v144rewR(3).lick(:,2:2:end),2);
    matchPCs(i).v144rewR(3).movMu = mean(matchPCs(i).v144rewR(3).mov(:,2:2:end),2);
    
    for j = [1,2]
        matchPCs(i).v144on(j).spkMap = mean(matchPCs(i).v144on(j).spk(:,matchPCs(i).v144rois(:,j),:),3);
        matchPCs(i).v144on(j).spkMu = mean(matchPCs(i).v144on(j).spkMap,2);
        matchPCs(i).v144on(j).lickMu = mean(matchPCs(i).v144on(j).lick,2);
        matchPCs(i).v144on(j).movMu = mean(matchPCs(i).v144on(j).mov,2);
        
        matchPCs(i).v144mov(j).spkMap = mean(matchPCs(i).v144mov(j).spk(:,matchPCs(i).v144rois(:,j),:),3);
        matchPCs(i).v144mov(j).spkMu = mean(matchPCs(i).v144mov(j).spkMap,2);
        matchPCs(i).v144mov(j).lickMu = mean(matchPCs(i).v144mov(j).lick,2);
        matchPCs(i).v144mov(j).movMu = mean(matchPCs(i).v144mov(j).mov,2);
        
        matchPCs(i).v144rewT(j).spkMap = mean(matchPCs(i).v144rewT(j).spk(:,matchPCs(i).v144rois(:,j),:),3);
        matchPCs(i).v144rewT(j).spkMu = mean(matchPCs(i).v144rewT(j).spkMap,2);
        matchPCs(i).v144rewT(j).lickMu = mean(matchPCs(i).v144rewT(j).lick,2);
        matchPCs(i).v144rewT(j).movMu = mean(matchPCs(i).v144rewT(j).mov,2);
        
        matchPCs(i).v144rewC(j).spkMap = mean(matchPCs(i).v144rewC(j).spk(:,matchPCs(i).v144rois(:,j),:),3);
        matchPCs(i).v144rewC(j).spkMu = mean(matchPCs(i).v144rewC(j).spkMap,2);
        matchPCs(i).v144rewC(j).lickMu = mean(matchPCs(i).v144rewC(j).lick,2);
        matchPCs(i).v144rewC(j).movMu = mean(matchPCs(i).v144rewC(j).mov,2);
        
        matchPCs(i).v144rewR(j).spkMap = mean(matchPCs(i).v144rewR(j).spk(:,matchPCs(i).v144rois(:,j),:),3);
        matchPCs(i).v144rewR(j).spkMu = mean(matchPCs(i).v144rewR(j).spkMap,2);
        matchPCs(i).v144rewR(j).lickMu = mean(matchPCs(i).v144rewR(j).lick,2);
        matchPCs(i).v144rewR(j).movMu = mean(matchPCs(i).v144rewR(j).mov,2);

    end
    
%     for j = 1:3
%         
%     matchPCs(i).on.spkMap = cat(1,matchPCs(i).v144on.spkMap(end/2-29:end/2+30,:));
    
end

%% Combine traces to get plottable vectors

v144on = cat(2,matchPCs.v144on); spkMap0on = cat(2,v144on.spkMap0);
coeff = pca(spkMap0on(end/2-fs.imfs/2+1:end/2+fs.imfs/2,:));
[~, iSortOn] = sort(coeff(:,1),'descend');
cols0 = cat(1,matchPCs.col); cols0 = cols0(iSortOn,:);
colsOn = zeros([1 size(cols0)]); colsOn(1,:,:) = cols0; colsOn = permute(repmat(colsOn,5,1,1),[2 1 3]); clear cols0
v144onPlot = fig1_match_v145trim(v144on,iSortOn,fs,0);

v144mov = cat(2,matchPCs.v144mov); spkMap0mov = cat(2,v144mov.spkMap0);
coeff = pca(spkMap0mov(end/2-fs.imfs/2+1:end/2+fs.imfs/2,:));
[~, iSortOn] = sort(coeff(:,1),'descend');
cols0 = cat(1,matchPCs.col); cols0 = cols0(iSortOn,:);
colsMov = zeros([1 size(cols0)]); colsMov(1,:,:) = cols0; colsMov = permute(repmat(colsMov,5,1,1),[2 1 3]); clear cols0
v144movPlot = fig1_match_v145trim(v144mov,iSortOn,fs,0);

v144rewT = cat(2,matchPCs.v144rewT); spkMap0rewT = cat(2,v144rewT.spkMap0);
coeff = pca(spkMap0rewT([end/2-fs.imfs/2+1:end/2+fs.imfs/2]+fs.rewDelay*fs.imfs,:));
[~, iSortOn] = sort(coeff(:,1),'descend');
cols0 = cat(1,matchPCs.col); cols0 = cols0(iSortOn,:);
colsRewT = zeros([1 size(cols0)]); colsRewT(1,:,:) = cols0; colsRewT = permute(repmat(colsRewT,5,1,1),[2 1 3]); clear cols0
v144rewTPlot = fig1_match_v145trim(v144rewT,iSortOn,fs,fs.rewDelay);

v144rewC = cat(2,matchPCs.v144rewC); spkMap0rewC = cat(2,v144rewC.spkMap0);
coeff = pca(spkMap0rewC(end/2-fs.imfs/2+1:end/2+fs.imfs/2,:));
[~, iSortOn] = sort(coeff(:,1),'descend');
cols0 = cat(1,matchPCs.col); cols0 = cols0(iSortOn,:);
colsRewC = zeros([1 size(cols0)]); colsRewC(1,:,:) = cols0; colsRewC = permute(repmat(colsRewC,5,1,1),[2 1 3]); clear cols0
v144rewCPlot = fig1_match_v145trim(v144rewC,iSortOn,fs,0);

v144rewR = cat(2,matchPCs.v144rewR); spkMap0rewR = cat(2,v144rewR.spkMap0);
coeff = pca(spkMap0rewR(end/2-fs.imfs/2+1:end/2+fs.imfs/2,:));
[~, iSortOn] = sort(coeff(:,1),'descend');
cols0 = cat(1,matchPCs.col); cols0 = cols0(iSortOn,:);
colsRewR = zeros([1 size(cols0)]); colsRewR(1,:,:) = cols0; colsRewR = permute(repmat(colsRewR,5,1,1),[2 1 3]); clear cols0
v144rewRPlot = fig1_match_v145trim(v144rewR,iSortOn,fs,0);


%% Plot heatmaps and mean responses for behaviour

nRois = size(colsOn,1);
t2p = ((1:fs.imfs*6)/fs.imfs)-1/(fs.imfs*2);
tvr = ((1:fs.behavfs*6)/fs.behavfs)-1/(fs.behavfs*2);

figout(1) = figure;
subplot(6,15,[1:4,16:19,31:34]);
imagesc(t2p,[1:nRois],v144onPlot.spkMap',[0 10]);
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks(0:.5:6); xticklabels({[],[],'Naive',[],[],[],'Learning',[],[],[],'Trained',[],[]}); xtickangle(0);
yticks(0.5:50:nRois); yticklabels([1,50:50:nRois]);
hold on;
xline(1,'--w','LineWidth',1); xline(3,'--w','LineWidth',1); xline(5,'--w','LineWidth',1)
xline(2,'-w','LineWidth',2); xline(4,'-w','LineWidth',2);
title('Trial onset');
ylabel([fovname,' PC dendrite number']);
subplot(6,15,[5 20 35]);
imagesc(colsOn); axis image; axis off;


subplot(6,15,[1:4,16:19,31:34]+5);
imagesc(t2p,[1:nRois],v144movPlot.spkMap',[0 10]);
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks(0:.5:6); xticklabels({[],[],'Naive',[],[],[],'Learning',[],[],[],'Trained',[],[]}); xtickangle(0);
yticks(0.5:50:nRois); yticklabels([1,50:50:nRois]);
hold on;
xline(1,'--w','LineWidth',1); xline(3,'--w','LineWidth',1); xline(5,'--w','LineWidth',1)
xline(2,'-w','LineWidth',2); xline(4,'-w','LineWidth',2);
title('Movement onset');
subplot(6,15,[5 20 35]+5);
imagesc(colsMov); axis image; axis off;


subplot(6,15,[1:4,16:19,31:34]+10);
imagesc(t2p,[1:nRois],v144rewTPlot.spkMap',[0 10]);
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks(0:.5:6); xticklabels({[],[],'Naive',[],[],[],'Learning',[],[],[],'Trained',[],[]}); xtickangle(0);
yticks(0.5:50:nRois); yticklabels([1,50:50:nRois]);
hold on;
xline(1,'--w','LineWidth',1); xline(3,'--w','LineWidth',1); xline(5,'--w','LineWidth',1)
xline(2,'-w','LineWidth',3); xline(4,'-w','LineWidth',3);
title('Trial reward');
subplot(6,15,[5 20 35]+10);
imagesc(colsRewT); axis image; axis off;

subplot(6,15,[46:49]);
h = boundedline(t2p,mean(v144onPlot.spkMu,2), std(v144onPlot.spkMu,[],2),'-k');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(5,'--b','LineWidth',1)
xline(2,'-b','LineWidth',2); xline(4,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 6 0 12]);
xticks(0:.5:6); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:3:12]);
ylabel('Pop. response (spk/s)');

subplot(6,15,[46:49]+5);
h = boundedline(t2p,mean(v144movPlot.spkMu,2), std(v144movPlot.spkMu,[],2),'-k');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(5,'--b','LineWidth',1)
xline(2,'-b','LineWidth',2); xline(4,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 6 0 12]);
xticks(0:.5:6); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:3:12]);

subplot(6,15,[46:49]+10);
h = boundedline(t2p,mean(v144rewTPlot.spkMu,2), std(v144rewTPlot.spkMu,[],2),'-k');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(5,'--b','LineWidth',1)
xline(2,'-b','LineWidth',2); xline(4,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 6 0 12]);
xticks(0:.5:6); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
xticks([-1:0.5:1]); xticklabels({});
yticks([0:3:12]);

subplot(6,15,[61:64]);
h = boundedline(tvr,mean(v144onPlot.lickMu,2), std(v144onPlot.lickMu,[],2),'-k');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(5,'--b','LineWidth',1)
xline(2,'-b','LineWidth',2); xline(4,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 6 0 1]);
xticks(0:.5:6); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:.25:1]);
ylabel('Licking');

subplot(6,15,[61:64]+5);
h = boundedline(tvr,mean(v144movPlot.lickMu,2), std(v144movPlot.lickMu,[],2),'-k');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(5,'--b','LineWidth',1)
xline(2,'-b','LineWidth',2); xline(4,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 6 0 1]);
xticks(0:.5:6); xticklabels({});  % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:.25:1]);

subplot(6,15,[61:64]+10);
h = boundedline(tvr,mean(v144rewTPlot.lickMu,2), std(v144rewTPlot.lickMu,[],2),'-k');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(5,'--b','LineWidth',1)
xline(2,'-b','LineWidth',2); xline(4,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 6 0 1]);
xticks(0:.5:6); xticklabels({});  % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:.25:1]);

subplot(6,15,[76:79]);
h = boundedline(tvr,mean(v144onPlot.movMu,2), std(v144onPlot.movMu,[],2),'-k');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(5,'--b','LineWidth',1)
xline(2,'-b','LineWidth',2); xline(4,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 6 -2.5 12.5]);
xticks(0:.5:6); xticklabels({});  % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([-2.5:2.5:12.5]);
ylabel('Movement (cm/s)');

subplot(6,15,[76:79]+5);
h = boundedline(tvr,mean(v144movPlot.movMu,2), std(v144movPlot.movMu,[],2),'-k');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(5,'--b','LineWidth',1)
xline(2,'-b','LineWidth',2); xline(4,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 6 -2.5 12.5]);
xticks(0:.5:6); xticklabels(repmat([-1:.5:0.5],1,3)); xtickangle(0);
xlabel('Time (s)');
yticks([-2.5:2.5:12.5]);

subplot(6,15,[76:79]+10);
h = boundedline(tvr,mean(v144rewTPlot.movMu,2), std(v144rewTPlot.movMu,[],2),'-k');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(5,'--b','LineWidth',1)
xline(2,'-b','LineWidth',2); xline(4,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 6 -2.5 12.5]);
xticks(0:.5:6); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([-2.5:2.5:12.5]);

%% Plot heatmaps and mean responses for rewards

figout(2) = figure;
subplot(6,15,[1:4,16:19,31:34]);
imagesc(t2p,[1:nRois],v144rewRPlot.spkMap',[0 10]);
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks(0:.5:6); xticklabels({[],[],'Naive',[],[],[],'Learning',[],[],[],'Trained',[],[]}); xtickangle(0);
yticks(0.5:50:nRois); yticklabels([1,50:50:nRois]);
hold on;
xline(1,'--w','LineWidth',1); xline(3,'--w','LineWidth',1); xline(5,'--w','LineWidth',1)
xline(2,'-w','LineWidth',3); xline(4,'-w','LineWidth',3);
title('Random reward');
ylabel([fovname,' PC dendrite number']);
subplot(6,15,[5 20 35]);
imagesc(colsRewR); axis image; axis off;

subplot(6,15,[1:4,16:19,31:34]+5);
imagesc(t2p,[1:nRois],v144rewTPlot.spkMap',[0 10]);
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks(0:.5:6); xticklabels({[],[],'Naive',[],[],[],'Learning',[],[],[],'Trained',[],[]}); xtickangle(0);
yticks(0.5:50:nRois); yticklabels([1,50:50:nRois]);
hold on;
xline(1,'--w','LineWidth',1); xline(3,'--w','LineWidth',1); xline(5,'--w','LineWidth',1)
xline(2,'-w','LineWidth',3); xline(4,'-w','LineWidth',3);
title('Trial reward');
subplot(6,15,[5 20 35]+5);
imagesc(colsRewT); axis image; axis off;

subplot(6,15,[1:4,16:19,31:34]+10);
imagesc(t2p,[1:nRois],v144rewCPlot.spkMap',[0 10]);
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks(0:.5:6); xticklabels({[],[],'Naive',[],[],[],'Learning',[],[],[],'Trained',[],[]}); xtickangle(0);
yticks(0.5:50:nRois); yticklabels([1,50:50:nRois]);
hold on;
xline(1,'--w','LineWidth',1); xline(3,'--w','LineWidth',1); xline(5,'--w','LineWidth',1)
xline(2,'-w','LineWidth',3); xline(4,'-w','LineWidth',3);
title('Cued reward');
subplot(6,15,[5 20 35]+10);
imagesc(colsRewC); axis image; axis off;

subplot(6,15,[46:49]);
h = boundedline(t2p,mean(v144rewRPlot.spkMu,2), std(v144rewRPlot.spkMu,[],2),'-k');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(5,'--b','LineWidth',1)
xline(2,'-b','LineWidth',2); xline(4,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 6 0 12]);
xticks(0:.5:6); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:3:12]);
ylabel('Pop. response (spk/s)');

subplot(6,15,[46:49]+5);
h = boundedline(t2p,mean(v144rewTPlot.spkMu,2), std(v144rewTPlot.spkMu,[],2),'-k');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(5,'--b','LineWidth',1)
xline(2,'-b','LineWidth',2); xline(4,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 6 0 12]);
xticks(0:.5:6); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:3:12]);

subplot(6,15,[46:49]+10);
h = boundedline(t2p,mean(v144rewCPlot.spkMu,2), std(v144rewCPlot.spkMu,[],2),'-k');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(5,'--b','LineWidth',1)
xline(2,'-b','LineWidth',2); xline(4,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 6 0 12]);
xticks(0:.5:6); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
xticks([-1:0.5:1]); xticklabels({});
yticks([0:3:12]);

subplot(6,15,[61:64]);
h = boundedline(tvr,mean(v144rewRPlot.lickMu,2), std(v144rewRPlot.lickMu,[],2),'-k');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(5,'--b','LineWidth',1)
xline(2,'-b','LineWidth',2); xline(4,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 6 0 1]);
xticks(0:.5:6); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:.25:1]);
ylabel('Licking');

subplot(6,15,[61:64]+5);
h = boundedline(tvr,mean(v144rewTPlot.lickMu,2), std(v144rewTPlot.lickMu,[],2),'-k');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(5,'--b','LineWidth',1)
xline(2,'-b','LineWidth',2); xline(4,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 6 0 1]);
xticks(0:.5:6); xticklabels({});  % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:.25:1]);

subplot(6,15,[61:64]+10);
h = boundedline(tvr,mean(v144rewCPlot.lickMu,2), std(v144rewCPlot.lickMu,[],2),'-k');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(5,'--b','LineWidth',1)
xline(2,'-b','LineWidth',2); xline(4,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 6 0 1]);
xticks(0:.5:6); xticklabels({});  % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:.25:1]);

subplot(6,15,[76:79]);
h = boundedline(tvr,mean(v144rewRPlot.movMu,2), std(v144rewRPlot.movMu,[],2),'-k');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(5,'--b','LineWidth',1)
xline(2,'-b','LineWidth',2); xline(4,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 6 -2.5 12.5]);
xticks(0:.5:6); xticklabels({});  % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([-2.5:2.5:12.5]);
ylabel('Movement (cm/s)');

subplot(6,15,[76:79]+5);
h = boundedline(tvr,mean(v144rewTPlot.movMu,2), std(v144rewTPlot.movMu,[],2),'-k');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(5,'--b','LineWidth',1)
xline(2,'-b','LineWidth',2); xline(4,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 6 -2.5 12.5]);
xticks(0:.5:6); xticklabels(repmat([-1:.5:0.5],1,3));
xlabel('Time (s)');
yticks([-2.5:2.5:12.5]);

subplot(6,15,[76:79]+10);
h = boundedline(tvr,mean(v144rewCPlot.movMu,2), std(v144rewCPlot.movMu,[],2),'-k');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(5,'--b','LineWidth',1)
xline(2,'-b','LineWidth',2); xline(4,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 6 -2.5 12.5]);
xticks(0:.5:6); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([-2.5:2.5:12.5]);

%% Summary plots - behaviour

figout(3) = figure;
subplot(3,3,1); hold on;
h = boundedline(t2p(1:end/3),mean(v144onPlot.spkMu(1:end/3,:),2), std(v144onPlot.spkMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(t2p(1:end/3),mean(v144onPlot.spkMu(end/3*2+1:end,:),2), std(v144onPlot.spkMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(t2p(1:end/3),mean(v144onPlot.spkMu(end/3+1:end/3*2,:),2), std(v144onPlot.spkMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 0 12]);
xticks(0:.5:2); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:3:12]);
ylabel('Pop. response (spk/s)');
title('Trial onset');

subplot(3,3,2); hold on;
h = boundedline(t2p(1:end/3),mean(v144movPlot.spkMu(1:end/3,:),2), std(v144movPlot.spkMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(t2p(1:end/3),mean(v144movPlot.spkMu(end/3*2+1:end,:),2), std(v144movPlot.spkMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(t2p(1:end/3),mean(v144movPlot.spkMu(end/3+1:end/3*2,:),2), std(v144movPlot.spkMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 0 12]);
xticks(0:.5:2); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:3:12]);
title('Movement onset');

subplot(3,3,3); hold on;
h = boundedline(t2p(1:end/3),mean(v144rewTPlot.spkMu(1:end/3,:),2), std(v144rewTPlot.spkMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(t2p(1:end/3),mean(v144rewTPlot.spkMu(end/3*2+1:end,:),2), std(v144rewTPlot.spkMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(t2p(1:end/3),mean(v144rewTPlot.spkMu(end/3+1:end/3*2,:),2), std(v144rewTPlot.spkMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 0 12]);
xticks(0:.5:2); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:3:12]);
title('Trial reward');

subplot(3,3,4); hold on;
h = boundedline(tvr(1:end/3),mean(v144onPlot.lickMu(1:end/3,:),2), std(v144onPlot.lickMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v144onPlot.lickMu(end/3*2+1:end,:),2), std(v144onPlot.lickMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v144onPlot.lickMu(end/3+1:end/3*2,:),2), std(v144onPlot.lickMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 0 1]);
xticks(0:.5:2); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:.25:1]);
ylabel('Licking');

subplot(3,3,5); hold on;
h = boundedline(tvr(1:end/3),mean(v144movPlot.lickMu(1:end/3,:),2), std(v144movPlot.lickMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v144movPlot.lickMu(end/3*2+1:end,:),2), std(v144movPlot.lickMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v144movPlot.lickMu(end/3+1:end/3*2,:),2), std(v144movPlot.lickMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 0 1]);
xticks(0:.5:2); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:.25:1]);

subplot(3,3,6); hold on;
h = boundedline(tvr(1:end/3),mean(v144rewTPlot.lickMu(1:end/3,:),2), std(v144rewTPlot.lickMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v144rewTPlot.lickMu(end/3*2+1:end,:),2), std(v144rewTPlot.lickMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v144rewTPlot.lickMu(end/3+1:end/3*2,:),2), std(v144rewTPlot.lickMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 0 1]);
xticks(0:.5:2); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:.25:1]);

subplot(3,3,7); hold on;
h = boundedline(tvr(1:end/3),mean(v144onPlot.movMu(1:end/3,:),2), std(v144onPlot.movMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v144onPlot.movMu(end/3*2+1:end,:),2), std(v144onPlot.movMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v144onPlot.movMu(end/3+1:end/3*2,:),2), std(v144onPlot.movMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 -2.5 12.5]);
xticks(0:.5:2); xticklabels([-1:.5:1]);
yticks([-2.5:2.5:12.5]);
ylabel('Movement (cm/s)');

subplot(3,3,8); hold on;
h = boundedline(tvr(1:end/3),mean(v144movPlot.movMu(1:end/3,:),2), std(v144movPlot.movMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v144movPlot.movMu(end/3*2+1:end,:),2), std(v144movPlot.movMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v144movPlot.movMu(end/3+1:end/3*2,:),2), std(v144movPlot.movMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 -2.5 12.5]);
xticks(0:.5:2); xticklabels([-1:.5:1]);
yticks([-2.5:2.5:12.5]);
xlabel('Time (s)');

subplot(3,3,9); hold on;
h = boundedline(tvr(1:end/3),mean(v144rewTPlot.movMu(1:end/3,:),2), std(v144rewTPlot.movMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v144rewTPlot.movMu(end/3*2+1:end,:),2), std(v144rewTPlot.movMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v144rewTPlot.movMu(end/3+1:end/3*2,:),2), std(v144rewTPlot.movMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 -2.5 12.5]);
xticks(0:.5:2); xticklabels([-1:.5:1]);
yticks([-2.5:2.5:12.5]);

%% Summary plots - reward

figout(4) = figure;
subplot(3,3,1); hold on;
h = boundedline(t2p(1:end/3),mean(v144rewRPlot.spkMu(1:end/3,:),2), std(v144rewRPlot.spkMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(t2p(1:end/3),mean(v144rewRPlot.spkMu(end/3*2+1:end,:),2), std(v144rewRPlot.spkMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(t2p(1:end/3),mean(v144rewRPlot.spkMu(end/3+1:end/3*2,:),2), std(v144rewRPlot.spkMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 0 15]);
xticks(0:.5:2); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:3:15]);
ylabel('Pop. response (spk/s)');
title('Random reward');

subplot(3,3,2); hold on;
h = boundedline(t2p(1:end/3),mean(v144rewTPlot.spkMu(1:end/3,:),2), std(v144rewTPlot.spkMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(t2p(1:end/3),mean(v144rewTPlot.spkMu(end/3*2+1:end,:),2), std(v144rewTPlot.spkMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(t2p(1:end/3),mean(v144rewTPlot.spkMu(end/3+1:end/3*2,:),2), std(v144rewTPlot.spkMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 0 15]);
xticks(0:.5:2); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:3:15]);
title('Trial reward');

subplot(3,3,3); hold on;
h = boundedline(t2p(1:end/3),mean(v144rewCPlot.spkMu(1:end/3,:),2), std(v144rewCPlot.spkMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(t2p(1:end/3),mean(v144rewCPlot.spkMu(end/3*2+1:end,:),2), std(v144rewCPlot.spkMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(t2p(1:end/3),mean(v144rewCPlot.spkMu(end/3+1:end/3*2,:),2), std(v144rewCPlot.spkMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 0 15]);
xticks(0:.5:2); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:3:15]);
title('Cued reward');

subplot(3,3,4); hold on;
h = boundedline(tvr(1:end/3),mean(v144rewRPlot.lickMu(1:end/3,:),2), std(v144rewRPlot.lickMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v144rewRPlot.lickMu(end/3*2+1:end,:),2), std(v144rewRPlot.lickMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v144rewRPlot.lickMu(end/3+1:end/3*2,:),2), std(v144rewRPlot.lickMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 0 1]);
xticks(0:.5:2); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:.25:1]);
ylabel('Licking');

subplot(3,3,5); hold on;
h = boundedline(tvr(1:end/3),mean(v144rewTPlot.lickMu(1:end/3,:),2), std(v144rewTPlot.lickMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v144rewTPlot.lickMu(end/3*2+1:end,:),2), std(v144rewTPlot.lickMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v144rewTPlot.lickMu(end/3+1:end/3*2,:),2), std(v144rewTPlot.lickMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 0 1]);
xticks(0:.5:2); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:.25:1]);

subplot(3,3,6); hold on;
h = boundedline(tvr(1:end/3),mean(v144rewCPlot.lickMu(1:end/3,:),2), std(v144rewCPlot.lickMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v144rewCPlot.lickMu(end/3*2+1:end,:),2), std(v144rewCPlot.lickMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v144rewCPlot.lickMu(end/3+1:end/3*2,:),2), std(v144rewCPlot.lickMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 0 1]);
xticks(0:.5:2); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:.25:1]);

subplot(3,3,7); hold on;
h = boundedline(tvr(1:end/3),mean(v144rewRPlot.movMu(1:end/3,:),2), std(v144rewRPlot.movMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v144rewRPlot.movMu(end/3*2+1:end,:),2), std(v144rewRPlot.movMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v144rewRPlot.movMu(end/3+1:end/3*2,:),2), std(v144rewRPlot.movMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 -2.5 12.5]);
xticks(0:.5:2); xticklabels([-1:.5:1]);
yticks([-2.5:2.5:12.5]);
ylabel('Movement (cm/s)');

subplot(3,3,8); hold on;
h = boundedline(tvr(1:end/3),mean(v144rewTPlot.movMu(1:end/3,:),2), std(v144rewTPlot.movMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v144rewTPlot.movMu(end/3*2+1:end,:),2), std(v144rewTPlot.movMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v144rewTPlot.movMu(end/3+1:end/3*2,:),2), std(v144rewTPlot.movMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 -2.5 12.5]);
xticks(0:.5:2); xticklabels([-1:.5:1]);
yticks([-2.5:2.5:12.5]);
xlabel('Time (s)');

subplot(3,3,9); hold on;
h = boundedline(tvr(1:end/3),mean(v144rewCPlot.movMu(1:end/3,:),2), std(v144rewCPlot.movMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v144rewCPlot.movMu(end/3*2+1:end,:),2), std(v144rewCPlot.movMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v144rewCPlot.movMu(end/3+1:end/3*2,:),2), std(v144rewCPlot.movMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 -2.5 12.5]);
xticks(0:.5:2); xticklabels([-1:.5:1]);
yticks([-2.5:2.5:12.5]);

end