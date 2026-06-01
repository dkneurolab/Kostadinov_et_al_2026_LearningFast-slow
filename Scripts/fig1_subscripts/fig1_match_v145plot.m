function [figout,matchPCs] = fig1_match_v145plot(matchFOVs,fovname)
%% Loop through and assemble struct with matching rois and x-validated v4 trials

matchPCs = rmfield(matchFOVs,{'v1date','v4_Ndate','v5date','ilobv','isim2','v14rois','v45rois'});

for i = 1:numel(matchFOVs)
    matchPCs(i).col = repmat(matchPCs(i).col,size(matchPCs(i).v145rois,1),1);
    
    matchPCs(i).v145on(2).spkMap0 = mean(matchPCs(i).v145on(2).spk(:,matchPCs(i).v145rois(:,2),1:2:end),3);
    matchPCs(i).v145on(2).spkMap = mean(matchPCs(i).v145on(2).spk(:,matchPCs(i).v145rois(:,2),2:2:end),3);
    matchPCs(i).v145on(2).spkMu = mean(matchPCs(i).v145on(2).spkMap,2);
    matchPCs(i).v145on(2).lickMu = mean(matchPCs(i).v145on(2).lick(:,2:2:end),2);
    matchPCs(i).v145on(2).movMu = mean(matchPCs(i).v145on(2).mov(:,2:2:end),2);
    
    matchPCs(i).v145mov(2).spkMap0 = mean(matchPCs(i).v145mov(2).spk(:,matchPCs(i).v145rois(:,2),1:2:end),3);
    matchPCs(i).v145mov(2).spkMap = mean(matchPCs(i).v145mov(2).spk(:,matchPCs(i).v145rois(:,2),2:2:end),3);
    matchPCs(i).v145mov(2).spkMu = mean(matchPCs(i).v145mov(2).spkMap,2);
    matchPCs(i).v145mov(2).lickMu = mean(matchPCs(i).v145mov(2).lick(:,2:2:end),2);
    matchPCs(i).v145mov(2).movMu = mean(matchPCs(i).v145mov(2).mov(:,2:2:end),2);
    
    matchPCs(i).v145rewT(2).spkMap0 = mean(matchPCs(i).v145rewT(2).spk(:,matchPCs(i).v145rois(:,2),1:2:end),3);
    matchPCs(i).v145rewT(2).spkMap = mean(matchPCs(i).v145rewT(2).spk(:,matchPCs(i).v145rois(:,2),2:2:end),3);
    matchPCs(i).v145rewT(2).spkMu = mean(matchPCs(i).v145rewT(2).spkMap,2);
    matchPCs(i).v145rewT(2).lickMu = mean(matchPCs(i).v145rewT(2).lick(:,2:2:end),2);
    matchPCs(i).v145rewT(2).movMu = mean(matchPCs(i).v145rewT(2).mov(:,2:2:end),2);
    
    matchPCs(i).v145rewC(2).spkMap0 = mean(matchPCs(i).v145rewC(2).spk(:,matchPCs(i).v145rois(:,2),1:2:end),3);
    matchPCs(i).v145rewC(2).spkMap = mean(matchPCs(i).v145rewC(2).spk(:,matchPCs(i).v145rois(:,2),2:2:end),3);
    matchPCs(i).v145rewC(2).spkMu = mean(matchPCs(i).v145rewC(2).spkMap,2);
    matchPCs(i).v145rewC(2).lickMu = mean(matchPCs(i).v145rewC(2).lick(:,2:2:end),2);
    matchPCs(i).v145rewC(2).movMu = mean(matchPCs(i).v145rewC(2).mov(:,2:2:end),2);
    
    matchPCs(i).v145rewR(2).spkMap0 = mean(matchPCs(i).v145rewR(2).spk(:,matchPCs(i).v145rois(:,2),1:2:end),3);
    matchPCs(i).v145rewR(2).spkMap = mean(matchPCs(i).v145rewR(2).spk(:,matchPCs(i).v145rois(:,2),2:2:end),3);
    matchPCs(i).v145rewR(2).spkMu = mean(matchPCs(i).v145rewR(2).spkMap,2);
    matchPCs(i).v145rewR(2).lickMu = mean(matchPCs(i).v145rewR(2).lick(:,2:2:end),2);
    matchPCs(i).v145rewR(2).movMu = mean(matchPCs(i).v145rewR(2).mov(:,2:2:end),2);
    
    for j = [1,3]
        matchPCs(i).v145on(j).spkMap = mean(matchPCs(i).v145on(j).spk(:,matchPCs(i).v145rois(:,j),:),3);
        matchPCs(i).v145on(j).spkMu = mean(matchPCs(i).v145on(j).spkMap,2);
        matchPCs(i).v145on(j).lickMu = mean(matchPCs(i).v145on(j).lick,2);
        matchPCs(i).v145on(j).movMu = mean(matchPCs(i).v145on(j).mov,2);
        
        matchPCs(i).v145mov(j).spkMap = mean(matchPCs(i).v145mov(j).spk(:,matchPCs(i).v145rois(:,j),:),3);
        matchPCs(i).v145mov(j).spkMu = mean(matchPCs(i).v145mov(j).spkMap,2);
        matchPCs(i).v145mov(j).lickMu = mean(matchPCs(i).v145mov(j).lick,2);
        matchPCs(i).v145mov(j).movMu = mean(matchPCs(i).v145mov(j).mov,2);
        
        matchPCs(i).v145rewT(j).spkMap = mean(matchPCs(i).v145rewT(j).spk(:,matchPCs(i).v145rois(:,j),:),3);
        matchPCs(i).v145rewT(j).spkMu = mean(matchPCs(i).v145rewT(j).spkMap,2);
        matchPCs(i).v145rewT(j).lickMu = mean(matchPCs(i).v145rewT(j).lick,2);
        matchPCs(i).v145rewT(j).movMu = mean(matchPCs(i).v145rewT(j).mov,2);
        
        matchPCs(i).v145rewC(j).spkMap = mean(matchPCs(i).v145rewC(j).spk(:,matchPCs(i).v145rois(:,j),:),3);
        matchPCs(i).v145rewC(j).spkMu = mean(matchPCs(i).v145rewC(j).spkMap,2);
        matchPCs(i).v145rewC(j).lickMu = mean(matchPCs(i).v145rewC(j).lick,2);
        matchPCs(i).v145rewC(j).movMu = mean(matchPCs(i).v145rewC(j).mov,2);
        
        matchPCs(i).v145rewR(j).spkMap = mean(matchPCs(i).v145rewR(j).spk(:,matchPCs(i).v145rois(:,j),:),3);
        matchPCs(i).v145rewR(j).spkMu = mean(matchPCs(i).v145rewR(j).spkMap,2);
        matchPCs(i).v145rewR(j).lickMu = mean(matchPCs(i).v145rewR(j).lick,2);
        matchPCs(i).v145rewR(j).movMu = mean(matchPCs(i).v145rewR(j).mov,2);

    end
    
%     for j = 1:3
%         
%     matchPCs(i).on.spkMap = cat(1,matchPCs(i).v145on.spkMap(end/2-29:end/2+30,:));
    
end

%% Combine traces to get plottable vectors

v145on = cat(2,matchPCs.v145on); spkMap0on = cat(2,v145on.spkMap0);
coeff = pca(spkMap0on(end/2-14:end/2+15,:));
[~, iSortOn] = sort(coeff(:,1),'descend');
cols0 = cat(1,matchPCs.col); cols0 = cols0(iSortOn,:);
colsOn = zeros([1 size(cols0)]); colsOn(1,:,:) = cols0; colsOn = permute(repmat(colsOn,5,1,1),[2 1 3]); clear cols0
v145onPlot = fig1_match_v145trim(v145on,iSortOn);

v145mov = cat(2,matchPCs.v145mov); spkMap0mov = cat(2,v145mov.spkMap0);
coeff = pca(spkMap0mov(end/2-14:end/2+15,:));
[~, iSortOn] = sort(coeff(:,1),'descend');
cols0 = cat(1,matchPCs.col); cols0 = cols0(iSortOn,:);
colsMov = zeros([1 size(cols0)]); colsMov(1,:,:) = cols0; colsMov = permute(repmat(colsMov,5,1,1),[2 1 3]); clear cols0
v145movPlot = fig1_match_v145trim(v145mov,iSortOn);

v145rewT = cat(2,matchPCs.v145rewT); spkMap0rewT = cat(2,v145rewT.spkMap0);
coeff = pca(spkMap0rewT(end/2-14:end/2+15,:));
[~, iSortOn] = sort(coeff(:,1),'descend');
cols0 = cat(1,matchPCs.col); cols0 = cols0(iSortOn,:);
colsRewT = zeros([1 size(cols0)]); colsRewT(1,:,:) = cols0; colsRewT = permute(repmat(colsRewT,5,1,1),[2 1 3]); clear cols0
v145rewTPlot = fig1_match_v145trim(v145rewT,iSortOn);

v145rewC = cat(2,matchPCs.v145rewC); spkMap0rewC = cat(2,v145rewC.spkMap0);
coeff = pca(spkMap0rewC(end/2-14:end/2+15,:));
[~, iSortOn] = sort(coeff(:,1),'descend');
cols0 = cat(1,matchPCs.col); cols0 = cols0(iSortOn,:);
colsRewC = zeros([1 size(cols0)]); colsRewC(1,:,:) = cols0; colsRewC = permute(repmat(colsRewC,5,1,1),[2 1 3]); clear cols0
v145rewCPlot = fig1_match_v145trim(v145rewC,iSortOn);

v145rewR = cat(2,matchPCs.v145rewR); spkMap0rewR = cat(2,v145rewR.spkMap0);
coeff = pca(spkMap0rewR(end/2-14:end/2+15,:));
[~, iSortOn] = sort(coeff(:,1),'descend');
cols0 = cat(1,matchPCs.col); cols0 = cols0(iSortOn,:);
colsRewR = zeros([1 size(cols0)]); colsRewR(1,:,:) = cols0; colsRewR = permute(repmat(colsRewR,5,1,1),[2 1 3]); clear cols0
v145rewRPlot = fig1_match_v145trim(v145rewR,iSortOn);


%% Plot heatmaps and mean responses for behaviour

nRois = size(colsOn,1);
t2p = ((1:180)/30)-1/60;
tvr = ((1:600)/100)-1/200;

figout(1) = figure;
subplot(6,15,[1:4,16:19,31:34]);
imagesc(t2p,[1:nRois],v145onPlot.spkMap',[0 10]);
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks(0:.5:6); xticklabels({[],[],'Naive',[],[],[],'Trained',[],[],[],'Adaptation',[],[]}); xtickangle(0);
yticks(0.5:50:nRois); yticklabels([1,50:50:nRois]);
hold on;
xline(1,'--w','LineWidth',1); xline(3,'--w','LineWidth',1); xline(5,'--w','LineWidth',1)
xline(2,'-w','LineWidth',2); xline(4,'-w','LineWidth',2);
title('Trial onset');
ylabel([fovname,' PC dendrite number']);
subplot(6,15,[5 20 35]);
imagesc(colsOn); axis image; axis off;


subplot(6,15,[1:4,16:19,31:34]+5);
imagesc(t2p,[1:nRois],v145movPlot.spkMap',[0 10]);
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks(0:.5:6); xticklabels({[],[],'Naive',[],[],[],'Trained',[],[],[],'Adaptation',[],[]}); xtickangle(0);
yticks(0.5:50:nRois); yticklabels([1,50:50:nRois]);
hold on;
xline(1,'--w','LineWidth',1); xline(3,'--w','LineWidth',1); xline(5,'--w','LineWidth',1)
xline(2,'-w','LineWidth',2); xline(4,'-w','LineWidth',2);
title('Movement onset');
subplot(6,15,[5 20 35]+5);
imagesc(colsMov); axis image; axis off;


subplot(6,15,[1:4,16:19,31:34]+10);
imagesc(t2p,[1:nRois],v145rewTPlot.spkMap',[0 10]);
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks(0:.5:6); xticklabels({[],[],'Naive',[],[],[],'Trained',[],[],[],'Adaptation',[],[]}); xtickangle(0);
yticks(0.5:50:nRois); yticklabels([1,50:50:nRois]);
hold on;
xline(1,'--w','LineWidth',1); xline(3,'--w','LineWidth',1); xline(5,'--w','LineWidth',1)
xline(2,'-w','LineWidth',3); xline(4,'-w','LineWidth',3);
title('Trial reward');
subplot(6,15,[5 20 35]+10);
imagesc(colsRewT); axis image; axis off;

subplot(6,15,[46:49]);
h = boundedline(t2p,mean(v145onPlot.spkMu,2), std(v145onPlot.spkMu,[],2),'-k');
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
h = boundedline(t2p,mean(v145movPlot.spkMu,2), std(v145movPlot.spkMu,[],2),'-k');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(5,'--b','LineWidth',1)
xline(2,'-b','LineWidth',2); xline(4,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 6 0 12]);
xticks(0:.5:6); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:3:12]);

subplot(6,15,[46:49]+10);
h = boundedline(t2p,mean(v145rewTPlot.spkMu,2), std(v145rewTPlot.spkMu,[],2),'-k');
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
h = boundedline(tvr,mean(v145onPlot.lickMu,2), std(v145onPlot.lickMu,[],2),'-k');
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
h = boundedline(tvr,mean(v145movPlot.lickMu,2), std(v145movPlot.lickMu,[],2),'-k');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(5,'--b','LineWidth',1)
xline(2,'-b','LineWidth',2); xline(4,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 6 0 1]);
xticks(0:.5:6); xticklabels({});  % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:.25:1]);

subplot(6,15,[61:64]+10);
h = boundedline(tvr,mean(v145rewTPlot.lickMu,2), std(v145rewTPlot.lickMu,[],2),'-k');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(5,'--b','LineWidth',1)
xline(2,'-b','LineWidth',2); xline(4,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 6 0 1]);
xticks(0:.5:6); xticklabels({});  % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:.25:1]);

subplot(6,15,[76:79]);
h = boundedline(tvr,mean(v145onPlot.movMu,2), std(v145onPlot.movMu,[],2),'-k');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(5,'--b','LineWidth',1)
xline(2,'-b','LineWidth',2); xline(4,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 6 -10 10]);
xticks(0:.5:6); xticklabels({});  % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([-10:5:10]);
ylabel('Movement (cm/s)');

subplot(6,15,[76:79]+5);
h = boundedline(tvr,mean(v145movPlot.movMu,2), std(v145movPlot.movMu,[],2),'-k');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(5,'--b','LineWidth',1)
xline(2,'-b','LineWidth',2); xline(4,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 6 -10 10]);
xticks(0:.5:6); xticklabels(repmat([-1:.5:0.5],1,3));
xlabel('Time (s)');
yticks([-10:5:10]);

subplot(6,15,[76:79]+10);
h = boundedline(tvr,mean(v145rewTPlot.movMu,2), std(v145rewTPlot.movMu,[],2),'-k');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(5,'--b','LineWidth',1)
xline(2,'-b','LineWidth',2); xline(4,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 6 -10 10]);
xticks(0:.5:6); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([-10:5:10]);

%% Plot heatmaps and mean responses for rewards

figout(2) = figure;
subplot(6,15,[1:4,16:19,31:34]);
imagesc(t2p,[1:nRois],v145rewRPlot.spkMap',[0 10]);
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks(0:.5:6); xticklabels({[],[],'Naive',[],[],[],'Trained',[],[],[],'Adaptation',[],[]}); xtickangle(0);
yticks(0.5:50:nRois); yticklabels([1,50:50:nRois]);
hold on;
xline(1,'--w','LineWidth',1); xline(3,'--w','LineWidth',1); xline(5,'--w','LineWidth',1)
xline(2,'-w','LineWidth',3); xline(4,'-w','LineWidth',3);
title('Random reward');
ylabel([fovname,' PC dendrite number']);
subplot(6,15,[5 20 35]);
imagesc(colsRewR); axis image; axis off;

subplot(6,15,[1:4,16:19,31:34]+5);
imagesc(t2p,[1:nRois],v145rewTPlot.spkMap',[0 10]);
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks(0:.5:6); xticklabels({[],[],'Naive',[],[],[],'Trained',[],[],[],'Adaptation',[],[]}); xtickangle(0);
yticks(0.5:50:nRois); yticklabels([1,50:50:nRois]);
hold on;
xline(1,'--w','LineWidth',1); xline(3,'--w','LineWidth',1); xline(5,'--w','LineWidth',1)
xline(2,'-w','LineWidth',3); xline(4,'-w','LineWidth',3);
title('Trial reward');
subplot(6,15,[5 20 35]+5);
imagesc(colsRewT); axis image; axis off;

subplot(6,15,[1:4,16:19,31:34]+10);
imagesc(t2p,[1:nRois],v145rewCPlot.spkMap',[0 10]);
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks(0:.5:6); xticklabels({[],[],'Naive',[],[],[],'Trained',[],[],[],'Adaptation',[],[]}); xtickangle(0);
yticks(0.5:50:nRois); yticklabels([1,50:50:nRois]);
hold on;
xline(1,'--w','LineWidth',1); xline(3,'--w','LineWidth',1); xline(5,'--w','LineWidth',1)
xline(2,'-w','LineWidth',3); xline(4,'-w','LineWidth',3);
title('Cued reward');
subplot(6,15,[5 20 35]+10);
imagesc(colsRewC); axis image; axis off;

subplot(6,15,[46:49]);
h = boundedline(t2p,mean(v145rewRPlot.spkMu,2), std(v145rewRPlot.spkMu,[],2),'-k');
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
h = boundedline(t2p,mean(v145rewTPlot.spkMu,2), std(v145rewTPlot.spkMu,[],2),'-k');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(5,'--b','LineWidth',1)
xline(2,'-b','LineWidth',2); xline(4,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 6 0 12]);
xticks(0:.5:6); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:3:12]);

subplot(6,15,[46:49]+10);
h = boundedline(t2p,mean(v145rewCPlot.spkMu,2), std(v145rewCPlot.spkMu,[],2),'-k');
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
h = boundedline(tvr,mean(v145rewRPlot.lickMu,2), std(v145rewRPlot.lickMu,[],2),'-k');
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
h = boundedline(tvr,mean(v145rewTPlot.lickMu,2), std(v145rewTPlot.lickMu,[],2),'-k');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(5,'--b','LineWidth',1)
xline(2,'-b','LineWidth',2); xline(4,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 6 0 1]);
xticks(0:.5:6); xticklabels({});  % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:.25:1]);

subplot(6,15,[61:64]+10);
h = boundedline(tvr,mean(v145rewCPlot.lickMu,2), std(v145rewCPlot.lickMu,[],2),'-k');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(5,'--b','LineWidth',1)
xline(2,'-b','LineWidth',2); xline(4,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 6 0 1]);
xticks(0:.5:6); xticklabels({});  % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:.25:1]);

subplot(6,15,[76:79]);
h = boundedline(tvr,mean(v145rewRPlot.movMu,2), std(v145rewRPlot.movMu,[],2),'-k');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(5,'--b','LineWidth',1)
xline(2,'-b','LineWidth',2); xline(4,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 6 -10 10]);
xticks(0:.5:6); xticklabels({});  % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([-10:5:10]);
ylabel('Movement (cm/s)');

subplot(6,15,[76:79]+5);
h = boundedline(tvr,mean(v145rewTPlot.movMu,2), std(v145rewTPlot.movMu,[],2),'-k');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(5,'--b','LineWidth',1)
xline(2,'-b','LineWidth',2); xline(4,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 6 -10 10]);
xticks(0:.5:6); xticklabels(repmat([-1:.5:0.5],1,3));
xlabel('Time (s)');
yticks([-10:5:10]);

subplot(6,15,[76:79]+10);
h = boundedline(tvr,mean(v145rewCPlot.movMu,2), std(v145rewCPlot.movMu,[],2),'-k');
set(h,'linewidth',1.5);
hold on;
xline(1,'--b','LineWidth',1); xline(3,'--b','LineWidth',1); xline(5,'--b','LineWidth',1)
xline(2,'-b','LineWidth',2); xline(4,'-b','LineWidth',2);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 6 -10 10]);
xticks(0:.5:6); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([-10:5:10]);

%% Summary plots - behaviour

figout(3) = figure;
subplot(3,3,1); hold on;
h = boundedline(t2p(1:end/3),mean(v145onPlot.spkMu(1:end/3,:),2), std(v145onPlot.spkMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(t2p(1:end/3),mean(v145onPlot.spkMu(end/3*2+1:end,:),2), std(v145onPlot.spkMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(t2p(1:end/3),mean(v145onPlot.spkMu(end/3+1:end/3*2,:),2), std(v145onPlot.spkMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 0 12]);
xticks(0:.5:2); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:3:12]);
ylabel('Pop. response (spk/s)');
title('Trial onset');

subplot(3,3,2); hold on;
h = boundedline(t2p(1:end/3),mean(v145movPlot.spkMu(1:end/3,:),2), std(v145movPlot.spkMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(t2p(1:end/3),mean(v145movPlot.spkMu(end/3*2+1:end,:),2), std(v145movPlot.spkMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(t2p(1:end/3),mean(v145movPlot.spkMu(end/3+1:end/3*2,:),2), std(v145movPlot.spkMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 0 12]);
xticks(0:.5:2); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:3:12]);
title('Movement onset');

subplot(3,3,3); hold on;
h = boundedline(t2p(1:end/3),mean(v145rewTPlot.spkMu(1:end/3,:),2), std(v145rewTPlot.spkMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(t2p(1:end/3),mean(v145rewTPlot.spkMu(end/3*2+1:end,:),2), std(v145rewTPlot.spkMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(t2p(1:end/3),mean(v145rewTPlot.spkMu(end/3+1:end/3*2,:),2), std(v145rewTPlot.spkMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 0 12]);
xticks(0:.5:2); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:3:12]);
title('Trial reward');

subplot(3,3,4); hold on;
h = boundedline(tvr(1:end/3),mean(v145onPlot.lickMu(1:end/3,:),2), std(v145onPlot.lickMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v145onPlot.lickMu(end/3*2+1:end,:),2), std(v145onPlot.lickMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v145onPlot.lickMu(end/3+1:end/3*2,:),2), std(v145onPlot.lickMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 0 1]);
xticks(0:.5:2); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:.25:1]);
ylabel('Licking');

subplot(3,3,5); hold on;
h = boundedline(tvr(1:end/3),mean(v145movPlot.lickMu(1:end/3,:),2), std(v145movPlot.lickMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v145movPlot.lickMu(end/3*2+1:end,:),2), std(v145movPlot.lickMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v145movPlot.lickMu(end/3+1:end/3*2,:),2), std(v145movPlot.lickMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 0 1]);
xticks(0:.5:2); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:.25:1]);

subplot(3,3,6); hold on;
h = boundedline(tvr(1:end/3),mean(v145rewTPlot.lickMu(1:end/3,:),2), std(v145rewTPlot.lickMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v145rewTPlot.lickMu(end/3*2+1:end,:),2), std(v145rewTPlot.lickMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v145rewTPlot.lickMu(end/3+1:end/3*2,:),2), std(v145rewTPlot.lickMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 0 1]);
xticks(0:.5:2); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:.25:1]);

subplot(3,3,7); hold on;
h = boundedline(tvr(1:end/3),mean(v145onPlot.movMu(1:end/3,:),2), std(v145onPlot.movMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v145onPlot.movMu(end/3*2+1:end,:),2), std(v145onPlot.movMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v145onPlot.movMu(end/3+1:end/3*2,:),2), std(v145onPlot.movMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 -10 10]);
xticks(0:.5:2); xticklabels([-1:.5:1]);
yticks([-10:5:10]);
ylabel('Movement (cm/s)');

subplot(3,3,8); hold on;
h = boundedline(tvr(1:end/3),mean(v145movPlot.movMu(1:end/3,:),2), std(v145movPlot.movMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v145movPlot.movMu(end/3*2+1:end,:),2), std(v145movPlot.movMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v145movPlot.movMu(end/3+1:end/3*2,:),2), std(v145movPlot.movMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 -10 10]);
xticks(0:.5:2); xticklabels([-1:.5:1]);
yticks([-10:5:10]);
xlabel('Time (s)');

subplot(3,3,9); hold on;
h = boundedline(tvr(1:end/3),mean(v145rewTPlot.movMu(1:end/3,:),2), std(v145rewTPlot.movMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v145rewTPlot.movMu(end/3*2+1:end,:),2), std(v145rewTPlot.movMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v145rewTPlot.movMu(end/3+1:end/3*2,:),2), std(v145rewTPlot.movMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 -10 10]);
xticks(0:.5:2); xticklabels([-1:.5:1]);
yticks([-10:5:10]);

%% Summary plots - reward

figout(4) = figure;
subplot(3,3,1); hold on;
h = boundedline(t2p(1:end/3),mean(v145rewRPlot.spkMu(1:end/3,:),2), std(v145rewRPlot.spkMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(t2p(1:end/3),mean(v145rewRPlot.spkMu(end/3*2+1:end,:),2), std(v145rewRPlot.spkMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(t2p(1:end/3),mean(v145rewRPlot.spkMu(end/3+1:end/3*2,:),2), std(v145rewRPlot.spkMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 0 15]);
xticks(0:.5:2); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:3:15]);
ylabel('Pop. response (spk/s)');
title('Random reward');

subplot(3,3,2); hold on;
h = boundedline(t2p(1:end/3),mean(v145rewTPlot.spkMu(1:end/3,:),2), std(v145rewTPlot.spkMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(t2p(1:end/3),mean(v145rewTPlot.spkMu(end/3*2+1:end,:),2), std(v145rewTPlot.spkMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(t2p(1:end/3),mean(v145rewTPlot.spkMu(end/3+1:end/3*2,:),2), std(v145rewTPlot.spkMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 0 15]);
xticks(0:.5:2); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:3:15]);
title('Trial reward');

subplot(3,3,3); hold on;
h = boundedline(t2p(1:end/3),mean(v145rewCPlot.spkMu(1:end/3,:),2), std(v145rewCPlot.spkMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(t2p(1:end/3),mean(v145rewCPlot.spkMu(end/3*2+1:end,:),2), std(v145rewCPlot.spkMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(t2p(1:end/3),mean(v145rewCPlot.spkMu(end/3+1:end/3*2,:),2), std(v145rewCPlot.spkMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 0 15]);
xticks(0:.5:2); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:3:15]);
title('Cued reward');

subplot(3,3,4); hold on;
h = boundedline(tvr(1:end/3),mean(v145rewRPlot.lickMu(1:end/3,:),2), std(v145rewRPlot.lickMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v145rewRPlot.lickMu(end/3*2+1:end,:),2), std(v145rewRPlot.lickMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v145rewRPlot.lickMu(end/3+1:end/3*2,:),2), std(v145rewRPlot.lickMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 0 1]);
xticks(0:.5:2); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:.25:1]);
ylabel('Licking');

subplot(3,3,5); hold on;
h = boundedline(tvr(1:end/3),mean(v145rewTPlot.lickMu(1:end/3,:),2), std(v145rewTPlot.lickMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v145rewTPlot.lickMu(end/3*2+1:end,:),2), std(v145rewTPlot.lickMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v145rewTPlot.lickMu(end/3+1:end/3*2,:),2), std(v145rewTPlot.lickMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 0 1]);
xticks(0:.5:2); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:.25:1]);

subplot(3,3,6); hold on;
h = boundedline(tvr(1:end/3),mean(v145rewCPlot.lickMu(1:end/3,:),2), std(v145rewCPlot.lickMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v145rewCPlot.lickMu(end/3*2+1:end,:),2), std(v145rewCPlot.lickMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v145rewCPlot.lickMu(end/3+1:end/3*2,:),2), std(v145rewCPlot.lickMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 0 1]);
xticks(0:.5:2); xticklabels({}); % xticklabels(repmat([-1:.5:0.5],1,3));
yticks([0:.25:1]);

subplot(3,3,7); hold on;
h = boundedline(tvr(1:end/3),mean(v145rewRPlot.movMu(1:end/3,:),2), std(v145rewRPlot.movMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v145rewRPlot.movMu(end/3*2+1:end,:),2), std(v145rewRPlot.movMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v145rewRPlot.movMu(end/3+1:end/3*2,:),2), std(v145rewRPlot.movMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 -10 10]);
xticks(0:.5:2); xticklabels([-1:.5:1]);
yticks([-10:5:10]);
ylabel('Movement (cm/s)');

subplot(3,3,8); hold on;
h = boundedline(tvr(1:end/3),mean(v145rewTPlot.movMu(1:end/3,:),2), std(v145rewTPlot.movMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v145rewTPlot.movMu(end/3*2+1:end,:),2), std(v145rewTPlot.movMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v145rewTPlot.movMu(end/3+1:end/3*2,:),2), std(v145rewTPlot.movMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 -10 10]);
xticks(0:.5:2); xticklabels([-1:.5:1]);
yticks([-10:5:10]);
xlabel('Time (s)');

subplot(3,3,9); hold on;
h = boundedline(tvr(1:end/3),mean(v145rewCPlot.movMu(1:end/3,:),2), std(v145rewCPlot.movMu(1:end/3,:),[],2),'-b','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v145rewCPlot.movMu(end/3*2+1:end,:),2), std(v145rewCPlot.movMu(end/3*2+1:end,:),[],2),'-r','alpha'); set(h,'linewidth',1.5);
h = boundedline(tvr(1:end/3),mean(v145rewCPlot.movMu(end/3+1:end/3*2,:),2), std(v145rewCPlot.movMu(end/3+1:end/3*2,:),[],2),'-k','alpha'); set(h,'linewidth',1.5);
xline(1,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
axis([0 2 -10 10]);
xticks(0:.5:2); xticklabels([-1:.5:1]);
yticks([-10:5:10]);

end