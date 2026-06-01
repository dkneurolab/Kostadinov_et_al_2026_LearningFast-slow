function fig4_matchPercSig(trainSesh,testSesh,trialType,paperFold,savebool)

figFold = fullfile(paperFold,'Figures','Fig4','Sub2_CrossDayPredictions');

GLMheatmap = struct;

%% Loop through train and test sessions and populate heatmap
a = 0;
for iTrain = 1:numel(trainSesh)
    for iTest = 1:numel(testSesh)
        if ~strcmpi(trainSesh{iTrain},testSesh{iTest})
            a = a + 1;
            GLMheatmap(a).heatPos = [iTrain iTest];
            GLMout = load(fullfile(figFold,sprintf('GLM_%s_%s_%s.mat',testSesh{iTest},trainSesh{iTrain},trialType))); GLMout = GLMout.GLMout;
            GLMlobv = GLMout(cat(1,GLMout.ilobv)); GLMsim2 = GLMout(cat(1,GLMout.isim2)); clear GLMout;
            % Process lobule V data
            onHeat = zeros(size(GLMlobv))'; movHeat = onHeat; rewHeat = onHeat; lickHeat = onHeat;
            for iFOVs = 1:numel(GLMlobv)
                glmStr = GLMlobv(iFOVs).GLMstruct; glmStr = glmStr(logical(cat(1,glmStr.sigV4)));
                boolVar = cat(1,glmStr.sigBool_shufVis);
                onHeat(iFOVs) = sum(boolVar,'omitnan')/sum(~isnan(boolVar))*100;
                boolVar = cat(1,glmStr.sigBool_shufMov);
                movHeat(iFOVs) = sum(boolVar,'omitnan')/sum(~isnan(boolVar))*100;
                boolVar = cat(1,glmStr.sigBool_shufRew);
                rewHeat(iFOVs) = sum(boolVar,'omitnan')/sum(~isnan(boolVar))*100;
                boolVar = cat(1,glmStr.sigBool_shufLick);
                lickHeat(iFOVs) = sum(boolVar,'omitnan')/sum(~isnan(boolVar))*100; 
            end
            GLMheatmap(a).onLobv = mean(onHeat);
            GLMheatmap(a).movLobv = mean(movHeat);
            GLMheatmap(a).rewLobv = mean(rewHeat);
            GLMheatmap(a).lickLobv = mean(lickHeat);
            
            % Process lobule V data
            onHeat = zeros(size(GLMsim2))'; movHeat = onHeat; rewHeat = onHeat; lickHeat = onHeat;
            for iFOVs = 1:numel(GLMsim2)
                glmStr = GLMsim2(iFOVs).GLMstruct; glmStr = glmStr(logical(cat(1,glmStr.sigV4)));
                boolVar = cat(1,glmStr.sigBool_shufVis);
                onHeat(iFOVs) = sum(boolVar,'omitnan')/sum(~isnan(boolVar))*100;
                boolVar = cat(1,glmStr.sigBool_shufMov);
                movHeat(iFOVs) = sum(boolVar,'omitnan')/sum(~isnan(boolVar))*100;
                boolVar = cat(1,glmStr.sigBool_shufRew);
                rewHeat(iFOVs) = sum(boolVar,'omitnan')/sum(~isnan(boolVar))*100;
                boolVar = cat(1,glmStr.sigBool_shufLick);
                lickHeat(iFOVs) = sum(boolVar,'omitnan')/sum(~isnan(boolVar))*100; 
            end
            GLMheatmap(a).onSim2 = mean(onHeat);
            GLMheatmap(a).movSim2 = mean(movHeat);
            GLMheatmap(a).rewSim2 = mean(rewHeat);
            GLMheatmap(a).lickSim2 = mean(lickHeat);
            
            % Keep all datapoints JIC
            GLMheatmap(a).onLobvAll = onHeat;
            GLMheatmap(a).movLobvAll = movHeat;
            GLMheatmap(a).rewLobvAll = rewHeat;
            GLMheatmap(a).lickLobvAll = lickHeat;
            
            GLMheatmap(a).onSim2All = onHeat;
            GLMheatmap(a).movSim2All = movHeat;
            GLMheatmap(a).rewSim2All = rewHeat;
            GLMheatmap(a).lickSim2All = lickHeat;
            
        end
    end
end

clear onHeat movHeat rewHeat lickHeat iFOVs glmStr boolVar a

%% Heatmap

iPos = cat(1,GLMheatmap.heatPos);
iPos = sub2ind([iTrain iTest],iPos(:,1),iPos(:,2));
cMax = 50;

fHeat = figure;
% Lobule V - onset
heatMat = zeros(iTrain,iTest)*nan;
heatMat(iPos) = cat(1,GLMheatmap.onLobv);
subplot(2,4,1);
imagesc(heatMat,[0 cMax]);
xticks(1:iTest); yticks(1:iTrain);
xticklabels(testSesh); yticklabels(trainSesh);
xlabel('Testing sessions');
ylabel('Training sessions');
title('Lobule V: object');
colormap viridis; box off; axis image;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');

% Lobule V - movement
heatMat = zeros(iTrain,iTest)*nan;
heatMat(iPos) = cat(1,GLMheatmap.movLobv);
subplot(2,4,2);
imagesc(heatMat,[0 cMax]);
xticks(1:iTest); yticks(1:iTrain);
xticklabels(testSesh); yticklabels(trainSesh);
title('Lobule V: movement');
colormap viridis; box off; axis image;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');

% Lobule V - reward
heatMat = zeros(iTrain,iTest)*nan;
heatMat(iPos) = cat(1,GLMheatmap.rewLobv);
subplot(2,4,3);
imagesc(heatMat,[0 cMax]);
xticks(1:iTest); yticks(1:iTrain);
xticklabels(testSesh); yticklabels(trainSesh);
title('Lobule V: reward');
colormap viridis; box off; axis image;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');

% Lobule V - licking
heatMat = zeros(iTrain,iTest)*nan;
heatMat(iPos) = cat(1,GLMheatmap.lickLobv);
subplot(2,4,4);
imagesc(heatMat,[0 cMax]);
xticks(1:iTest); yticks(1:iTrain);
xticklabels(testSesh); yticklabels(trainSesh);
title('Lobule V: licking');
colormap viridis; box off; axis image;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');

% Simplex - onset
heatMat = zeros(iTrain,iTest)*nan;
heatMat(iPos) = cat(1,GLMheatmap.onSim2);
subplot(2,4,5);
imagesc(heatMat,[0 cMax]);
xticks(1:iTest); yticks(1:iTrain);
xticklabels(testSesh); yticklabels(trainSesh);
xlabel('Testing sessions');
ylabel('Training sessions');
title('Simplex: object');
colormap viridis; box off; axis image;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');

% Simplex - movement
heatMat = zeros(iTrain,iTest)*nan;
heatMat(iPos) = cat(1,GLMheatmap.movSim2);
subplot(2,4,6);
imagesc(heatMat,[0 cMax]);
xticks(1:iTest); yticks(1:iTrain);
xticklabels(testSesh); yticklabels(trainSesh);
title('Simplex: movement');
colormap viridis; box off; axis image;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');

% Simplex - reward
heatMat = zeros(iTrain,iTest)*nan;
heatMat(iPos) = cat(1,GLMheatmap.rewSim2);
subplot(2,4,7);
imagesc(heatMat,[0 cMax]);
xticks(1:iTest); yticks(1:iTrain);
xticklabels(testSesh); yticklabels(trainSesh);
title('Simplex: reward');
colormap viridis; box off; axis image;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');

% Simplex - licking
heatMat = zeros(iTrain,iTest)*nan;
heatMat(iPos) = cat(1,GLMheatmap.lickSim2);
subplot(2,4,8);
imagesc(heatMat,[0 cMax]);
xticks(1:iTest); yticks(1:iTrain);
xticklabels(testSesh); yticklabels(trainSesh);
title('Simplex: licking');
colormap viridis; box off; axis image;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');

if savebool
    savename = 'LobVvsSim2_matchPercHeat';
    savefig(fHeat,fullfile(figFold,[savename,'.fig']));
    saveas(fHeat,fullfile(figFold,[savename,'.png']));
    print(fHeat,fullfile(figFold,[savename,'.eps']), '-depsc', '-painters');
    save(fullfile(figFold,[savename,'.mat']),'GLMheatmap');
    
end

