function fig5_adaptRatiovsDist(seshv4, paperFold, figtype, savebool)
%% Default parameters:
imfs = 30;
fovum = 670;
fovpx = 512;

%% Get adaptation data and format for convenience
% Load adaptSummary and concatenate:
adaptFold = fullfile(paperFold,'Figures','Fig5','Sub2_behavSummary');
adaptSummary = load(fullfile(adaptFold,'v5_adaptSummary.mat'));
[adaptSummary.lobvdata.fov] = deal('lobv');
[adaptSummary.sim2data.fov] = deal('sim2');
adaptSummary = cat(2,adaptSummary.lobvdata,adaptSummary.sim2data);

% Get current field order and make 'fov' the second one (after 'name')
flds = fieldnames(adaptSummary);
newOrder = [flds(1); 'fov'; flds(2:end-1)];
adaptSummary = orderfields(adaptSummary, newOrder);
keepFields = fieldnames(adaptSummary);
keepFields = keepFields(1:6);
adaptDistStruct = rmfield(adaptSummary, setdiff(fieldnames(adaptSummary), keepFields));


%% Loop through and get matching xCoordPx and compute adaptation ratios
for isesh = 1:numel(adaptSummary)
    % Get xCoordPx for all matched ROIs
    xCtrPx = [];
    for jseshv4 = 1:numel(seshv4)
        if (strcmpi(adaptSummary(isesh).name,seshv4(jseshv4).name) && strcmpi(adaptSummary(isesh).fov, seshv4(jseshv4).fov))
            xCtrPx = seshv4(jseshv4).xCtrsPx(adaptSummary(isesh).iROIv4);
        end
    end
    adaptDistStruct(isesh).xCtrPx = xCtrPx;
        
    % Deal with movement first:
    adaptDistStruct(isesh).xCtrMovUm = xCtrPx(adaptDistStruct(isesh).sigBoolMov)*fovum/fovpx;
    adaptDistStruct(isesh).spkMovBase = mean(adaptSummary(isesh).spkMovBase(end/2-.3*imfs+1:end/2,:))';
    adaptDistStruct(isesh).spkMovMid = mean(adaptSummary(isesh).spkMovMid(end/2-.3*imfs+1:end/2,:))';
    adaptDistStruct(isesh).spkMovRatio = adaptDistStruct(isesh).spkMovMid./adaptDistStruct(isesh).spkMovBase;
    adaptDistStruct(isesh).spkMovDirection = adaptDistStruct(isesh).spkMovRatio > 1;
    iGoodMov = (adaptDistStruct(isesh).spkMovRatio ~= 0) & isfinite(adaptDistStruct(isesh).spkMovRatio) & ~isnan(adaptDistStruct(isesh).spkMovRatio);
    adaptDistStruct(isesh).iGoodMov = iGoodMov;
    matchDist0 = adaptDistStruct(isesh).xCtrMovUm(iGoodMov) - adaptDistStruct(isesh).xCtrMovUm(iGoodMov).';    
    adaptDistStruct(isesh).movMatch_Dist = abs(matchDist0(triu(true(size(matchDist0)),1)));
    
    % Deal with reward:
    adaptDistStruct(isesh).xCtrRewUm = xCtrPx(adaptDistStruct(isesh).sigBoolRew)*fovum/fovpx;
    adaptDistStruct(isesh).spkRewBase = mean(adaptSummary(isesh).spkRewBase(end/2-.3*imfs+1:end/2,:))';
    adaptDistStruct(isesh).spkRewMid = mean(adaptSummary(isesh).spkRewMid(end/2-.3*imfs+1:end/2,:))';
    adaptDistStruct(isesh).spkRewRatio = adaptDistStruct(isesh).spkRewMid./adaptDistStruct(isesh).spkRewBase;
    adaptDistStruct(isesh).spkRewDirection = adaptDistStruct(isesh).spkRewRatio > 1;    
    iGoodRew = (adaptDistStruct(isesh).spkRewRatio ~= 0) & isfinite(adaptDistStruct(isesh).spkRewRatio) & ~isnan(adaptDistStruct(isesh).spkRewRatio);
    adaptDistStruct(isesh).iGoodRew = iGoodRew;
    matchDist0 = adaptDistStruct(isesh).xCtrRewUm(iGoodRew) - adaptDistStruct(isesh).xCtrRewUm(iGoodRew).';
    adaptDistStruct(isesh).rewMatch_Dist = abs(matchDist0(triu(true(size(matchDist0)),1)));



    % Common switch case for type of processing we will do:
    switch figtype
        case 'meandiff'
            % Movement
            matchSimilarity0 = adaptDistStruct(isesh).spkMovRatio(iGoodMov) - adaptDistStruct(isesh).spkMovRatio(iGoodMov).';
            adaptDistStruct(isesh).movMatch_Similarity = abs(matchSimilarity0(triu(true(size(matchSimilarity0)),1)));
            % Reward
            matchSimilarity0 = adaptDistStruct(isesh).spkRewRatio(iGoodRew) - adaptDistStruct(isesh).spkRewRatio(iGoodRew).';
            adaptDistStruct(isesh).rewMatch_Similarity = abs(matchSimilarity0(triu(true(size(matchSimilarity0)),1)));            
            
            ylimPlot = 3;

        case 'bothup'
            % Movement
            matchSimilarity0 = adaptDistStruct(isesh).spkMovDirection(iGoodMov) + adaptDistStruct(isesh).spkMovDirection(iGoodMov).';
            movMatch_Similarity = matchSimilarity0(triu(true(size(matchSimilarity0)),1));
            adaptDistStruct(isesh).movMatch_Similarity = movMatch_Similarity > 1;            
            % Reward
            matchSimilarity0 = adaptDistStruct(isesh).spkRewDirection(iGoodRew) + adaptDistStruct(isesh).spkRewDirection(iGoodRew).';
            rewMatch_Similarity = matchSimilarity0(triu(true(size(matchSimilarity0)),1));
            adaptDistStruct(isesh).rewMatch_Similarity = rewMatch_Similarity > 1;            

            ylimPlot = .5;

        case 'samedir'
            % Movement
            matchSimilarity0 = adaptDistStruct(isesh).spkMovDirection(iGoodMov) - adaptDistStruct(isesh).spkMovDirection(iGoodMov).';
            adaptDistStruct(isesh).movMatch_Similarity = ~(matchSimilarity0(triu(true(size(matchSimilarity0)),1)) == 1);            
            % Reward
            matchSimilarity0 = adaptDistStruct(isesh).spkRewDirection(iGoodRew) - adaptDistStruct(isesh).spkRewDirection(iGoodRew).';
            adaptDistStruct(isesh).rewMatch_Similarity = ~(matchSimilarity0(triu(true(size(matchSimilarity0)),1)) == 1);            

            ylimPlot = 1;

        case 'abslog'
            % Movement
            matchSimilarity0 = log2(adaptDistStruct(isesh).spkMovRatio(iGoodMov) ./ adaptDistStruct(isesh).spkMovRatio(iGoodMov).');
            adaptDistStruct(isesh).movMatch_Similarity = abs(matchSimilarity0(triu(true(size(matchSimilarity0)),1)));
            % Reward
            matchSimilarity0 = log2(adaptDistStruct(isesh).spkRewRatio(iGoodRew) ./ adaptDistStruct(isesh).spkRewRatio(iGoodRew).');
            adaptDistStruct(isesh).rewMatch_Similarity = abs(matchSimilarity0(triu(true(size(matchSimilarity0)),1)));            

            ylimPlot = 2;

        case 'absdiv'
            % Movement
            matchSimilarity0 = log2(adaptDistStruct(isesh).spkMovRatio(iGoodMov) ./ adaptDistStruct(isesh).spkMovRatio(iGoodMov).');
            adaptDistStruct(isesh).movMatch_Similarity = 2.^abs(matchSimilarity0(triu(true(size(matchSimilarity0)),1)));
            % Reward
            matchSimilarity0 = log2(adaptDistStruct(isesh).spkRewRatio(iGoodRew) ./ adaptDistStruct(isesh).spkRewRatio(iGoodRew).');
            adaptDistStruct(isesh).rewMatch_Similarity = 2.^abs(matchSimilarity0(triu(true(size(matchSimilarity0)),1)));            

            ylimPlot = 5;


    end
    % fprintf('%i. movsim: %f\n',isesh, mean(adaptDistStruct(isesh).movMatch_Similarity));
    % fprintf('%i. rewsim: %f\n',isesh, mean(adaptDistStruct(isesh).rewMatch_Similarity));
    
end

%% Concatenate and do stats:
movDistAll = cat(1,adaptDistStruct.movMatch_Dist);
movSimAll = cat(1,adaptDistStruct.movMatch_Similarity);
rewDistAll = cat(1,adaptDistStruct.rewMatch_Dist);
rewSimAll = cat(1,adaptDistStruct.rewMatch_Similarity);

f1 = figure;
subplot(2,2,1);
scatter(movDistAll,movSimAll, 'filled', MarkerEdgeColor='none',MarkerFaceColor = 'b', MarkerFaceAlpha = 0.1); axis square; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlim([-25 650]); ylim([-0.05 1.05]*ylimPlot*2)
ylabel('Difference in adaptation ratio')
xlabel('Intercellular distance (um)')
title('Fast learning: movement')
subplot(2,2,2);
scatter(rewDistAll,rewSimAll, 'filled', MarkerEdgeColor='none',MarkerFaceColor = 'b', MarkerFaceAlpha = 0.1); axis square; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlim([-25 650]); ylim([-0.05 1.05]*ylimPlot*2)
ylabel('Difference in adaptation ratio')
xlabel('Intercellular distance (um)')
title('Fast learning: reward')

iMov_under100 = find(movDistAll < 100);
iMov_over100 = find(movDistAll >= 100);
iRew_under100 = find(rewDistAll < 100);
iRew_over100 = find(rewDistAll >= 100);

subplot(2,2,3); hold on;
xplot = randn(numel(iMov_under100),1)*.2;
xplot(abs(xplot) > 0.5) = 0;
scatter(xplot+1,movSimAll(iMov_under100),'filled', MarkerEdgeColor='none',MarkerFaceColor = 'k', MarkerFaceAlpha = 0.1); axis square; box off;
xplot = randn(numel(iMov_over100),1)*.2;
xplot(abs(xplot) > 0.5) = 0;
scatter(xplot+3,movSimAll(iMov_over100),'filled', MarkerEdgeColor='none',MarkerFaceColor = [0.5 .1 .1], MarkerFaceAlpha = 0.1); axis square; box off;
xlim([0 4]); ylim([-0.05 1.05]*ylimPlot*2)
axis square; box off
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([1,3]); xticklabels({'<100','>=100'})
xlabel('Intercellular distance (um)')
ylabel('Difference in adaptation ratio')
a = ranksum(movSimAll(iMov_under100),movSimAll(iMov_over100));
text(2,1.5,sprintf('p = %.3f',a));

subplot(2,2,4); hold on;
xplot = randn(numel(iRew_under100),1)*.2;
xplot(abs(xplot) > 0.5) = 0;
scatter(xplot+1,rewSimAll(iRew_under100),'filled', MarkerEdgeColor='none',MarkerFaceColor = 'k', MarkerFaceAlpha = 0.3); axis square; box off;
xplot = randn(numel(iRew_over100),1)*.2;
xplot(abs(xplot) > 0.5) = 0;
scatter(xplot+3,rewSimAll(iRew_over100),'filled', MarkerEdgeColor='none',MarkerFaceColor = [0.5 .1 .1], MarkerFaceAlpha = 0.3); axis square; box off;
xlim([0 4]); ylim([-0.05 1.05]*ylimPlot*2)
axis square; box off
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([1,3]); xticklabels({'<100','>=100'})
xlabel('Intercellular distance (um)')
ylabel('Difference in adaptation ratio')
a = ranksum(rewSimAll(iRew_under100),rewSimAll(iRew_over100));
text(2,1.5,sprintf('p = %.3f',a));

suptitle('Distance-dependence of adaptation ratios');

%% Concatenate and do stats: lobule V
ilobv = find(strcmpi({adaptDistStruct.fov}, 'lobv'));
movDistAll_lobv = cat(1,adaptDistStruct(ilobv).movMatch_Dist);
movSimAll_lobv = cat(1,adaptDistStruct(ilobv).movMatch_Similarity);
rewDistAll_lobv = cat(1,adaptDistStruct(ilobv).rewMatch_Dist);
rewSimAll_lobv = cat(1,adaptDistStruct(ilobv).rewMatch_Similarity);

f2 = figure;
subplot(2,2,1);
scatter(movDistAll_lobv,movSimAll_lobv, 'filled', MarkerEdgeColor='none',MarkerFaceColor = 'b', MarkerFaceAlpha = 0.3); axis square; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlim([-25 650]); ylim([-0.05 1.05]*ylimPlot*2)
ylabel('Difference in adaptation ratio')
xlabel('Intercellular distance (um)')
title('Fast learning: movement (lobV)')
subplot(2,2,2);
scatter(rewDistAll_lobv,rewSimAll_lobv, 'filled', MarkerEdgeColor='none',MarkerFaceColor = 'b', MarkerFaceAlpha = 0.3); axis square; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlim([-25 650]); ylim([-0.05 1.05]*ylimPlot*2)
ylabel('Difference in adaptation ratio')
xlabel('Intercellular distance (um)')
title('Fast learning: reward (lobV)')

iMov_under100_lobv = find(movDistAll_lobv < 100);
iMov_over100_lobv = find(movDistAll_lobv >= 100);
iRew_under100_lobv = find(rewDistAll_lobv < 100);
iRew_over100_lobv = find(rewDistAll_lobv >= 100);

subplot(2,2,3); hold on;
xplot = randn(numel(iMov_under100_lobv),1)*.2;
xplot(abs(xplot) > 0.5) = 0;
scatter(xplot+1,movSimAll_lobv(iMov_under100_lobv),'filled', MarkerEdgeColor='none',MarkerFaceColor = 'k', MarkerFaceAlpha = 0.3); axis square; box off;
xplot = randn(numel(iMov_over100_lobv),1)*.2;
xplot(abs(xplot) > 0.5) = 0;
scatter(xplot+3,movSimAll_lobv(iMov_over100_lobv),'filled', MarkerEdgeColor='none',MarkerFaceColor = [0.5 .1 .1], MarkerFaceAlpha = 0.3); axis square; box off;
xlim([0 4]); ylim([-0.05 1.05]*ylimPlot*2)
axis square; box off
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([1,3]); xticklabels({'<100','>=100'})
xlabel('Intercellular distance (um)')
ylabel('Difference in adaptation ratio')
a = ranksum(movSimAll_lobv(iMov_under100_lobv),movSimAll_lobv(iMov_over100_lobv));
text(2,1.5,sprintf('p = %.3f',a));

subplot(2,2,4); hold on;
xplot = randn(numel(iRew_under100_lobv),1)*.2;
xplot(abs(xplot) > 0.5) = 0;
scatter(xplot+1,rewSimAll_lobv(iRew_under100_lobv),'filled', MarkerEdgeColor='none',MarkerFaceColor = 'k', MarkerFaceAlpha = 0.3); axis square; box off;
xplot = randn(numel(iRew_over100_lobv),1)*.2;
xplot(abs(xplot) > 0.5) = 0;
scatter(xplot+3,rewSimAll_lobv(iRew_over100_lobv),'filled', MarkerEdgeColor='none',MarkerFaceColor = [0.5 .1 .1], MarkerFaceAlpha = 0.3); axis square; box off;
xlim([0 4]); ylim([-0.05 1.05]*ylimPlot*2)
axis square; box off
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([1,3]); xticklabels({'<100','>=100'})
xlabel('Intercellular distance (um)')
ylabel('Difference in adaptation ratio')
a = ranksum(rewSimAll_lobv(iRew_under100_lobv),rewSimAll_lobv(iRew_over100_lobv));
text(2,1.5,sprintf('p = %.3f',a));

suptitle('Distance-dependence of adaptation ratios');

%% Concatenate and do stats:
isim2 = find(strcmpi({adaptDistStruct.fov}, 'sim2'));
movDistAll_sim2 = cat(1,adaptDistStruct(isim2).movMatch_Dist);
movSimAll_sim2 = cat(1,adaptDistStruct(isim2).movMatch_Similarity);
rewDistAll_sim2 = cat(1,adaptDistStruct(isim2).rewMatch_Dist);
rewSimAll_sim2 = cat(1,adaptDistStruct(isim2).rewMatch_Similarity);

f3 = figure;
subplot(2,2,1);
scatter(movDistAll_sim2,movSimAll_sim2, 'filled', MarkerEdgeColor='none',MarkerFaceColor = 'b', MarkerFaceAlpha = 0.3); axis square; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlim([-25 650]); ylim([-0.05 1.05]*ylimPlot*2)
ylabel('Difference in adaptation ratio')
xlabel('Intercellular distance (um)')
title('Fast learning: movement (sim2)')
subplot(2,2,2);
scatter(rewDistAll_sim2,rewSimAll_sim2, 'filled', MarkerEdgeColor='none',MarkerFaceColor = 'b', MarkerFaceAlpha = 0.3); axis square; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlim([-25 650]); ylim([-0.05 1.05]*ylimPlot*2)
ylabel('Difference in adaptation ratio')
xlabel('Intercellular distance (um)')
title('Fast learning: reward (sim2)')

iMov_under100_sim2 = find(movDistAll_sim2 < 100);
iMov_over100_sim2 = find(movDistAll_sim2 >= 100);
iRew_under100_sim2 = find(rewDistAll_sim2 < 100);
iRew_over100_sim2 = find(rewDistAll_sim2 >= 100);

subplot(2,2,3); hold on;
xplot = randn(numel(iMov_under100_sim2),1)*.2;
xplot(abs(xplot) > 0.5) = 0;
scatter(xplot+1,movSimAll_sim2(iMov_under100_sim2),'filled', MarkerEdgeColor='none',MarkerFaceColor = 'k', MarkerFaceAlpha = 0.3); axis square; box off;
xplot = randn(numel(iMov_over100_sim2),1)*.2;
xplot(abs(xplot) > 0.5) = 0;
scatter(xplot+3,movSimAll_sim2(iMov_over100_sim2),'filled', MarkerEdgeColor='none',MarkerFaceColor = [0.5 .1 .1], MarkerFaceAlpha = 0.3); axis square; box off;
xlim([0 4]); ylim([-0.05 1.05]*ylimPlot*2)
axis square; box off
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([1,3]); xticklabels({'<100','>=100'})
xlabel('Intercellular distance (um)')
ylabel('Difference in adaptation ratio')
a = ranksum(movSimAll_sim2(iMov_under100_sim2),movSimAll_sim2(iMov_over100_sim2));
text(2,1.5,sprintf('p = %.3f',a));

subplot(2,2,4); hold on;
xplot = randn(numel(iRew_under100_sim2),1)*.2;
xplot(abs(xplot) > 0.5) = 0;
scatter(xplot+1,rewSimAll_sim2(iRew_under100_sim2),'filled', MarkerEdgeColor='none',MarkerFaceColor = 'k', MarkerFaceAlpha = 0.3); axis square; box off;
xplot = randn(numel(iRew_over100_sim2),1)*.2;
xplot(abs(xplot) > 0.5) = 0;
scatter(xplot+3,rewSimAll_sim2(iRew_over100_sim2),'filled', MarkerEdgeColor='none',MarkerFaceColor = [0.5 .1 .1], MarkerFaceAlpha = 0.3); axis square; box off;
xlim([0 4]); ylim([-0.05 1.05]*ylimPlot*2)
axis square; box off
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([1,3]); xticklabels({'<100','>=100'})
xlabel('Intercellular distance (um)')
ylabel('Difference in adaptation ratio')
a = ranksum(rewSimAll_sim2(iRew_under100_sim2),rewSimAll_sim2(iRew_over100_sim2));
text(2,1.5,sprintf('p = %.3f',a));

suptitle('Distance-dependence of adaptation ratios');

%% Summary bar graphs only

f4 = figure;
subplot(2,3,1)
e1 = errorbar([1,3],[mean(movSimAll(iMov_under100)), mean(movSimAll(iMov_over100))], [std(movSimAll(iMov_under100))/sqrt(numel(iMov_under100)), std(movSimAll(iMov_over100))/sqrt(numel(iMov_over100))]);
% e1 = errorbar([1,3],[mean(movSimAll(iMov_under100)), mean(movSimAll(iMov_over100))], [std(movSimAll(iMov_under100)), std(movSimAll(iMov_over100))]);
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
hold on; axis square; box off
scatter([1,3],[mean(movSimAll(iMov_under100)), mean(movSimAll(iMov_over100))],60,'MarkerEdgeColor','none','MarkerFaceColor','k');
xlim([0 4]); ylim([0 1.05]*ylimPlot)
% xlim([0 4]); ylim([-0.05 1.05])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([1,3]); xticklabels({'<100','>=100'})
xlabel('Intercellular distance (um)')
ylabel('Difference in adaptation ratio')
title('Combined FOVs, movement')
a = ranksum(movSimAll(iMov_under100),movSimAll(iMov_over100));
text(2,.5,sprintf('p = %.3f',a));

subplot(2,3,4)
e1 = errorbar([1,3],[mean(rewSimAll(iRew_under100)), mean(rewSimAll(iRew_over100))], [std(rewSimAll(iRew_under100))/sqrt(numel(iRew_under100)), std(rewSimAll(iRew_over100))/sqrt(numel(iRew_over100))]);
% e1 = errorbar([1,3],[mean(rewSimAll(iRew_under100)), mean(rewSimAll(iRew_over100))], [std(rewSimAll(iRew_under100)), std(rewSimAll(iRew_over100))]);
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
hold on; axis square; box off
scatter([1,3],[mean(rewSimAll(iRew_under100)), mean(rewSimAll(iRew_over100))],60,'MarkerEdgeColor','none','MarkerFaceColor','k');
xlim([0 4]); ylim([0 1.05]*ylimPlot)
% xlim([0 4]); ylim([-0.05 1.05])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([1,3]); xticklabels({'<100','>=100'})
xlabel('Intercellular distance (um)')
ylabel('Difference in adaptation ratio')
title('Combined FOVs, reward')
a = ranksum(rewSimAll(iRew_under100),rewSimAll(iRew_over100));
text(2,.5,sprintf('p = %.3f',a));

subplot(2,3,2)
e1 = errorbar([1,3],[mean(movSimAll_lobv(iMov_under100_lobv)), mean(movSimAll_lobv(iMov_over100_lobv))], [std(movSimAll_lobv(iMov_under100_lobv))/sqrt(numel(iMov_under100_lobv)), std(movSimAll_lobv(iMov_over100_lobv))/sqrt(numel(iMov_over100_lobv))]);
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
hold on; axis square; box off
scatter([1,3],[mean(movSimAll_lobv(iMov_under100_lobv)), mean(movSimAll_lobv(iMov_over100_lobv))],60,'MarkerEdgeColor','none','MarkerFaceColor','k');
xlim([0 4]); ylim([0 1.05]*ylimPlot)
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([1,3]); xticklabels({'<100','>=100'})
xlabel('Intercellular distance (um)')
ylabel('Difference in adaptation ratio')
title('Lobule V, movement')
a = ranksum(movSimAll_lobv(iMov_under100_lobv),movSimAll_lobv(iMov_over100_lobv));
text(2,.5,sprintf('p = %.3f',a));

subplot(2,3,5)
e1 = errorbar([1,3],[mean(rewSimAll_lobv(iRew_under100_lobv)), mean(rewSimAll_lobv(iRew_over100_lobv))], [std(rewSimAll_lobv(iRew_under100_lobv))/sqrt(numel(iRew_over100_lobv)), std(rewSimAll_lobv(iRew_over100_lobv))/sqrt(numel(iRew_over100_lobv))]);
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
hold on; axis square; box off
scatter([1,3],[mean(rewSimAll_lobv(iRew_under100_lobv)), mean(rewSimAll_lobv(iRew_over100_lobv))],60,'MarkerEdgeColor','none','MarkerFaceColor','k');
xlim([0 4]); ylim([0 1.05]*ylimPlot)
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([1,3]); xticklabels({'<100','>=100'})
xlabel('Intercellular distance (um)')
ylabel('Difference in adaptation ratio')
title('Lobule V, reward')
a = ranksum(rewSimAll_lobv(iRew_under100_lobv),rewSimAll_lobv(iRew_over100_lobv));
text(2,.5,sprintf('p = %.3f',a));

subplot(2,3,3)
e1 = errorbar([1,3],[mean(movSimAll_sim2(iMov_under100_sim2)), mean(movSimAll_sim2(iMov_over100_sim2))], [std(movSimAll_sim2(iMov_under100_sim2))/sqrt(numel(iMov_under100_sim2)), std(movSimAll_sim2(iMov_over100_sim2))/sqrt(numel(iMov_over100_sim2))]);
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
hold on; axis square; box off
scatter([1,3],[mean(movSimAll_sim2(iMov_under100_sim2)), mean(movSimAll_sim2(iMov_over100_sim2))],60,'MarkerEdgeColor','none','MarkerFaceColor','k');
xlim([0 4]); ylim([0 1.05]*ylimPlot)
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([1,3]); xticklabels({'<100','>=100'})
xlabel('Intercellular distance (um)')
ylabel('Difference in adaptation ratio')
title('Simplex, movement')
a = ranksum(movSimAll_sim2(iMov_under100_sim2),movSimAll_sim2(iMov_over100_sim2));
text(2,.5,sprintf('p = %.3f',a));

subplot(2,3,6)
e1 = errorbar([1,3],[mean(rewSimAll_sim2(iRew_under100_sim2)), mean(rewSimAll_sim2(iRew_over100_sim2))], [std(rewSimAll_sim2(iRew_under100_sim2))/sqrt(numel(iRew_under100_sim2)), std(rewSimAll_sim2(iRew_over100_sim2))/sqrt(numel(iRew_over100_sim2))]);
e1.CapSize = 0; e1.LineWidth = 1.5; e1.Marker = 'none'; e1.Color = 'k';
hold on; axis square; box off
scatter([1,3],[mean(rewSimAll_sim2(iRew_under100_sim2)), mean(rewSimAll_sim2(iRew_over100_sim2))],60,'MarkerEdgeColor','none','MarkerFaceColor','k');
xlim([0 4]); ylim([0 1.05]*ylimPlot)
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([1,3]); xticklabels({'<100','>=100'})
xlabel('Intercellular distance (um)')
ylabel('Difference in adaptation ratio')
title('Simplex, reward')
a = ranksum(rewSimAll_sim2(iRew_under100_sim2),rewSimAll_sim2(iRew_over100_sim2));
text(2,.5,sprintf('p = %.3f',a));

suptitle('Distance-dependence of adaptation ratios');

%% Save maybe
if savebool
    figname = ['adaptRatiovsDist_all_',figtype];
    savefig(f1, fullfile(adaptFold,[figname,'.fig']));
    saveas(f1,fullfile(adaptFold,[figname,'.png']));
    print(f1,fullfile(adaptFold,[figname,'.eps']), '-depsc', '-vector');

    figname = ['adaptRatiovsDist_lobv_',figtype];
    savefig(f2, fullfile(adaptFold,[figname,'.fig']));
    saveas(f2,fullfile(adaptFold,[figname,'.png']));
    print(f2,fullfile(adaptFold,[figname,'.eps']), '-depsc', '-vector');

    figname = ['adaptRatiovsDist_sim2_',figtype];
    savefig(f3, fullfile(adaptFold,[figname,'.fig']));
    saveas(f3,fullfile(adaptFold,[figname,'.png']));
    print(f3,fullfile(adaptFold,[figname,'.eps']), '-depsc', '-vector');

    figname = ['adaptRatiovsDist_summary_',figtype];
    savefig(f4, fullfile(adaptFold,[figname,'.fig']));
    saveas(f4,fullfile(adaptFold,[figname,'.png']));
    print(f4,fullfile(adaptFold,[figname,'.eps']), '-depsc', '-vector');

    save(fullfile(adaptFold,sprintf('adaptDistSummary_%s.mat',figtype)),'adaptDistStruct');

end