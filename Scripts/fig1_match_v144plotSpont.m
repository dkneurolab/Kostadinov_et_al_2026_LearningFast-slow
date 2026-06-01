function lobData = fig1_match_v144plotSpont(lobData, dataFold, paperFold, savebool)

for i = 1:numel(lobData)
    nMatched = size(lobData(i).v144rois,1);
    for j = 1:3
        lobData(i).v144spkRate{j} = lobData(i).v144spkRate{j}(lobData(i).v144rois(:,j));
        lobData(i).v144corrMat{j} = lobData(i).v144corrMat{j}(lobData(i).v144rois(:,j),lobData(i).v144rois(:,j));
        lobData(i).v144corrMat{j}(1:nMatched+1:end) = nan;
    end

    lobData(i).v144spkRate = cell2mat(lobData(i).v144spkRate');
    lobData(i).v144corrMat = reshape(cell2mat(lobData(i).v144corrMat)',[nMatched, nMatched, 3]);
        
    Rbool = logical(triu(ones(nMatched,nMatched),1));
    R1lin = lobData(i).v144corrMat(:,:,1); R1lin = R1lin(Rbool);
    R4lin = lobData(i).v144corrMat(:,:,2); R4lin = R4lin(Rbool);
    R5lin = lobData(i).v144corrMat(:,:,3); R5lin = R5lin(Rbool);
    
    R1lin_1off = circshift(R1lin,1);
    R5lin_1off = circshift(R5lin,1);
    rng(69);
    R1lin_rand = R1lin(randperm(numel(R1lin)));
    rng(69);
    R5lin_rand = R5lin(randperm(numel(R5lin)));
    
    corr14 = corrcoef(R4lin,R1lin);
    corr45 = corrcoef(R4lin,R5lin);
    corr14off = corrcoef(R4lin,R1lin_1off);
    corr45off = corrcoef(R4lin,R5lin_1off);
    corr14rand = corrcoef(R4lin,R1lin_rand);
    corr45rand = corrcoef(R4lin,R5lin_rand);
    
    lobData(i).corrPlot = [corr14(1,2) corr14off(1,2) corr14rand(1,2); corr45(1,2) corr45off(1,2) corr45rand(1,2)]'; 
    
    clear R1lin_1off R5lin_1off R1lin_rand R5lin_rand corr14 corr45 corr14off corr45off corr14rand corr45rand

    figout(i) = figure;
    subplot(2,3,1);
    imagesc(lobData(i).v144corrMat(:,:,1)); axis square; colormap viridis; box off;
    caxis([prctile(R1lin,1) prctile(R1lin,99)]);
    title(sprintf('%s, %s, Naive',upper(lobData(i).name),lobData(i).fov));
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
    xticks([1,20:20:nMatched]); yticks([1,20:20:nMatched]);
    ylabel('PC dendrite number');
    
    subplot(2,3,2);
    imagesc(lobData(i).v144corrMat(:,:,2)); axis square; colormap viridis; box off;
    caxis([prctile(R4lin,1) prctile(R4lin,99)]);
    title(sprintf('%s, %s, Learning',upper(lobData(i).name),lobData(i).fov));
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
    xticks([1,20:20:nMatched]); yticks([1,20:20:nMatched]);
    ylabel('PC dendrite number');
    
    subplot(2,3,3);
    imagesc(lobData(i).v144corrMat(:,:,3)); axis square; colormap viridis; box off;
    caxis([prctile(R5lin,1) prctile(R5lin,99)]);
    title(sprintf('%s, %s, Trained',upper(lobData(i).name),lobData(i).fov));
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
    xticks([1,20:20:nMatched]); yticks([1,20:20:nMatched]);
    ylabel('PC dendrite number');

    subplot(2,2,3);
    hold on;
    bar(1,mean(lobData(i).v144spkRate(:,1)),0.6,'EdgeColor','none','FaceColor',[.6 .6 .6],'FaceAlpha',0.5);
    errorbar(1, mean(lobData(i).v144spkRate(:,1)), std(lobData(i).v144spkRate(:,1)),'k','CapSize',0,'LineWidth',2);
    bar(4,mean(lobData(i).v144spkRate(:,2)),0.6,'EdgeColor','none','FaceColor',[.6 .6 .6],'FaceAlpha',0.5);
    errorbar(4, mean(lobData(i).v144spkRate(:,2)), std(lobData(i).v144spkRate(:,2)),'k','CapSize',0,'LineWidth',2);
    bar(7,mean(lobData(i).v144spkRate(:,3)),0.6,'EdgeColor','none','FaceColor',[.6 .6 .6],'FaceAlpha',0.5);
    errorbar(7, mean(lobData(i).v144spkRate(:,3)), std(lobData(i).v144spkRate(:,3)),'k','CapSize',0,'LineWidth',2);
    plot([2*ones(1,nMatched); 3*ones(1,nMatched)], [lobData(i).v144spkRate(:,1) lobData(i).v144spkRate(:,2)]','LineWidth',1,'Color',[0 0 0 0.5]);
    plot([5*ones(1,nMatched); 6*ones(1,nMatched)], [lobData(i).v144spkRate(:,2) lobData(i).v144spkRate(:,3)]','LineWidth',1,'Color',[0 0 0 0.5]);
    
    ylim([0 2.5]); yticks([0:0.5:2.5]);
    xlim([0 8]); xticks([1,4,7]); xticklabels({'Naive','Learning','Trained'})
    ylabel('Firing rate (spks/s)');
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
    axis square;
    
    subplot(2,2,4);
    hold on;
    plot([1 4 7],lobData(i).corrPlot(:,1),'b','LineWidth',2)
    plot([1 4 7],lobData(i).corrPlot(:,2),'r','LineWidth',2)
    ylim([-.05 1]); yticks([0:0.2:1]);
    xlim([0 8]); xticks([1,4,7]); xticklabels({'Matched','Neighbor','Random'}); xtickangle(0);
    ylabel('Correlation matrix similarity');
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
    axis square;
    legend({'Naive-Trained','Learning-Trained'})
    
    if savebool
        figFold = fullfile(paperFold,'Figures','Fig4','Sub4b_summarySpont');
        if ~exist(figFold,'dir'); mkdir(figFold); end
        figFold2 = fullfile(dataFold,'paperFigs','Fig4','Sub4b_summarySpont');
        if ~exist(figFold2,'dir'); mkdir(figFold2); end
        figname = sprintf('%s_%s_v144corrs',lobData(i).name,lobData(i).fov);
        savefig(figout(i), fullfile(figFold,[figname,'.fig']));
        saveas(figout(i),fullfile(figFold,[figname,'.png']));
        print(figout(i),fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
        savefig(figout(i), fullfile(figFold2,[figname,'.fig']));
        saveas(figout(i),fullfile(figFold2,[figname,'.png']));
        print(figout(i),fullfile(figFold2,[figname,'.eps']), '-depsc', '-painters');

    end
    
end


end