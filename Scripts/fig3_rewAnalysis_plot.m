function figOut = fig3_rewAnalysis_plot(spkMat,xMat,vMat,imfs,vrfs,seshInfo, colMat)

nPCs = size(spkMat,2);
t2p = ((.5:1:4*imfs)/imfs-2)';
tVR = ((.5:1:4*vrfs)/vrfs-2)';

figOut = figure;
plotOffset = 0;
for iPlots = 1:size(spkMat,3)
    % Heatmap
    subplot(2,8,plotOffset+1)
    imagesc(spkMat(end/2-imfs+1:end/2+imfs,:,iPlots)',[0 10]);    
    colormap viridis; box off
    xline(imfs+.5,'--w','LineWidth',1);
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    xticks([0.5:imfs/2:2*imfs+.5]);
    yticks([0.5,19.5:20:nPCs+.5]);
    xticklabels({});
    
    if plotOffset == 0
        ylabel('PC dendrite number');
        yticklabels([1,20:20:nPCs]);
    else
        yticklabels({});
    end

    % Position
    subplot(4,8,17+plotOffset);
    h1 = boundedline(tVR,mean(xMat{iPlots},2),std(xMat{iPlots},[],2)/sqrt(size(xMat{iPlots},2)),'cmap',colMat(iPlots,:),'alpha');
    set(h1,'linewidth',1.5);
    xline(0,'--k','LineWidth',1);
    rectangle('Position',[-2 2/3 4 2/3],'FaceColor',[0 .5 .5 .2],'EdgeColor','none');
    xticks([-2:.5:2]); xticklabels({});
    yticks([0:2/3:2]); yticklabels({});
    axis([-1 1 0 2]);
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    if plotOffset == 0; ylabel('Object position'); end
    
    % Velocity
    subplot(4,8,25+plotOffset);
    h2 = boundedline(tVR,mean(vMat{iPlots},2),std(vMat{iPlots},[],2)/sqrt(size(vMat{iPlots},2)),'cmap',colMat(iPlots,:),'alpha');
    set(h2,'linewidth',1.5);
    xline(0,'--k','LineWidth',1);    
    xticks([-2:.5:2]); xticklabels({});
    yticks([-2.5:2.5:17.5]);
    axis([-1 1 -2.5 17.5]);
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    if plotOffset == 0
        ylabel({'Wheel velocity'; '(cm/s)'});
        yticklabels({'',0,'',5,'',10,'',15,''});
        xticklabels([-2:.5:2]);
        xtickangle(0);
        xlabel('Time (s)');
    else
        yticklabels({});
    end
    plotOffset = plotOffset+1;
    
    
    % Mean responses
    subplot(2,8,[6:8]);
    hold on;
    h3 = boundedline(t2p,mean(spkMat(:,:,iPlots),2),std(spkMat(:,:,iPlots),[],2)/sqrt(nPCs),'cmap',colMat(iPlots,:),'alpha');
    set(h3,'linewidth',1.5);
    subplot(4,8,[6:8]+16);
    hold on;
    h4 = boundedline(tVR,mean(xMat{iPlots},2),std(xMat{iPlots},[],2)/sqrt(size(xMat{iPlots},2)),'cmap',colMat(iPlots,:),'alpha');
    set(h4,'linewidth',1.5);
    subplot(4,8,[6:8]+24);
    hold on;
    h5 = boundedline(tVR,mean(vMat{iPlots},2),std(vMat{iPlots},[],2)/sqrt(size(vMat{iPlots},2)),'cmap',colMat(iPlots,:),'alpha');
    set(h5,'linewidth',1.5);
    
    if iPlots == size(spkMat,3)
        subplot(2,8,[6:8]);
        xline(0,'--k','LineWidth',1);
        axis([-1 .5 0 15]);
        xticks([-2:.5:2]); xticklabels({});
        yticks([0:3:15]);
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
%         ylabel('Mean field response (spk/s)');
        text(.75,8,seshInfo.range,'FontSize',12,'HorizontalAlignment','center');
        title(sprintf('%s - %s, sorted by %s',seshInfo.name,seshInfo.fov,seshInfo.splitType));
        
        subplot(4,8,[6:8]+16);
        xline(0,'--k','LineWidth',1);
        rectangle('Position',[-2 2/3 4 2/3],'FaceColor',[0 .5 .5 .2],'EdgeColor','none');
        xticks([-2:.5:2]); xticklabels({});
        yticks([0:2/3:2]); yticklabels({});
        axis([-1 .5 0 2]);
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
%         ylabel('Object position');
        
        subplot(4,8,[6:8]+24);
        xline(0,'--k','LineWidth',1);    
        xticks([-2:.5:2]);
        yticks([-2.5:2.5:17.5]); yticklabels({'',0,'',5,'',10,'',15,''});
        axis([-1 .5 -2.5 17.5]);
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
%         ylabel({'Wheel velocity (cm/s)'});
    end
    
end



end