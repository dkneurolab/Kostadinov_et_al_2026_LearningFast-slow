function matchSpont = fig2_roiMatchSpont(matchData,dataFold,paperFold,matchColors, savebool)
%%
alignFold = fullfile(dataFold,'01_Behav+imaging');
figFold = fullfile(paperFold,'Figures','Fig1','Sub3b_spontFigs');
if ~exist(figFold,'dir'); mkdir(figFold); end
figFold2 = fullfile(dataFold,'paperFigs','Fig1','Sub3b_spontFigs');
if ~exist(figFold2,'dir'); mkdir(figFold2); end

%% Load in data (separate from figure-making so that we don't have to go into big files all the time)
if ~exist(fullfile(figFold2,'matchDataSpont.mat'),'file')
    % If this has not already been done, iterate through session pairs and
    % extract spont zF and spks
    
    matchSpont = struct;
    for i = 1:numel(matchData)
        
        if isempty(matchData(i).name)
            continue
        else
            matchSpont(i).name = matchData(i).name;
            matchSpont(i).fov = matchData(i).fov;
            matchSpont(i).matchColor = matchData(i).matchColor;
            matchSpont(i).date1 = matchData(i).date1;
            matchSpont(i).date2 = matchData(i).date2;
            matchSpont(i).v1 = matchData(i).v1;
            matchSpont(i).v2 = matchData(i).v2;
            matchSpont(i).iROIs1 = matchData(i).iROIs1;
            matchSpont(i).iROIs2 = matchData(i).iROIs2;
            
            % Spont data for fov1
            sets1 = load(fullfile(alignFold,sprintf('Version_%i',matchData(i).v1),matchSpont(i).name,...
                sprintf('%s_%s',matchData(i).date1,matchData(i).fov),'setsonly.mat')); 
            matchSpont(i).sets1 = sets1.sets;
            matchSpont(i).spk1All = reshape(permute(matchData(i).on1.spk(1:end/2,:,:),[1,3,2]),[],numel(sets1.sets.pcaindclustall));
            matchSpont(i).zF1All = sets1.sets.Fspont;
%             % Load in spont data for fov1
%             eventsdata1 = load(fullfile(alignFold,sprintf('Version_%i',matchData(i).v1),matchSpont(i).name,...
%                 sprintf('%s_%s',matchData(i).date1,matchData(i).fov),'eventsdata.mat'));
%             zF1All = []; spk1All = [];
%             if isfield(eventsdata1.savedata.data_trials.F,'zFONL')
%                 zF1All = cat(3,zF1All,eventsdata1.savedata.data_trials.F.zFONL);
%                 spk1All = cat(3,spk1All,eventsdata1.savedata.data_trials.spk.spkONL);
%             end
%             if isfield(eventsdata1.savedata.data_trials.F,'zFONR')
%                 zF1All = cat(3,zF1All,eventsdata1.savedata.data_trials.F.zFONR);
%                 spk1All = cat(3,spk1All,eventsdata1.savedata.data_trials.spk.spkONR);
%             end
%             sets1 = eventsdata1.sets;
%             zF1All = reshape(permute(zF1All(end/2-2*sets1.imfs+1:end/2,:,:),[1,3,2]),[],size(zF1All,2));
%             spk1All = reshape(permute(spk1All(end/2-2*sets1.imfs+1:end/2,:,:),[1,3,2]),[],size(spk1All,2));
            
%             matchSpont(i).sets1 = sets1;
%             matchSpont(i).zF1All = zF1All;
%             matchSpont(i).spk1All = spk1All;
% %             matchSpont(i).zF1 = zF1All(:,matchSpont(i).iROIs1);
% %             matchSpont(i).spk1 = spk1All(:,matchSpont(i).iROIs1);
            clear eventsdata1 sets1 zF1All spk1All
            
            % Spont data for fov2
            sets2 = load(fullfile(alignFold,sprintf('Version_%i',matchData(i).v2),matchSpont(i).name,...
                sprintf('%s_%s',matchData(i).date2,matchData(i).fov),'setsonly.mat')); 
            matchSpont(i).sets2 = sets2.sets;
            matchSpont(i).spk2All = reshape(permute(matchData(i).on2.spk(1:end/2,:,:),[1,3,2]),[],numel(sets2.sets.pcaindclustall));
            matchSpont(i).zF2All = sets2.sets.Fspont;
%             % Load in spont data for fov2
%             eventsdata2 = load(fullfile(alignFold,sprintf('Version_%i',matchData(i).v2),matchSpont(i).name,...
%                 sprintf('%s_%s',matchData(i).date2,matchData(i).fov),'eventsdata.mat'));
%             zF2All = []; spk2All = [];
%             if isfield(eventsdata2.savedata.data_trials.F,'zFONL')
%                 zF2All = cat(3,zF2All,eventsdata2.savedata.data_trials.F.zFONL);
%                 spk2All = cat(3,spk2All,eventsdata2.savedata.data_trials.spk.spkONL);
%             end
%             if isfield(eventsdata2.savedata.data_trials.F,'zFONR')
%                 zF2All = cat(3,zF2All,eventsdata2.savedata.data_trials.F.zFONR);
%                 spk2All = cat(3,spk2All,eventsdata2.savedata.data_trials.spk.spkONR);
%             end
%             sets2 = eventsdata2.sets;
%             zF2All = reshape(permute(zF2All(end/2-2*sets2.imfs+1:end/2,:,:),[1,3,2]),[],size(zF2All,2));
%             spk2All = reshape(permute(spk2All(end/2-2*sets2.imfs+1:end/2,:,:),[1,3,2]),[],size(spk2All,2));
%             
%             matchSpont(i).sets2 = sets2;
%             matchSpont(i).zF2All = zF2All;
%             matchSpont(i).spk2All = spk2All;
%             %     matchSpont(i).zF2 = zF2All(:,matchSpont(i).iROIs2);
%             %     matchSpont(i).spk2 = spk2All(:,matchSpont(i).iROIs2);
            
            clear eventsdata2 sets2 zF2All spk2All
            
            %
            
            matchSpont(i).spkRate1 = sum(logical(matchSpont(i).spk1All))'/(size(matchSpont(i).spk1All,1)/matchSpont(i).sets1.imfs);
            matchSpont(i).spkRate2 = sum(logical(matchSpont(i).spk2All))'/(size(matchSpont(i).spk2All,1)/matchSpont(i).sets2.imfs);
            
            if (matchSpont(i).v1 == 1 && matchSpont(i).v2 == 4) || (matchSpont(i).v1 == 4 && matchSpont(i).v2 == 4)
                [~,iSort] = sort(matchSpont(i).iROIs2,'ascend');
            elseif matchSpont(i).v1 == 4 && matchSpont(i).v2 == 5
                [~,iSort] = sort(matchSpont(i).iROIs1,'ascend');
            end
            matchSpont(i).R1 = corrcoef(matchSpont(i).zF1All(:,matchSpont(i).iROIs1(iSort))); matchSpont(i).R1(1:size(matchSpont(i).R1,1)+1:end) = nan;
            matchSpont(i).R2 = corrcoef(matchSpont(i).zF2All(:,matchSpont(i).iROIs2(iSort))); matchSpont(i).R2(1:size(matchSpont(i).R2,1)+1:end) = nan;
            Rbool = logical(triu(ones(size(matchSpont(i).R1)),1));
            
            R1lin = matchSpont(i).R1(Rbool);
            R2lin = matchSpont(i).R2(Rbool);
            R2lin_1off = circshift(R2lin,1);
            rng(69);
            R2lin_rand = R2lin(randperm(numel(R2lin)));
            
            corr12 = corrcoef(R1lin,R2lin);
            corr12off = corrcoef(R1lin,R2lin_1off);
            corr12rand = corrcoef(R1lin,R2lin_rand);
            matchSpont(i).Rmat(i,1:3) = [corr12(1,2) corr12off(1,2) corr12rand(1,2)];            
            clear corr12 corr12off corr12rand
            
        end
        
    end
    
    if savebool
        save(fullfile(figFold2,'matchDataSpont.mat'),'matchSpont','-v7.3');
    end
    
else
    matchSpont = load(fullfile(figFold2,'matchDataSpont.mat')); matchSpont = matchSpont.matchSpont;
end

%% Plotting figures

corrSummary = struct;
for i = 1:numel(matchSpont)
    
    if isempty(matchSpont(i).name)
        continue
    else
        nRois = size(matchSpont(i).R1,1);
        Rbool = logical(triu(ones(size(matchSpont(i).R1)),1));
        R1lin = matchSpont(i).R1(Rbool);
        R2lin = matchSpont(i).R2(Rbool);
        
        corrSummary(i).delDays = abs(days(datetime(matchSpont(i).date2)-datetime(matchSpont(i).date1)));
        if strcmpi(matchSpont(i).fov,'lobv'); corrSummary(i).ilobv = true; else corrSummary(i).ilobv = false; end
        if strcmpi(matchSpont(i).fov,'sim2'); corrSummary(i).isim2 = true; else corrSummary(i).isim2 = false; end
        corrSummary(i).Rmat = matchSpont(i).Rmat(i,:);
        corrSummary(i).spkDel = mean(matchSpont(i).spkRate2(matchSpont(i).iROIs2)) - mean(matchSpont(i).spkRate1(matchSpont(i).iROIs1));
        
        
        fig1 = figure;
        subplot(2,2,1);
        imagesc(matchSpont(i).R1); axis square; colormap viridis; box off;
        caxis([prctile(R1lin,1) prctile(R1lin,99)]);
        title(sprintf('%s, %s, version %i on %s',upper(matchSpont(i).name),matchSpont(i).fov,matchSpont(i).v1,matchSpont(i).date1));
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
        xticks([1,20:20:nRois]); yticks([1,20:20:nRois]);
        
        subplot(2,2,3);
        imagesc(matchSpont(i).R2); axis square; colormap viridis; box off;
        caxis([prctile(R2lin,1) prctile(R2lin,99)]);
        title(sprintf('%s, %s, version %i on %s',upper(matchSpont(i).name),matchSpont(i).fov,matchSpont(i).v2,matchSpont(i).date2));
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
        xticks([1,20:20:nRois]); yticks([1,20:20:nRois]);
        
        subplot(2,2,2);
        hold on;
        plot([-.2 1],[-.2 1],'--b','LineWidth',1)
        scatter(R1lin,R2lin,30,'MarkerFaceColor','k','MarkerFaceAlpha',.1,'MarkerEdgeColor','none');
        xlabel(sprintf('Corr. coeffs: Version %i',matchSpont(i).v1));
        ylabel(sprintf('Corr. coeffs: Version %i',matchSpont(i).v2));
        axis([-.2 1 -.2 1]); axis square;
        box off;
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
        xticks([-.2:.2:1]); yticks([-.2:.2:1]);
        
        subplot(2,2,4);
        hold on;
        bar(1,mean(matchSpont(i).spkRate1(matchSpont(i).iROIs1)),0.6,'EdgeColor','none','FaceColor',[.6 .6 .6],'FaceAlpha',0.5);
        errorbar(1, mean(matchSpont(i).spkRate1(matchSpont(i).iROIs1)), std(matchSpont(i).spkRate1(matchSpont(i).iROIs1)),'k','CapSize',0,'LineWidth',2);
        bar(4,mean(matchSpont(i).spkRate2(matchSpont(i).iROIs2)),0.6,'EdgeColor','none','FaceColor',[.6 .6 .6],'FaceAlpha',0.5);
        errorbar(4, mean(matchSpont(i).spkRate2(matchSpont(i).iROIs2)), std(matchSpont(i).spkRate2(matchSpont(i).iROIs2)),'k','CapSize',0,'LineWidth',2);
        plot([2*ones(1,nRois); 3*ones(1,nRois)], [matchSpont(i).spkRate1(matchSpont(i).iROIs1) matchSpont(i).spkRate2(matchSpont(i).iROIs2)]','LineWidth',1,'Color',[0 0 0 0.5]);
        ylim([0 2.5]); yticks([0:0.5:2.5]);
        xlim([0 5]); xticks([1,4]); xticklabels({sprintf('Version %i',matchSpont(i).v1),sprintf('Version %i',matchSpont(i).v2)})
        ylabel('Firing rate (spks/s)');
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
        axis square;
        
        
        
        
        if savebool
            figname = sprintf('%s_%s_v%i_v%i_corrs',matchData(i).name,matchData(i).fov,matchData(i).v1,matchData(i).v2);
            figFoldfov1 = fullfile(figFold,sprintf('%s_%s_v%i_v%i',matchData(i).name,matchData(i).fov,matchData(i).v1,matchData(i).v2));
            if ~exist(figFoldfov1,'dir'); mkdir(figFoldfov1); end
            figFoldfov2 = fullfile(figFold2,sprintf('%s_%s_v%i_v%i',matchData(i).name,matchData(i).fov,matchData(i).v1,matchData(i).v2));
            if ~exist(figFoldfov2,'dir'); mkdir(figFoldfov2); end
            
            savefig(fig1, fullfile(figFoldfov1,[figname,'.fig']));
            saveas(fig1,fullfile(figFoldfov1,[figname,'.png']));
            print(fig1,fullfile(figFoldfov1,[figname,'.eps']), '-depsc', '-painters');
            savefig(fig1, fullfile(figFoldfov2,[figname,'.fig']));
            saveas(fig1,fullfile(figFoldfov2,[figname,'.png']));
            print(fig1,fullfile(figFoldfov2,[figname,'.eps']), '-depsc', '-painters');
            
        end
    end
    
end

%%
delDays = cat(1,corrSummary.delDays);
spkDel = cat(1,corrSummary.spkDel);
ilobv = logical(cat(1,corrSummary.ilobv));
isim2 = logical(cat(1,corrSummary.isim2));
Rmat = cat(1,corrSummary.Rmat);

% Linear fits - spike change:
oSpkLV = regress_hd(delDays(ilobv),spkDel(ilobv),1);
oSpkLS = regress_hd(delDays(isim2),spkDel(isim2),1);
oSpk = regress_hd(delDays,spkDel,1);

% Linear fits - corr change:
oCorrLV = regress_hd(delDays(ilobv),Rmat(ilobv,1),1);
oCorrLS = regress_hd(delDays(isim2),Rmat(isim2,1),1);
oCorr = regress_hd(delDays,Rmat(:,1),1);

fig2 = figure;
subplot(2,3,1); hold on;
h1 = boundedline(oSpkLV.x_fit,oSpkLV.y_fit,oSpkLV.delta,'cmap',[.5 0 .5],'alpha');
h2 = boundedline(oSpkLS.x_fit,oSpkLS.y_fit,oSpkLS.delta,'cmap',[0 .5 .5],'alpha');
set(h1,'linewidth',2); set(h2,'linewidth',1.5);
s1 = scatter(delDays(ilobv),spkDel(ilobv),108,[.5 0 .5],'filled');
s1.MarkerEdgeColor = 'none'; s1.MarkerFaceAlpha = 0.5;
s2 = scatter(delDays(isim2),spkDel(isim2),108,[0 .5 .5],'filled');
s2.MarkerEdgeColor = 'none'; s2.MarkerFaceAlpha = 0.5;
axis square;
xlim([-1.5 75]); ylim([-1 1]);
xticks([0:15:75]); yticks([-1:0.25:1]);
ylabel({'Mean firing rate change','Early higher   Later higher'});
text(60,.75,{sprintf('Rsq LV = %0.3f, p = %0.3f',oSpkLV.rsq,oSpkLV.corr_p),sprintf('Rsq LS = %0.3f, p = %0.3f',oSpkLS.rsq,oSpkLS.corr_p)})
xlabel('Days between sessions');
title(['\color[rgb]{.5 0 .5} LobV','\color{black} and','\color[rgb]{0 .5 .5} Simplex','\color{black} fields of view'])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');

subplot(2,3,4); hold on;
h1 = boundedline(oSpk.x_fit,oSpk.y_fit,oSpk.delta,'cmap',[0 0 0],'alpha');
set(h1,'linewidth',2);
s1 = scatter(delDays,spkDel,108,[0 0 0],'filled');
s1.MarkerEdgeColor = 'none'; s1.MarkerFaceAlpha = 0.5;
axis square;
xlim([-1.5 75]); ylim([-1 1]);
xticks([0:15:75]); yticks([-1:0.25:1]);
ylabel({'Mean firing rate change','Early higher   Later higher'});
text(60,.75,sprintf('Rsq = %0.3f, p = %0.3f',oSpk.rsq,oSpk.corr_p));
xlabel('Days between sessions');
title(['Combined fields of view'])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');

subplot(2,3,2); hold on;
h1 = boundedline(oCorrLV.x_fit,oCorrLV.y_fit,oCorrLV.delta,'cmap',[.5 0 .5],'alpha');
h2 = boundedline(oCorrLS.x_fit,oCorrLS.y_fit,oCorrLS.delta,'cmap',[0 .5 .5],'alpha');
set(h1,'linewidth',2); set(h2,'linewidth',1.5);
s1 = scatter(delDays(ilobv),Rmat(ilobv,1),108,[.5 0 .5],'filled');
s1.MarkerEdgeColor = 'none'; s1.MarkerFaceAlpha = 0.5;
s2 = scatter(delDays(isim2),Rmat(isim2,1),108,[0 .5 .5],'filled');
s2.MarkerEdgeColor = 'none'; s2.MarkerFaceAlpha = 0.5;
axis square;
xlim([-1.5 75]); ylim([-.05 1]);
xticks([0:15:75]); yticks([0:0.2:1]);
ylabel('Correlation matrix similarity');
text(60,.75,{sprintf('Rsq LV = %0.3f, p = %0.3f',oCorrLV.rsq,oCorrLV.corr_p),sprintf('Rsq LS = %0.3f, p = %0.3f',oCorrLS.rsq,oCorrLS.corr_p)})
xlabel('Days between sessions');
title(['\color[rgb]{.5 0 .5} LobV','\color{black} and','\color[rgb]{0 .5 .5} Simplex','\color{black} fields of view'])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');

subplot(2,3,5); hold on;
h1 = boundedline(oCorr.x_fit,oCorr.y_fit,oCorr.delta,'cmap',[0 0 0],'alpha');
set(h1,'linewidth',2);
s1 = scatter(delDays,Rmat(:,1),108,[0 0 0],'filled');
s1.MarkerEdgeColor = 'none'; s1.MarkerFaceAlpha = 0.5;
axis square;
xlim([-1.5 75]); ylim([-.05 1]);
xticks([0:15:75]); yticks([0:0.2:1]);
ylabel('Correlation matrix similarity');
text(60,.75,sprintf('Rsq = %0.3f, p = %0.3f',oCorr.rsq,oCorr.corr_p));
xlabel('Days between sessions');
title(['Combined fields of view'])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');

subplot(2,3,3); hold on;
s1 = scatter(delDays(ilobv),Rmat(ilobv,1),108,[.5 0 .5],'filled');
s1.MarkerEdgeColor = 'none'; s1.MarkerFaceAlpha = 0.5;
s2 = scatter(delDays(ilobv),Rmat(ilobv,2),108,[.5 .5 .5],'filled');
s2.MarkerEdgeColor = 'none'; s2.MarkerFaceAlpha = 0.5;
s3 = scatter(delDays(ilobv),Rmat(ilobv,3),108,[0 0 0],'filled');
s3.MarkerEdgeColor = 'none'; s3.MarkerFaceAlpha = 0.5;
axis square;
xlim([-1.5 75]); ylim([-.05 1]);
xticks([0:15:75]); yticks([0:0.2:1]);
ylabel('Correlation matrix similarity');
xlabel('Days between sessions');
title('Lobule V fields of view');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');

subplot(2,3,6); hold on;
s1 = scatter(delDays(isim2),Rmat(isim2,1),108,[0 .5 .5],'filled');
s1.MarkerEdgeColor = 'none'; s1.MarkerFaceAlpha = 0.5;
s2 = scatter(delDays(isim2),Rmat(isim2,2),108,[.5 .5 .5],'filled');
s2.MarkerEdgeColor = 'none'; s2.MarkerFaceAlpha = 0.5;
s3 = scatter(delDays(isim2),Rmat(isim2,3),108,[0 0 0],'filled');
s3.MarkerEdgeColor = 'none'; s3.MarkerFaceAlpha = 0.5;
axis square;
xlim([-1.5 75]); ylim([-.05 1]);
xticks([0:15:75]); yticks([0:0.2:1]);
ylabel('Correlation matrix similarity');
xlabel('Days between sessions');
title('Simplex fields of view');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');

if savebool
    figname = ['corrSummary'];
    savefig(fig2, fullfile(figFold,[figname,'.fig']));
    saveas(fig2,fullfile(figFold,[figname,'.png']));
    print(fig2,fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    savefig(fig2, fullfile(figFold2,[figname,'.fig']));
    saveas(fig2,fullfile(figFold2,[figname,'.png']));
    print(fig2,fullfile(figFold2,[figname,'.eps']), '-depsc', '-painters');
    save(fullfile(figFold2,'corrSummary.mat'),'corrSummary');
    
end