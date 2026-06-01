function fig3_v4UCO(seshv4,dataFold,paperFold,savebool)
% Plot UCO trials for each v4n FOV with error bars and then plot mean response across mice

%% Settings
figFold = fullfile(paperFold,'Figures','Fig3','Sub1_behavUCO');
if ~exist(figFold,'dir'); mkdir(figFold); end
seshv4 = rmfield(seshv4,{'vrnum','paqnum','data2p','dataVid','fpblock','numplanes','version'});



%% Load in each UCO file
for iSesh = 1:numel(seshv4)
    if seshv4(iSesh).good
        % Load in dataset
        dataPath = fullfile(dataFold,'Version_4',seshv4(iSesh).name,sprintf('%s_%s',seshv4(iSesh).date,seshv4(iSesh).fov));
        trialsOnly = load(fullfile(dataPath,'trialsonly.mat'));
        trialData = trialsOnly.trialstructs;
        trialSets = trialsOnly.sets; clear trialsOnly
        imfs = trialSets.imfs;
        
        % Extract relevant data
        if strcmpi(seshv4(iSesh).fov,'lobv')
            seshv4(iSesh).iLobv = true;
        else
            seshv4(iSesh).iLobv = false;
        end
        
        seshv4(iSesh).nRoi = size(trialSets.Fspont,2);
        seshv4(iSesh).mov.zfUnder = cat(3, trialData.LUstruct.zFLmov);
        seshv4(iSesh).mov.zfCorr = cat(3, trialData.LCstruct.zFLmov);
        seshv4(iSesh).mov.zfOver = cat(3, trialData.LOstruct.zFLmov);
        seshv4(iSesh).mov.spkUnder = cat(3, trialData.LUstruct.spkLmov);
        seshv4(iSesh).mov.spkCorr = cat(3, trialData.LCstruct.spkLmov);
        seshv4(iSesh).mov.spkOver = cat(3, trialData.LOstruct.spkLmov);
        seshv4(iSesh).off.zfUnder = cat(3, trialData.LUstruct.zFLoff);
        seshv4(iSesh).off.zfCorr = cat(3, trialData.LCstruct.zFLoff);
        seshv4(iSesh).off.zfOver = cat(3, trialData.LOstruct.zFLoff);
        seshv4(iSesh).off.spkUnder = cat(3, trialData.LUstruct.spkLoff);
        seshv4(iSesh).off.spkCorr = cat(3, trialData.LCstruct.spkLoff);
        seshv4(iSesh).off.spkOver = cat(3, trialData.LOstruct.spkLoff);
        
        % Average across trials and cells and clip to +/- 1s
        seshv4(iSesh).zfMovUCOmu = [mean(mean(seshv4(iSesh).mov.zfUnder(end/2-imfs+1:end/2+imfs,:,:),2),3) ...
            mean(mean(seshv4(iSesh).mov.zfCorr(end/2-imfs+1:end/2+imfs,:,:),2),3) ...
            mean(mean(seshv4(iSesh).mov.zfOver(end/2-imfs+1:end/2+imfs,:,:),2),3)];
        seshv4(iSesh).zfMovUCOsd = [std(mean(seshv4(iSesh).mov.zfUnder(end/2-imfs+1:end/2+imfs,:,:),2),[],3) ...
            std(mean(seshv4(iSesh).mov.zfCorr(end/2-imfs+1:end/2+imfs,:,:),2),[],3) ...
            std(mean(seshv4(iSesh).mov.zfOver(end/2-imfs+1:end/2+imfs,:,:),2),[],3)];
        seshv4(iSesh).spkMovUCOmu = imfs*[mean(mean(seshv4(iSesh).mov.spkUnder(end/2-imfs+1:end/2+imfs,:,:),2),3) ...
            mean(mean(seshv4(iSesh).mov.spkCorr(end/2-imfs+1:end/2+imfs,:,:),2),3) ...
            mean(mean(seshv4(iSesh).mov.spkOver(end/2-imfs+1:end/2+imfs,:,:),2),3)];
        seshv4(iSesh).spkMovUCOsd = imfs*[std(mean(seshv4(iSesh).mov.spkUnder(end/2-imfs+1:end/2+imfs,:,:),2),[],3) ...
            std(mean(seshv4(iSesh).mov.spkCorr(end/2-imfs+1:end/2+imfs,:,:),2),[],3) ...
            std(mean(seshv4(iSesh).mov.spkOver(end/2-imfs+1:end/2+imfs,:,:),2),[],3)];
        
        seshv4(iSesh).zfOffUCOmu = [mean(mean(seshv4(iSesh).off.zfUnder(end/2-imfs+1:end/2+imfs,:,:),2),3) ...
            mean(mean(seshv4(iSesh).off.zfCorr(end/2-imfs+1:end/2+imfs,:,:),2),3) ...
            mean(mean(seshv4(iSesh).off.zfOver(end/2-imfs+1:end/2+imfs,:,:),2),3)];
        seshv4(iSesh).zfOffUCOsd = [std(mean(seshv4(iSesh).off.zfUnder(end/2-imfs+1:end/2+imfs,:,:),2),[],3) ...
            std(mean(seshv4(iSesh).off.zfCorr(end/2-imfs+1:end/2+imfs,:,:),2),[],3) ...
            std(mean(seshv4(iSesh).off.zfOver(end/2-imfs+1:end/2+imfs,:,:),2),[],3)];
        seshv4(iSesh).spkOffUCOmu = imfs*[mean(mean(seshv4(iSesh).off.spkUnder(end/2-imfs+1:end/2+imfs,:,:),2),3) ...
            mean(mean(seshv4(iSesh).off.spkCorr(end/2-imfs+1:end/2+imfs,:,:),2),3) ...
            mean(mean(seshv4(iSesh).off.spkOver(end/2-imfs+1:end/2+imfs,:,:),2),3)];
        seshv4(iSesh).spkOffUCOsd = imfs*[std(mean(seshv4(iSesh).off.spkUnder(end/2-imfs+1:end/2+imfs,:,:),2),[],3) ...
            std(mean(seshv4(iSesh).off.spkCorr(end/2-imfs+1:end/2+imfs,:,:),2),[],3) ...
            std(mean(seshv4(iSesh).off.spkOver(end/2-imfs+1:end/2+imfs,:,:),2),[],3)];
        
        tPlot = [((1:2*imfs)/imfs-1/imfs/2)-1]';
        % Plot!
        f1 = figure;
        subplot(2,2,1); % zF, mov
        h1 = boundedline(tPlot,seshv4(iSesh).zfMovUCOmu(:,1),seshv4(iSesh).zfMovUCOsd(:,1)/sqrt(seshv4(iSesh).nRoi),'cmap',[0 .5 0],'alpha');
        set(h1,'linewidth',1.5);
        h2 = boundedline(tPlot,seshv4(iSesh).zfMovUCOmu(:,3),seshv4(iSesh).zfMovUCOsd(:,3)/sqrt(seshv4(iSesh).nRoi),'cmap',[.8 .5 0],'alpha');
        set(h2,'linewidth',1.5);
        h3 = boundedline(tPlot,seshv4(iSesh).zfMovUCOmu(:,2),seshv4(iSesh).zfMovUCOsd(:,2)/sqrt(seshv4(iSesh).nRoi),'cmap',[0 0 0],'alpha');
        set(h3,'linewidth',1.5);
        axis([-1 1 -.5 1.5]); axis square;
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
        title('Fluo, movement')
        
        subplot(2,2,2); % zF, off
        h1 = boundedline(tPlot,seshv4(iSesh).zfOffUCOmu(:,1),seshv4(iSesh).zfOffUCOsd(:,1)/sqrt(seshv4(iSesh).nRoi),'cmap',[0 .5 0],'alpha');
        set(h1,'linewidth',1.5);
        h2 = boundedline(tPlot,seshv4(iSesh).zfOffUCOmu(:,3),seshv4(iSesh).zfOffUCOsd(:,3)/sqrt(seshv4(iSesh).nRoi),'cmap',[.8 .5 0],'alpha');
        set(h2,'linewidth',1.5);
        h3 = boundedline(tPlot,seshv4(iSesh).zfOffUCOmu(:,2),seshv4(iSesh).zfOffUCOsd(:,2)/sqrt(seshv4(iSesh).nRoi),'cmap',[0 0 0],'alpha');
        set(h3,'linewidth',1.5);
        axis([-1 1 -.5 1.5]); axis square;
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
        title('Fluo, offset')
        
        subplot(2,2,3); % spk, mov
        h1 = boundedline(tPlot,seshv4(iSesh).spkMovUCOmu(:,1),seshv4(iSesh).spkMovUCOsd(:,1)/sqrt(seshv4(iSesh).nRoi),'cmap',[0 .5 0],'alpha');
        set(h1,'linewidth',1.5);
        h2 = boundedline(tPlot,seshv4(iSesh).spkMovUCOmu(:,3),seshv4(iSesh).spkMovUCOsd(:,3)/sqrt(seshv4(iSesh).nRoi),'cmap',[.8 .5 0],'alpha');
        set(h2,'linewidth',1.5);
        h3 = boundedline(tPlot,seshv4(iSesh).spkMovUCOmu(:,2),seshv4(iSesh).spkMovUCOsd(:,2)/sqrt(seshv4(iSesh).nRoi),'cmap',[0 0 0],'alpha');
        set(h3,'linewidth',1.5);
        axis([-1 1 -2 10]); axis square;
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
        title('Spks, movement')
        
        subplot(2,2,4); % spk, off
        h1 = boundedline(tPlot,seshv4(iSesh).spkOffUCOmu(:,1),seshv4(iSesh).spkOffUCOsd(:,1)/sqrt(seshv4(iSesh).nRoi),'cmap',[0 .5 0],'alpha');
        set(h1,'linewidth',1.5);
        h2 = boundedline(tPlot,seshv4(iSesh).spkOffUCOmu(:,3),seshv4(iSesh).spkOffUCOsd(:,3)/sqrt(seshv4(iSesh).nRoi),'cmap',[.8 .5 0],'alpha');
        set(h2,'linewidth',1.5);
        h3 = boundedline(tPlot,seshv4(iSesh).spkOffUCOmu(:,2),seshv4(iSesh).spkOffUCOsd(:,2)/sqrt(seshv4(iSesh).nRoi),'cmap',[0 0 0],'alpha');
        set(h3,'linewidth',1.5);
        axis([-1 1 -2 10]); axis square;
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
        title('Spks, offset')
        
        suptitle(sprintf('Mouse: %s, FOV: %s, Date: %s',seshv4(iSesh).name, seshv4(iSesh).fov, seshv4(iSesh).date));
        
        
        if savebool
            fName = sprintf('%s_%s_%s',seshv4(iSesh).fov, seshv4(iSesh).name, seshv4(iSesh).date);
            savefig(f1, fullfile(figFold,[fName,'.fig']));
            saveas(f1,fullfile(figFold,[fName,'.png']));
            print(f1,fullfile(figFold,[fName,'.eps']), '-depsc', '-painters');
        end
    end
end

%% Plot summary data
iLobv = cat(1,seshv4.iLobv);
iSim2 = ~iLobv;
zfMovUCO = cat(3,seshv4.zfMovUCOmu);
zfOffUCO = cat(3,seshv4.zfOffUCOmu);
spkMovUCO = cat(3,seshv4.spkMovUCOmu);
spkOffUCO = cat(3,seshv4.spkOffUCOmu);

% Lobule V
f2 = figure;
subplot(2,2,1); % zF, mov
h1 = boundedline(tPlot,mean(zfMovUCO(:,1,iLobv),3),std(zfMovUCO(:,1,iLobv),[],3)/sqrt(sum(iLobv)),'cmap',[0 .5 0],'alpha');
set(h1,'linewidth',1.5);
h2 = boundedline(tPlot,mean(zfMovUCO(:,3,iLobv),3),std(zfMovUCO(:,3,iLobv),[],3)/sqrt(sum(iLobv)),'cmap',[.8 .5 0],'alpha');
set(h2,'linewidth',1.5);
h3 = boundedline(tPlot,mean(zfMovUCO(:,2,iLobv),3),std(zfMovUCO(:,2,iLobv),[],3)/sqrt(sum(iLobv)),'cmap',[0 0 0],'alpha');
set(h3,'linewidth',1.5);
axis([-1 1 -.5 1.5]); axis square;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Fluo, movement')
xline(0,'--b')

subplot(2,2,2); % zF, off
h1 = boundedline(tPlot,mean(zfOffUCO(:,1,iLobv),3),std(zfOffUCO(:,1,iLobv),[],3)/sqrt(sum(iLobv)),'cmap',[0 .5 0],'alpha');
set(h1,'linewidth',1.5);
h2 = boundedline(tPlot,mean(zfOffUCO(:,3,iLobv),3),std(zfOffUCO(:,3,iLobv),[],3)/sqrt(sum(iLobv)),'cmap',[.8 .5 0],'alpha');
set(h2,'linewidth',1.5);
h3 = boundedline(tPlot,mean(zfOffUCO(:,2,iLobv),3),std(zfOffUCO(:,2,iLobv),[],3)/sqrt(sum(iLobv)),'cmap',[0 0 0],'alpha');
set(h3,'linewidth',1.5);
axis([-1 1 -.5 1.5]); axis square;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Fluo, offset')
xline(0,'--b')

subplot(2,2,3); % spk, mov
h1 = boundedline(tPlot,mean(spkMovUCO(:,1,iLobv),3),std(spkMovUCO(:,1,iLobv),[],3)/sqrt(sum(iLobv)),'cmap',[0 .5 0],'alpha');
set(h1,'linewidth',1.5);
h2 = boundedline(tPlot,mean(spkMovUCO(:,3,iLobv),3),std(spkMovUCO(:,3,iLobv),[],3)/sqrt(sum(iLobv)),'cmap',[.8 .5 0],'alpha');
set(h2,'linewidth',1.5);
h3 = boundedline(tPlot,mean(spkMovUCO(:,2,iLobv),3),std(spkMovUCO(:,2,iLobv),[],3)/sqrt(sum(iLobv)),'cmap',[0 0 0],'alpha');
set(h3,'linewidth',1.5);
axis([-1 1 0 10]); axis square;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Fluo, movement')
xline(0,'--b')

subplot(2,2,4); % spk, off
h1 = boundedline(tPlot,mean(spkOffUCO(:,1,iLobv),3),std(spkOffUCO(:,1,iLobv),[],3)/sqrt(sum(iLobv)),'cmap',[0 .5 0],'alpha');
set(h1,'linewidth',1.5);
h2 = boundedline(tPlot,mean(spkOffUCO(:,3,iLobv),3),std(spkOffUCO(:,3,iLobv),[],3)/sqrt(sum(iLobv)),'cmap',[.8 .5 0],'alpha');
set(h2,'linewidth',1.5);
h3 = boundedline(tPlot,mean(spkOffUCO(:,2,iLobv),3),std(spkOffUCO(:,2,iLobv),[],3)/sqrt(sum(iLobv)),'cmap',[0 0 0],'alpha');
set(h3,'linewidth',1.5);
axis([-1 1 0 10]); axis square;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Fluo, movement')
xline(0,'--b')

suptitle('Lobule V fields of view');

% Simplex
f3 = figure;
subplot(2,2,1); % zF, mov
h1 = boundedline(tPlot,mean(zfMovUCO(:,1,iSim2),3),std(zfMovUCO(:,1,iSim2),[],3)/sqrt(sum(iSim2)),'cmap',[0 .5 0],'alpha');
set(h1,'linewidth',1.5);
h2 = boundedline(tPlot,mean(zfMovUCO(:,3,iSim2),3),std(zfMovUCO(:,3,iSim2),[],3)/sqrt(sum(iSim2)),'cmap',[.8 .5 0],'alpha');
set(h2,'linewidth',1.5);
h3 = boundedline(tPlot,mean(zfMovUCO(:,2,iSim2),3),std(zfMovUCO(:,2,iSim2),[],3)/sqrt(sum(iSim2)),'cmap',[0 0 0],'alpha');
set(h3,'linewidth',1.5);
axis([-1 1 -.5 1.5]); axis square;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Fluo, movement')
xline(0,'--b')

subplot(2,2,2); % zF, off
h1 = boundedline(tPlot,mean(zfOffUCO(:,1,iSim2),3),std(zfOffUCO(:,1,iSim2),[],3)/sqrt(sum(iSim2)),'cmap',[0 .5 0],'alpha');
set(h1,'linewidth',1.5);
h2 = boundedline(tPlot,mean(zfOffUCO(:,3,iSim2),3),std(zfOffUCO(:,3,iSim2),[],3)/sqrt(sum(iSim2)),'cmap',[.8 .5 0],'alpha');
set(h2,'linewidth',1.5);
h3 = boundedline(tPlot,mean(zfOffUCO(:,2,iSim2),3),std(zfOffUCO(:,2,iSim2),[],3)/sqrt(sum(iSim2)),'cmap',[0 0 0],'alpha');
set(h3,'linewidth',1.5);
axis([-1 1 -.5 1.5]); axis square;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Fluo, offset')
xline(0,'--b')

subplot(2,2,3); % spk, mov
h1 = boundedline(tPlot,mean(spkMovUCO(:,1,iSim2),3),std(spkMovUCO(:,1,iSim2),[],3)/sqrt(sum(iSim2)),'cmap',[0 .5 0],'alpha');
set(h1,'linewidth',1.5);
h2 = boundedline(tPlot,mean(spkMovUCO(:,3,iSim2),3),std(spkMovUCO(:,3,iSim2),[],3)/sqrt(sum(iSim2)),'cmap',[.8 .5 0],'alpha');
set(h2,'linewidth',1.5);
h3 = boundedline(tPlot,mean(spkMovUCO(:,2,iSim2),3),std(spkMovUCO(:,2,iSim2),[],3)/sqrt(sum(iSim2)),'cmap',[0 0 0],'alpha');
set(h3,'linewidth',1.5);
axis([-1 1 0 10]); axis square;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Fluo, movement')
xline(0,'--b')

subplot(2,2,4); % spk, off
h1 = boundedline(tPlot,mean(spkOffUCO(:,1,iSim2),3),std(spkOffUCO(:,1,iSim2),[],3)/sqrt(sum(iSim2)),'cmap',[0 .5 0],'alpha');
set(h1,'linewidth',1.5);
h2 = boundedline(tPlot,mean(spkOffUCO(:,3,iSim2),3),std(spkOffUCO(:,3,iSim2),[],3)/sqrt(sum(iSim2)),'cmap',[.8 .5 0],'alpha');
set(h2,'linewidth',1.5);
h3 = boundedline(tPlot,mean(spkOffUCO(:,2,iSim2),3),std(spkOffUCO(:,2,iSim2),[],3)/sqrt(sum(iSim2)),'cmap',[0 0 0],'alpha');
set(h3,'linewidth',1.5);
axis([-1 1 0 10]); axis square;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Fluo, movement')
xline(0,'--b')

suptitle('Lobule simplex fields of view');

% Save maybe
if savebool
    fName = '01_lobvUCO_summary';
    savefig(f2, fullfile(figFold,[fName,'.fig']));
    saveas(f2,fullfile(figFold,[fName,'.png']));
    print(f2,fullfile(figFold,[fName,'.eps']), '-depsc', '-painters');
    
    fName = '02_sim2UCO_summary';
    savefig(f3, fullfile(figFold,[fName,'.fig']));
    saveas(f3,fullfile(figFold,[fName,'.png']));
    print(f3,fullfile(figFold,[fName,'.eps']), '-depsc', '-painters');
end
