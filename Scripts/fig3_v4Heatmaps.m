function iV4n = fig3_v4Heatmaps(matchData, matchFOVs, dataFold, paperFold, savebool)
%% Get expert sessions:
% Use mov2 and rewT2 data (which should always have V4n data) to make new
% structures that have what we need:
figFold = fullfile(paperFold,'Figures','Fig3','Sub1_behavSummary');
figFold2 = fullfile(dataFold,'paperFigs','Fig3','Sub1_behavSummary');
imfs = 30;
vrfs = 100;
rewDelay = 0.4;


matchV4 = struct;
dataNames = cellstr(cat(1,matchData.name));
dataFOVs = cellstr(cat(1,matchData.fov));
dataDates1 = cellstr(cat(1,matchData.date1));
dataDates2 = cellstr(cat(1,matchData.date2));
ilobv = false(size(matchFOVs));

iV4n = zeros(numel(matchFOVs),2);

for iFOVs = 1:numel(matchFOVs)
    
    iNameMatch = cellfun( @(x) strcmpi(x,matchFOVs(iFOVs).name),dataNames);
    iFOVMatch = cellfun( @(x) strcmpi(x,matchFOVs(iFOVs).fov),dataFOVs);    
    % First try with date1:
    iDate1Match = cellfun( @(x) strcmpi(x,matchFOVs(iFOVs).v4_Ndate),dataDates1);
    iDate1 = find(all([iNameMatch iFOVMatch iDate1Match],2),1);
    
    if ~isempty(iDate1)
        iMov = 'mov1';
        iRew = 'rewT1';
        iDate = iDate1;
        iV4n(iFOVs,2) = 1;
    else
        % If date1 doesn't have any matches, use date2
        iDate2Match = cellfun( @(x) strcmpi(x,matchFOVs(iFOVs).v4_Ndate),dataDates2);
        iDate = find(all([iNameMatch iFOVMatch iDate2Match],2),1); % no need to name, it is the only one left
        iMov = 'mov2';
        iRew = 'rewT2';
        iV4n(iFOVs,2) = 2;
    end
    matchV4(iFOVs).name = matchFOVs(iFOVs).name;
    matchV4(iFOVs).date = matchFOVs(iFOVs).v4_Ndate;
    ilobv(iFOVs) = matchFOVs(iFOVs).ilobv;
    matchV4(iFOVs).mov = matchData(iDate).(iMov);
    matchV4(iFOVs).rew = matchData(iDate).(iRew);
    iV4n(iFOVs,1) = iDate;
end

%% Average across trials so that we can concatenate
matchLobv = matchV4(ilobv);
for i = 1:numel(matchLobv)
    matchLobv(i).movSpk = imfs*squeeze(mean(matchLobv(i).mov.spk(end/2-imfs+1:end/2+imfs,:,:),3,'omitnan'));
    matchLobv(i).movMov = mean(matchLobv(i).mov.mov(end/2-vrfs+1:end/2+vrfs,:),2,'omitnan');
    matchLobv(i).rewSpk = imfs*squeeze(mean(matchLobv(i).rew.spk(rewDelay*imfs+[end/2-imfs+1:end/2+imfs],:,:),3,'omitnan'));
    matchLobv(i).rewMov = mean(matchLobv(i).rew.mov(rewDelay*vrfs+[end/2-vrfs+1:end/2+vrfs],:),2,'omitnan');
end

matchSim2 = matchV4(~ilobv);
for i = 1:numel(matchSim2)
    matchSim2(i).movSpk = imfs*squeeze(mean(matchSim2(i).mov.spk(end/2-imfs+1:end/2+imfs,:,:),3,'omitnan'));
    matchSim2(i).movMov = mean(matchSim2(i).mov.mov(end/2-vrfs+1:end/2+vrfs,:),2,'omitnan');    
    matchSim2(i).rewSpk = imfs*squeeze(mean(matchSim2(i).rew.spk(rewDelay*imfs+[end/2-imfs+1:end/2+imfs],:,:),3,'omitnan'));
    matchSim2(i).rewMov = mean(matchSim2(i).rew.mov(rewDelay*vrfs+[end/2-vrfs+1:end/2+vrfs],:),2,'omitnan');
end

%% Make plotting matrices
% lobv mov
lobVmovSpk = cat(2,matchLobv.movSpk);
lobVmovMov = cat(2,matchLobv.movMov);
coeff = pca(lobVmovSpk(end/2-imfs/2+1:end/2+imfs/2,:));
[~, iSort] = sort(coeff(:,1),'descend');
lobVmovSpk = lobVmovSpk(:,iSort);
lobVmovSpkMu = mean(lobVmovSpk,2);
lobVmovSpkSem = std(lobVmovSpk,[],2)/sqrt(size(lobVmovSpk,2));
lobVmovMovMu = mean(lobVmovMov,2);
lobVmovMovSem = std(lobVmovMov,[],2)/sqrt(size(lobVmovMov,2));

% lobV rew
lobVrewSpk = cat(2,matchLobv.rewSpk);
lobVrewMov = cat(2,matchLobv.rewMov);
coeff = pca(lobVrewSpk(end/2-imfs/2+1:end/2+imfs/2,:));
[~, iSort] = sort(coeff(:,1),'descend');
lobVrewSpk = lobVrewSpk(:,iSort);
lobVrewSpkMu = mean(lobVrewSpk,2);
lobVrewSpkSem = std(lobVrewSpk,[],2)/sqrt(size(lobVrewSpk,2));
lobVrewMovMu = mean(lobVrewMov,2);
lobVrewMovSem = std(lobVrewMov,[],2)/sqrt(size(lobVrewMov,2));

% sim2 mov
sim2movSpk = cat(2,matchSim2.movSpk);
sim2movMov = cat(2,matchSim2.movMov);
coeff = pca(sim2movSpk(end/2-imfs/2+1:end/2+imfs/2,:));
[~, iSort] = sort(coeff(:,1),'descend');
sim2movSpk = sim2movSpk(:,iSort);
sim2movSpkMu = mean(sim2movSpk,2);
sim2movSpkSem = std(sim2movSpk,[],2)/sqrt(size(sim2movSpk,2));
sim2movMovMu = mean(sim2movMov,2);
sim2movMovSem = std(sim2movMov,[],2)/sqrt(size(sim2movMov,2));

% sim2 rew
sim2rewSpk = cat(2,matchSim2.rewSpk);
sim2rewMov = cat(2,matchSim2.rewMov);
coeff = pca(sim2rewSpk(end/2-imfs/2+1:end/2+imfs/2,:));
[~, iSort] = sort(coeff(:,1),'descend');
sim2rewSpk = sim2rewSpk(:,iSort);
sim2rewSpkMu = mean(sim2rewSpk,2);
sim2rewSpkSem = std(sim2rewSpk,[],2)/sqrt(size(sim2rewSpk,2));
sim2rewMovMu = mean(sim2rewMov,2);
sim2rewMovSem = std(sim2rewMov,[],2)/sqrt(size(sim2rewMov,2));

%% Plot

f1 = figure;
subplot(2,4,1);
imagesc(lobVmovSpk',[0 8]);
colormap viridis; box off
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([0.5:imfs/2:2*imfs]);
xline(imfs+.5,'--w','LineWidth',1);
xticklabels({});
ylabel('LobV - mov');

subplot(2,4,2);
imagesc(lobVrewSpk',[0 8]);
colormap viridis; box off
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([0.5:imfs/2:2*imfs]);
xline(imfs+.5,'--w','LineWidth',1);
xticklabels({});yticklabels({});
ylabel('LobV - reward');

subplot(2,4,3);
imagesc(sim2movSpk',[0 8]);
colormap viridis; box off
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([0.5:imfs/2:2*imfs]);
xline(imfs+.5,'--w','LineWidth',1);
xticklabels({});
ylabel('Sim2 - mov');

subplot(2,4,4);
imagesc(sim2rewSpk',[0 8]);
colormap viridis; box off
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticks([0.5:imfs/2:2*imfs]);
xline(imfs+.5,'--w','LineWidth',1);
xticklabels({});yticklabels({});
ylabel('Sim2 - reward');

tSpk = (.5:1:59.5)';
subplot(4,4,9);
h = boundedline(tSpk,lobVmovSpkMu,lobVmovSpkSem,'cmap',[0 0 0],'alpha');
set(h,'linewidth',2);
ylim([0 6]); yticks(0:1:6); xticks(0:imfs/2:2*imfs);
xline(imfs,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticklabels({'-1','','0','','1'});yticklabels({'0','','2','','4','','6'})

subplot(4,4,10);
h = boundedline(tSpk,lobVrewSpkMu,lobVrewSpkSem,'cmap',[0 0 0],'alpha');
set(h,'linewidth',2);
ylim([0 6]); yticks(0:1:6); xticks(0:imfs/2:2*imfs);
xline(imfs,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticklabels({});yticklabels({});

subplot(4,4,11);
h = boundedline(tSpk,sim2movSpkMu,sim2movSpkSem,'cmap',[0 0 0],'alpha');
set(h,'linewidth',2);
ylim([0 6]); yticks(0:1:6); xticks(0:imfs/2:2*imfs);
xline(imfs,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticklabels({});yticklabels({});

subplot(4,4,12);
h = boundedline(tSpk,sim2rewSpkMu,sim2rewSpkSem,'cmap',[0 0 0],'alpha');
set(h,'linewidth',2);
ylim([0 6]); yticks(0:1:6); xticks(0:imfs/2:2*imfs);
xline(imfs,'--k','LineWidth',1);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticklabels({});yticklabels({});

tMov = (.5:1:2*vrfs)';
subplot(4,4,13);
h = boundedline(tMov,lobVmovMovMu,lobVmovMovSem,'cmap',[0 0 0],'alpha');
set(h,'linewidth',2);
xline(vrfs,'--k','LineWidth',1);
ylim([-2.5 12.5]); yticks(-2.5:2.5:10); xticks(0:vrfs/4:2*vrfs);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticklabels({'-1','','','','0','','','','1'});yticklabels({'','0','','5','','10',''});

subplot(4,4,14);
h = boundedline(tMov,lobVrewMovMu,lobVrewMovSem,'cmap',[0 0 0],'alpha');
set(h,'linewidth',2);
xline(vrfs,'--k','LineWidth',1);
ylim([-2.5 12.5]); yticks(-2.5:2.5:10); xticks(0:vrfs/4:2*vrfs);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticklabels({});yticklabels({});

subplot(4,4,15);
h = boundedline(tMov,sim2movMovMu,sim2movMovSem,'cmap',[0 0 0],'alpha');
set(h,'linewidth',2);
xline(vrfs,'--k','LineWidth',1);
ylim([-2.5 12.5]); yticks(-2.5:2.5:10); xticks(0:vrfs/4:2*vrfs);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticklabels({});yticklabels({});

subplot(4,4,16);
h = boundedline(tMov,sim2rewMovMu,sim2rewMovSem,'cmap',[0 0 0],'alpha');
set(h,'linewidth',2);
xline(vrfs,'--k','LineWidth',1);
ylim([-2.5 12.5]); yticks(-2.5:2.5:10); xticks(0:vrfs/4:2*vrfs);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xticklabels({});yticklabels({});

suptitle('Expert heatmaps')

if savebool
    if ~exist(figFold,'dir'); mkdir(figFold); end
    if ~exist(figFold2,'dir'); mkdir(figFold2); end
    figname = 'v4Heatmaps';
    savefig(f1, fullfile(figFold,[figname,'.fig']));
    saveas(f1,fullfile(figFold,[figname,'.png']));
    print(f1,fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');    
end


end