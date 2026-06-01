function fig2_exampleFig(dataPaths,dataFold,paperFold,savebool)
lobVpath0 = fullfile(dataFold,sprintf('Version_%i',dataPaths(2).version),dataPaths(2).name,...
    sprintf('%s_%s',dataPaths(2).date,dataPaths(2).fov));
lobVpath = fullfile(lobVpath0,'GLM');
sim2path0 = fullfile(dataFold,sprintf('Version_%i',dataPaths(1).version),dataPaths(1).name,...
    sprintf('%s_%s',dataPaths(1).date,dataPaths(1).fov));
sim2path = fullfile(sim2path0,'GLM');

figFold = fullfile(paperFold,'Figures','Fig2','Sub1_BehavExample');

%% DK169 - lobV:

% lobVdm = load(fullfile(lobVpath,'DM+spikes_cells.mat'));
lobVglm = load(fullfile(lobVpath,'glmRfiles','GLMfull_small.mat')); lobVglm = lobVglm.GLMfull_small;
lobVbasis = load(fullfile(lobVpath,'basisstructs_cells.mat'));
expt = lobVbasis.trialbasis.expt;

% Trim out bad trials from expt!
trialstruct = load(fullfile(lobVpath0,'trialsonly.mat'));
trialsUCO = sort(cat(1,trialstruct.trialstructs.LUstruct.trialind,trialstruct.trialstructs.LCstruct.trialind,trialstruct.trialstructs.LOstruct.trialind));
trialsALL = cat(1,expt.trialnum);    
expt = expt(ismember(trialsALL,trialsUCO));

exTrials = [75:84]; % Real trial num 89-98
nTrials = numel(exTrials);

xJump = 0;
trimStart = 24;
trimEnd = 30;
trialBreak = 6;

tOnOff = [];
tRew = [];
tLick = [];
wheelPos = [];
wheelVel = [];
movOns = nan(nTrials,1);
rewOns = movOns;
spkData = [];
glmData = [];

for iTrial = 1:nTrials
    % Get single trial structure
    exptLocal = expt(exTrials(iTrial));    
    offset = exptLocal.start-1;
    exptLocal.trialon = exptLocal.trialon-offset;
    exptLocal.trialoff = exptLocal.trialoff-offset;
    exptLocal.reward = exptLocal.reward-offset;
    exptLocal.lick = exptLocal.lick-offset;
    exptLocal.movonset = exptLocal.movonset-offset;
    
    % Extract local vectors, etc.
    xPts = (1:exptLocal.duration-trimStart-trimEnd)'+xJump;
    tOnOffLocal = zeros(exptLocal.duration,1); tRewLocal = tOnOffLocal; tLickLocal = tOnOffLocal;
    tOnOffLocal(exptLocal.trialon:exptLocal.trialoff) = 1; tOnOffLocal = tOnOffLocal(trimStart+1:end-trimEnd);
    tRewLocal(exptLocal.reward) = 1; tRewLocal = tRewLocal(trimStart+1:end-trimEnd);
    tLickLocal(exptLocal.lick) = 1; tLickLocal = tLickLocal(trimStart+1:end-trimEnd);
    wheelPosLocal = exptLocal.wheelpos; wheelPosLocal = wheelPosLocal(trimStart+1:end-trimEnd);
    wheelPosLocal(end-30:end) = wheelPosLocal(end-31);
    wheelVelLocal = exptLocal.wheelvel/-10; wheelVelLocal = wheelVelLocal(trimStart+1:end-trimEnd);
    spkDataLocal = exptLocal.spkData(trimStart+1:end-trimEnd,:);
    glmDataLocal = lobVglm.DM(lobVglm.DMdata.endTrialIndices(exTrials(iTrial))+1:lobVglm.DMdata.endTrialIndices(exTrials(iTrial)+1),:);
    
    % Populate big vector
    tOnOff(xPts,1) = tOnOffLocal; tOnOff(xPts(end)+1:xPts(end)+ trialBreak) = nan;
    tRew(xPts,1) = tRewLocal; tRew(xPts(end)+1:xPts(end)+ trialBreak) = nan;
    tLick(xPts,1) = tLickLocal; tLick(xPts(end)+1:xPts(end)+ trialBreak) = nan;
    wheelPos(xPts,1) = wheelPosLocal; wheelPos(xPts(end)+1:xPts(end)+ trialBreak) = nan;
    wheelVel(xPts,1) = wheelVelLocal; wheelVel(xPts(end)+1:xPts(end)+ trialBreak) = nan;
    if ~isempty(exptLocal.movonset); movOns(iTrial) = exptLocal.movonset-trimStart+xJump; end
    if ~isempty(exptLocal.reward); rewOns(iTrial) = exptLocal.reward-trimStart+xJump; end
    xJump = numel(tOnOff);
    spkData = [spkData; spkDataLocal]; spkData(end+1:end+trialBreak,:) = nan;
    glmData = [glmData; glmDataLocal]; glmData(end+1:end+trialBreak,:) = nan;
end

xPts = 1:numel(tOnOff);
iRews = threshdetect_dk(tRew,.5,1);
iLicks = threshdetect_dk(tLick,.5,1);
iRewT = [zeros(1,numel(iRews)); ones(1,numel(iRews))]/2+.5;
iLickT = [zeros(1,numel(iLicks)); ones(1,numel(iLicks))]/2;

fLobV = figure;
subplot(3,1,1);
stairs(xPts,tOnOff+6.5,'-k','LineWidth',2);
hold on;
rectangle('Position',[1,3+2/3,xPts(end),1],'FaceColor',[.6 .0 .0],'EdgeColor','none');
rectangle('Position',[1,4+2/3,xPts(end),2/3],'FaceColor',[.8 .8 .8],'EdgeColor','none');
rectangle('Position',[1,5+1/3,xPts(end),1],'FaceColor',[.0 .6 .6],'EdgeColor','none');
plot(xPts,wheelPos+4,'-k','LineWidth',1.5);
plot(xPts,wheelVel+2,'-k','LineWidth',1.5);
plot([iRews';iRews'],iRewT,'-c','LineWidth',4);
plot([iLicks';iLicks'],iLickT,'-g','LineWidth',4);
% stairs(xPts,tRew,'-b','LineWidth',1.5);
% stairs(xPts,tLick,'-r','LineWidth',1.5);
title(sprintf('Lobule V, Trial range: %i-%i',exTrials(1),exTrials(end)));

for iMov = 1:numel(movOns)
    if ~isnan(movOns(iMov)); xline(movOns(iMov),'--','LineWidth',1.5); end
end
for iRew = 1:numel(rewOns)
    if ~isnan(rewOns(iRew)); xline(rewOns(iRew),'--','LineWidth',1.5); end
end
ylim([-.5 8]); box off; axis off;
xlim([.5 xPts(end)+.5]);
%     axis image;

hSpk = subplot(3,1,2);
imagesc(spkData',[-.0 2/3]); colormap(hSpk,viridis);
hold on; box off
for iMov = 1:numel(movOns)
    if ~isnan(movOns(iMov)); xline(movOns(iMov),'--w','LineWidth',1.5); end
    if ~isnan(rewOns(iMov)); xline(rewOns(iMov),'--w','LineWidth',1.5); end
end
xticks(.5:60:xPts(end)); xticklabels({});
yticks([.5:50:size(spkData,2)]); yticklabels({});
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');

hGLM = subplot(3,1,3);
imagesc(glmData',[-.1 1]); colormap(hGLM,flipud(gray));
hold on; box off
for iBases = 1:numel(lobVglm.DMdata.nBases)
    yline(cumsum(lobVglm.DMdata.nBases(1:iBases)),'--b','LineWidth',1);
end
yticks(cumsum(lobVglm.DMdata.nBases(1:iBases))); yticklabels({});
xticks(.5:60:xPts(end)); xticklabels({});
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');

if savebool
    figname = sprintf('%s_LobV_Trials_%i_to_%i',dataPaths(2).name,exTrials(1),exTrials(end));
    savefig(fLobV, fullfile(figFold,[figname,'.fig']));
    saveas(fLobV,fullfile(figFold,[figname,'.png']));
    print(fLobV,fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
end

%% Remake predictor matrix + DM correlation
fGLMcorr = figure;
hGLM = subplot(1,2,1);
imagesc(glmData',[-.1 1]); colormap(hGLM,flipud(gray));
hold on; box off
for iBases = 1:numel(lobVglm.DMdata.nBases)
    yline(cumsum(lobVglm.DMdata.nBases(1:iBases)),'--b','LineWidth',1);
end
yticks(cumsum(lobVglm.DMdata.nBases(1:iBases))); yticklabels({});
xticks(.5:60:xPts(end)); xticklabels({}); axis square;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');

Rdm = corrcoef(lobVglm.DM);
hCorr = subplot(1,2,2);
imagesc(Rdm,[-1 1]); colormap(hCorr,redblue);
hold on; box off; axis square;
for iBases = 1:numel(lobVglm.DMdata.nBases)
    xline(cumsum(lobVglm.DMdata.nBases(1:iBases)),'--b','LineWidth',1);
    yline(cumsum(lobVglm.DMdata.nBases(1:iBases)),'--b','LineWidth',1);
end
xticks(cumsum(lobVglm.DMdata.nBases(1:iBases))); xticklabels({});
yticks(cumsum(lobVglm.DMdata.nBases(1:iBases))); yticklabels({});
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');

if savebool
    figname = sprintf('%s_LobV_glmCorrs',dataPaths(2).name);
    savefig(fGLMcorr, fullfile(figFold,[figname,'.fig']));
    saveas(fGLMcorr,fullfile(figFold,[figname,'.png']));
    print(fGLMcorr,fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
end


%% DK169 - sim2:
% sim2dm = load(fullfile(sim2path,'DM+spikes_cells.mat'));
sim2glm = load(fullfile(sim2path,'glmRfiles','GLMfull_small.mat')); sim2glm = sim2glm.GLMfull_small;
sim2basis = load(fullfile(sim2path,'basisstructs_cells.mat'));
expt = sim2basis.trialbasis.expt;
exTrials = [141:150]; % Real trials 145:155 (missing 153)
nTrials = numel(exTrials);

% Trim out bad trials from expt!
trialstruct = load(fullfile(sim2path0,'trialsonly.mat'));
trialsUCO = sort(cat(1,trialstruct.trialstructs.LUstruct.trialind,trialstruct.trialstructs.LCstruct.trialind,trialstruct.trialstructs.LOstruct.trialind));
trialsALL = cat(1,expt.trialnum);    
expt = expt(ismember(trialsALL,trialsUCO));

xJump = 0;
trimStart = 24;
trimEnd = 30;
trialBreak = 6;

tOnOff = [];
tRew = [];
tLick = [];
wheelPos = [];
wheelVel = [];
movOns = nan(nTrials,1);
rewOns = movOns;
spkData = [];
glmData = [];

for iTrial = 1:nTrials
    % Get single trial structure
    exptLocal = expt(exTrials(iTrial));    
    offset = exptLocal.start-1;
    exptLocal.trialon = exptLocal.trialon-offset;
    exptLocal.trialoff = exptLocal.trialoff-offset;
    exptLocal.reward = exptLocal.reward-offset;
    exptLocal.lick = exptLocal.lick-offset;
    exptLocal.movonset = exptLocal.movonset-offset;
    
    % Extract local vectors, etc.
    xPts = (1:exptLocal.duration-trimStart-trimEnd)'+xJump;
    tOnOffLocal = zeros(exptLocal.duration,1); tRewLocal = tOnOffLocal; tLickLocal = tOnOffLocal;
    tOnOffLocal(exptLocal.trialon:exptLocal.trialoff) = 1; tOnOffLocal = tOnOffLocal(trimStart+1:end-trimEnd);
    tRewLocal(exptLocal.reward) = 1; tRewLocal = tRewLocal(trimStart+1:end-trimEnd);
    tLickLocal(exptLocal.lick) = 1; tLickLocal = tLickLocal(trimStart+1:end-trimEnd);
    wheelPosLocal = exptLocal.wheelpos; wheelPosLocal = wheelPosLocal(trimStart+1:end-trimEnd);
    wheelPosLocal(end-30:end) = wheelPosLocal(end-31);
    wheelVelLocal = exptLocal.wheelvel/-10; wheelVelLocal = wheelVelLocal(trimStart+1:end-trimEnd);
    spkDataLocal = exptLocal.spkData(trimStart+1:end-trimEnd,:);
    glmDataLocal = sim2glm.DM(sim2glm.DMdata.endTrialIndices(exTrials(iTrial))+1:sim2glm.DMdata.endTrialIndices(exTrials(iTrial)+1),:);
    
    % Populate big vector
    tOnOff(xPts,1) = tOnOffLocal; tOnOff(xPts(end)+1:xPts(end)+ trialBreak) = nan;
    tRew(xPts,1) = tRewLocal; tRew(xPts(end)+1:xPts(end)+ trialBreak) = nan;
    tLick(xPts,1) = tLickLocal; tLick(xPts(end)+1:xPts(end)+ trialBreak) = nan;
    wheelPos(xPts,1) = wheelPosLocal; wheelPos(xPts(end)+1:xPts(end)+ trialBreak) = nan;
    wheelVel(xPts,1) = wheelVelLocal; wheelVel(xPts(end)+1:xPts(end)+ trialBreak) = nan;
    if ~isempty(exptLocal.movonset); movOns(iTrial) = exptLocal.movonset-trimStart+xJump; end
    if ~isempty(exptLocal.reward); rewOns(iTrial) = exptLocal.reward-trimStart+xJump; end
    xJump = numel(tOnOff);
    spkData = [spkData; spkDataLocal]; spkData(end+1:end+trialBreak,:) = nan;
    glmData = [glmData; glmDataLocal]; glmData(end+1:end+trialBreak,:) = nan;    
end

xPts = 1:numel(tOnOff);
iRews = threshdetect_dk(tRew,.5,1);
iLicks = threshdetect_dk(tLick,.5,1);
iRewT = [zeros(1,numel(iRews)); ones(1,numel(iRews))]/2+.5;
iLickT = [zeros(1,numel(iLicks)); ones(1,numel(iLicks))]/2;

fsim2 = figure;
subplot(3,1,1);
stairs(xPts,tOnOff+6.5,'-k','LineWidth',2);
hold on;
rectangle('Position',[1,3+2/3,xPts(end),1],'FaceColor',[.6 .0 .0],'EdgeColor','none');
rectangle('Position',[1,4+2/3,xPts(end),2/3],'FaceColor',[.8 .8 .8],'EdgeColor','none');
rectangle('Position',[1,5+1/3,xPts(end),1],'FaceColor',[.0 .6 .6],'EdgeColor','none');
plot(xPts,wheelPos+4,'-k','LineWidth',1.5);
plot(xPts,wheelVel+2,'-k','LineWidth',1.5); 
plot([iRews';iRews'],iRewT,'-c','LineWidth',4);
plot([iLicks';iLicks'],iLickT,'-g','LineWidth',4);
% stairs(xPts,tRew,'-b','LineWidth',1.5);
% stairs(xPts,tLick,'-r','LineWidth',1.5);
title(sprintf('Lobule simplex, Trial range: %i-%i',exTrials(1),exTrials(end)));

for iMov = 1:numel(movOns)
    if ~isnan(movOns(iMov)); xline(movOns(iMov),'--','LineWidth',1.5); end
end
for iRew = 1:numel(rewOns)
    if ~isnan(rewOns(iRew)); xline(rewOns(iRew),'--','LineWidth',1.5); end
end
ylim([-.5 8]); box off; axis off;
xlim([.5 xPts(end)+.5]);
%     axis image;

hSpk = subplot(3,1,2);
imagesc(spkData',[0 2/3]); colormap(hSpk,viridis);
hold on; box off
for iMov = 1:numel(movOns)
    if ~isnan(movOns(iMov)); xline(movOns(iMov),'--w','LineWidth',1.5); end
    if ~isnan(rewOns(iMov)); xline(rewOns(iMov),'--w','LineWidth',1.5); end
end
xticks(.5:60:xPts(end)); xticklabels({});
yticks([.5:50:size(spkData,2)]); yticklabels({});
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');

hGLM = subplot(3,1,3);
imagesc(glmData',[-.1 1]); colormap(hGLM,flipud(gray));
hold on; box off
for iBases = 1:numel(lobVglm.DMdata.nBases)
    yline(cumsum(lobVglm.DMdata.nBases(1:iBases)),'--b','LineWidth',1);
end
yticks(cumsum(lobVglm.DMdata.nBases(1:iBases))); yticklabels({});
xticks(.5:60:xPts(end)); xticklabels({});
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');

% subplot(3,1,3);
% spkData2 = spkData; spkData2(spkData2>0) = 1;
% imagesc(spkData2',[-.2 1]); colormap viridis;
% hold on; box off
% for iMov = 1:numel(movOns)
%     if ~isnan(movOns(iMov)); xline(movOns(iMov),'--w','LineWidth',1.5); end
%     if ~isnan(rewOns(iMov)); xline(rewOns(iMov),'--w','LineWidth',1.5); end
% end
% xticks(.5:60:xPts(end)); xticklabels({});
% yticks([.5:50:size(spkData,2)]); yticklabels({});
% set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');

if savebool
    figname = sprintf('%s_sim2_Trials_%i_to_%i',dataPaths(2).name,exTrials(1),exTrials(end));
    savefig(fsim2, fullfile(figFold,[figname,'.fig']));
    saveas(fsim2,fullfile(figFold,[figname,'.png']));
    print(fsim2,fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
end

end