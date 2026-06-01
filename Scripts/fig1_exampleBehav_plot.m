function figOut = fig1_exampleBehav_plot(glmStruct,exptData,exTrials,mouseData)

%%
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
    exptLocal = exptData(exTrials(iTrial));    
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
    glmDataLocal = glmStruct.DM(glmStruct.DMdata.endTrialIndices(exTrials(iTrial))+1:glmStruct.DMdata.endTrialIndices(exTrials(iTrial)+1),:);
    
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

figOut = figure;
% subplot(3,1,1);
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
title(sprintf('%s: Version %i on %s, Trial range: %i-%i',mouseData.name,mouseData.version,mouseData.date,exTrials(1),exTrials(end)));

for iMov = 1:numel(movOns)
    if ~isnan(movOns(iMov)); xline(movOns(iMov),'--','LineWidth',1.5); end
end
for iRew = 1:numel(rewOns)
    if ~isnan(rewOns(iRew)); xline(rewOns(iRew),'--','LineWidth',1.5); end
end
ylim([-.5 8]); box off;
xlim([0 720]); xticks(0:60:720);

ylabel('Object/wheel endpoint');
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');


end
