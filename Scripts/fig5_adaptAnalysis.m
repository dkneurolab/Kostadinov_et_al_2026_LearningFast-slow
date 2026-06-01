function fig5_adaptAnalysis(matchData, glmParams, dataFold, paperFold, savebool)
%% Settings:
% Decimation:
nDec = 50;
% Behav smoothing:
msSmthing = 60;% 60 ms smoothing span
% Remove redundant fields in matchData:
matchData = rmfield(matchData,{'on1','on2','rewC1','rewC2','rewR1','rewR2'});

% Load in v4 sigSummary
glmFold0 = fullfile(paperFold,'Figures','Fig3','Sub3_SummaryBars');
glmFold = fullfile(glmFold0,sprintf('%s_%s_alpha%i',glmParams.dataType,glmParams.distr,glmParams.alpha*100));
glmSummary = load(fullfile(glmFold,'sigSummary.mat')); glmSummary = glmSummary.glmSummary;

% 
pathv5Behav = fullfile(paperFold,'Figures','Fig1','Sub1_percMovCorr');
sessionsv5 = load(fullfile(pathv5Behav,'percMovCorrect_v5.mat')); sessionsv5 = sessionsv5.sessionsv5;

%% Trim
glmGood = false(numel(glmSummary),1);
v5Good = false(numel(sessionsv5),1);
v5lobv = v5Good;
glmlobv = cat(1,glmSummary.iLobV);
% Find matching good sessions
for iGlm = 1:numel(glmGood)    
    for jV5 = 1:numel(v5Good)
        if strcmpi(glmSummary(iGlm).name, sessionsv5(jV5).name) && strcmpi(glmSummary(iGlm).fov, sessionsv5(jV5).fov)
            glmGood(iGlm) = true;
            v5Good(jV5) = true;
            if strcmpi(sessionsv5(jV5).fov,'lobv')
                v5lobv(jV5) = true;
            end
        end
    end
end
% Extract them
seshv5lobv = sessionsv5(v5Good & v5lobv);
seshv5sim2 = sessionsv5(v5Good & ~v5lobv);
glmv5lobv0 = glmSummary(glmGood & glmlobv);
glmv5sim20 = glmSummary(glmGood & ~glmlobv);

% Find v4-v5 match sessions:
iMatchlobv = false(numel(matchData),1);
iMatchsim2 = false(numel(matchData),1);
for iMatch = 1:numel(matchData)
    if matchData(iMatch).v1 == 4 && matchData(iMatch).v2 == 5
        if strcmpi(matchData(iMatch).fov,'lobv')
           iMatchlobv(iMatch) = true; 
        elseif strcmpi(matchData(iMatch).fov,'sim2')
            iMatchsim2(iMatch) = true;
        end
    end
end
% Extract them
matchV5lobv = matchData(iMatchlobv);
matchV5sim2 = matchData(iMatchsim2);

clear glmSummary sessionsv5 pathv5Behav glmFold0 glmGood glmlobv iGlm iMatch iMatchlobv iMatchsim2 jV5 v5Good v5lobv

%% Go through and dig out correct trials
dataPath0 = fullfile(dataFold,'01_Behav+imaging','Version_5');

% First LobV
glmv5lobv = struct;
for iFOV = 1:numel(glmv5lobv0)
    TS = load(fullfile(dataPath0,seshv5lobv(iFOV).name,[seshv5lobv(iFOV).date,'_',seshv5lobv(iFOV).fov],'trialsonly.mat'));
    roisv4 = matchV5lobv(iFOV).iROIs1;
    roisv5 = matchV5lobv(iFOV).iROIs2;
    
    % Take only good ROIs
    glmv5lobv(iFOV).iROIv4 = roisv4;
    glmv5lobv(iFOV).iROIv5 = roisv5;
    glmv5lobv(iFOV).sigBool = glmv5lobv0(iFOV).sigBool(roisv4);
    glmv5lobv(iFOV).sigBoolVis = glmv5lobv0(iFOV).sigBoolVis(roisv4);
    glmv5lobv(iFOV).sigBoolMov = glmv5lobv0(iFOV).sigBoolMov(roisv4);
    glmv5lobv(iFOV).sigBoolRew = glmv5lobv0(iFOV).sigBoolRew(roisv4);
    glmv5lobv(iFOV).sigBoolLick = glmv5lobv0(iFOV).sigBoolLick(roisv4);
    glmv5lobv(iFOV).gainUp = seshv5lobv(iFOV).sets.changetrial;
    glmv5lobv(iFOV).gainDown = seshv5lobv(iFOV).sets.changebacktrial;
    
    glmv5lobv(iFOV).iTrialsMov = cat(1,TS.trialstructs.LCstruct.trialind,TS.trialstructs.LUstruct.trialind,TS.trialstructs.LOstruct.trialind); 
    spkMov = cat(3,TS.trialstructs.LCstruct.spkLmov,TS.trialstructs.LUstruct.spkLmov,TS.trialstructs.LOstruct.spkLmov);
    xMov = cat(3,TS.trialstructs.LCstruct.VRLmov,TS.trialstructs.LUstruct.VRLmov,TS.trialstructs.LOstruct.VRLmov); xMov = squeeze(xMov(:,3,:));
    [~,iMovSort] = sort(glmv5lobv(iFOV).iTrialsMov,'ascend');
    glmv5lobv(iFOV).iTrialsMov = glmv5lobv(iFOV).iTrialsMov(iMovSort);
    spkMov = spkMov(:,:,iMovSort); xMov = xMov(:,iMovSort);
    
    glmv5lobv(iFOV).iTrialsRew = cat(1,TS.trialstructs.LCstruct.trialind);
    glmv5lobv(iFOV).iTrialsUnder = cat(1,TS.trialstructs.LUstruct.trialind);
    glmv5lobv(iFOV).iTrialsOver = cat(1,TS.trialstructs.LOstruct.trialind);
    spkRew = cat(3,TS.trialstructs.LCstruct.spkLoff);
    spkUnder = cat(3,TS.trialstructs.LUstruct.spkLoff);
    spkOver = cat(3,TS.trialstructs.LOstruct.spkLoff);
    xRew = cat(3,TS.trialstructs.LCstruct.VRLoff); xRew = squeeze(xRew(:,3,:));
    xUnder = cat(3,TS.trialstructs.LUstruct.VRLoff); xUnder = squeeze(xUnder(:,3,:));
    xOver = cat(3,TS.trialstructs.LOstruct.VRLoff); xOver = squeeze(xOver(:,3,:));
    
    % Initialize settings and arrays for behaviour:
    trialSets = TS.sets;    
    gainDiv = seshv5lobv(iFOV).sets.start1/360*4.5*pi/trialSets.vrgain;
    smspan = trialSets.fs*msSmthing/1000; 
    vrfsDec = trialSets.fs/nDec;
    imfs = trialSets.imfs;
    vMov = zeros(size(xMov,1)/nDec,size(xMov,2),2);
    vRew = zeros(size(xRew,1)/nDec,size(xRew,2),2);
    vUnder = zeros(size(xUnder,1)/nDec,size(xUnder,2),2);
    vOver = zeros(size(xOver,1)/nDec,size(xOver,2),2);
    
    % Iterate through movements - should be the same for mov and rew
    for jMov = 1:size(xMov,2)
        x0 = xMov(:,jMov);
        x1 = smooth(x0,smspan/numel(x0),'lowess');
        v1 = diff(x1)*trialSets.fs*-gainDiv; v1(end+1) = v1(end);
        x2 = -(x1-mean(x1((trialSets.imtime-1)*trialSets.fs+1:(trialSets.imtime-1)*trialSets.fs+trialSets.fs/2)));
        vMov(:,jMov,:) = [decimate(x1,nDec) decimate(v1,nDec)];
    end
    for jRew = 1:size(xRew,2)
        x0 = xRew(:,jRew);
        x1 = smooth(x0,smspan/numel(x0),'lowess');
        v1 = diff(x1)*trialSets.fs*-gainDiv; v1(end+1) = v1(end);
        x2 = -(x1-mean(x1((trialSets.imtime-1)*trialSets.fs+1:(trialSets.imtime-1)*trialSets.fs+trialSets.fs/2)));
        vRew(:,jRew,:) = [decimate(x1,nDec) decimate(v1,nDec)];
    end
    for jUnder = 1:size(xUnder,2)
        x0 = xUnder(:,jUnder);
        x1 = smooth(x0,smspan/numel(x0),'lowess');
        v1 = diff(x1)*trialSets.fs*-gainDiv; v1(end+1) = v1(end);
        x2 = -(x1-mean(x1((trialSets.imtime-1)*trialSets.fs+1:(trialSets.imtime-1)*trialSets.fs+trialSets.fs/2)));
        vUnder(:,jUnder,:) = [decimate(x1,nDec) decimate(v1,nDec)];
    end
    for jOver = 1:size(xOver,2)
        x0 = xOver(:,jOver);
        x1 = smooth(x0,smspan/numel(x0),'lowess');
        v1 = diff(x1)*trialSets.fs*-gainDiv; v1(end+1) = v1(end);
        x2 = -(x1-mean(x1((trialSets.imtime-1)*trialSets.fs+1:(trialSets.imtime-1)*trialSets.fs+trialSets.fs/2)));
        vOver(:,jOver,:) = [decimate(x1,nDec) decimate(v1,nDec)];
    end
    
    glmv5lobv(iFOV).spkMov = spkMov(:,roisv5,:);
    glmv5lobv(iFOV).posMov = vMov(:,:,1);
    glmv5lobv(iFOV).velMov = vMov(:,:,2);
    glmv5lobv(iFOV).spkCor = spkRew(:,roisv5,:);
    glmv5lobv(iFOV).posCor = vRew(:,:,1);
    glmv5lobv(iFOV).velCor = vRew(:,:,2);
    glmv5lobv(iFOV).spkUnder = spkUnder(:,roisv5,:);
    glmv5lobv(iFOV).posUnder = vUnder(:,:,1);
    glmv5lobv(iFOV).velUnder = vUnder(:,:,2);
    glmv5lobv(iFOV).spkOver = spkOver(:,roisv5,:);
    glmv5lobv(iFOV).posOver = vOver(:,:,1);
    glmv5lobv(iFOV).velOver = vOver(:,:,2);
end

% Now do sim2
glmv5sim2 = struct;
for iFOV = 1:numel(glmv5sim20)
    TS = load(fullfile(dataPath0,seshv5sim2(iFOV).name,[seshv5sim2(iFOV).date,'_',seshv5sim2(iFOV).fov],'trialsonly.mat'));
    roisv4 = matchV5sim2(iFOV).iROIs1;
    roisv5 = matchV5sim2(iFOV).iROIs2;
    
    % Take only good ROIs
    glmv5sim2(iFOV).iROIv4 = roisv4;
    glmv5sim2(iFOV).iROIv5 = roisv5;
    glmv5sim2(iFOV).sigBool = glmv5sim20(iFOV).sigBool(roisv4);
    glmv5sim2(iFOV).sigBoolVis = glmv5sim20(iFOV).sigBoolVis(roisv4);
    glmv5sim2(iFOV).sigBoolMov = glmv5sim20(iFOV).sigBoolMov(roisv4);
    glmv5sim2(iFOV).sigBoolRew = glmv5sim20(iFOV).sigBoolRew(roisv4);
    glmv5sim2(iFOV).sigBoolLick = glmv5sim20(iFOV).sigBoolLick(roisv4);
    glmv5sim2(iFOV).gainUp = seshv5sim2(iFOV).sets.changetrial;
    glmv5sim2(iFOV).gainDown = seshv5sim2(iFOV).sets.changebacktrial;
    
    glmv5sim2(iFOV).iTrialsMov = cat(1,TS.trialstructs.LCstruct.trialind,TS.trialstructs.LUstruct.trialind,TS.trialstructs.LOstruct.trialind); 
    spkMov = cat(3,TS.trialstructs.LCstruct.spkLmov,TS.trialstructs.LUstruct.spkLmov,TS.trialstructs.LOstruct.spkLmov);
    xMov = cat(3,TS.trialstructs.LCstruct.VRLmov,TS.trialstructs.LUstruct.VRLmov,TS.trialstructs.LOstruct.VRLmov); xMov = squeeze(xMov(:,3,:));
    [~,iMovSort] = sort(glmv5sim2(iFOV).iTrialsMov,'ascend');
    glmv5sim2(iFOV).iTrialsMov = glmv5sim2(iFOV).iTrialsMov(iMovSort);
    spkMov = spkMov(:,:,iMovSort); xMov = xMov(:,iMovSort);
    
    glmv5sim2(iFOV).iTrialsRew = cat(1,TS.trialstructs.LCstruct.trialind);
    glmv5sim2(iFOV).iTrialsUnder = cat(1,TS.trialstructs.LUstruct.trialind);
    glmv5sim2(iFOV).iTrialsOver = cat(1,TS.trialstructs.LOstruct.trialind);
    spkRew = cat(3,TS.trialstructs.LCstruct.spkLoff);    
    spkUnder = cat(3,TS.trialstructs.LUstruct.spkLoff);
    spkOver = cat(3,TS.trialstructs.LOstruct.spkLoff);
    xRew = cat(3,TS.trialstructs.LCstruct.VRLoff); xRew = squeeze(xRew(:,3,:));
    xUnder = cat(3,TS.trialstructs.LUstruct.VRLoff); xUnder = squeeze(xUnder(:,3,:));
    xOver = cat(3,TS.trialstructs.LOstruct.VRLoff); xOver = squeeze(xOver(:,3,:));

    % Initialize settings and arrays for behaviour:
    trialSets = TS.sets;    
    gainDiv = seshv5sim2(iFOV).sets.start1/360*4.5*pi/trialSets.vrgain;
    smspan = trialSets.fs*msSmthing/1000; 
    vrfsDec = trialSets.fs/nDec;
    vMov = zeros(size(xMov,1)/nDec,size(xMov,2),2);
    vRew = zeros(size(xRew,1)/nDec,size(xRew,2),2);
    vUnder = zeros(size(xRew,1)/nDec,size(xUnder,2),2);
    vOver = zeros(size(xRew,1)/nDec,size(xOver,2),2);
    
    % Iterate through movements - should be the same for mov and rew
    for jMov = 1:size(xMov,2)
        x0 = xMov(:,jMov);
        x1 = smooth(x0,smspan/numel(x0),'lowess');
        v1 = diff(x1)*trialSets.fs*-gainDiv; v1(end+1) = v1(end);
        x2 = -(x1-mean(x1((trialSets.imtime-1)*trialSets.fs+1:(trialSets.imtime-1)*trialSets.fs+trialSets.fs/2)));
        vMov(:,jMov,:) = [decimate(x1,nDec) decimate(v1,nDec)];
    end
    for jRew = 1:size(xRew,2)
        x0 = xRew(:,jRew);
        x1 = smooth(x0,smspan/numel(x0),'lowess');
        v1 = diff(x1)*trialSets.fs*-gainDiv; v1(end+1) = v1(end);
        x2 = -(x1-mean(x1((trialSets.imtime-1)*trialSets.fs+1:(trialSets.imtime-1)*trialSets.fs+trialSets.fs/2)));
        vRew(:,jRew,:) = [decimate(x1,nDec) decimate(v1,nDec)];
    end
    for jUnder = 1:size(xUnder,2)
        x0 = xUnder(:,jUnder);
        x1 = smooth(x0,smspan/numel(x0),'lowess');
        v1 = diff(x1)*trialSets.fs*-gainDiv; v1(end+1) = v1(end);
        x2 = -(x1-mean(x1((trialSets.imtime-1)*trialSets.fs+1:(trialSets.imtime-1)*trialSets.fs+trialSets.fs/2)));
        vUnder(:,jUnder,:) = [decimate(x1,nDec) decimate(v1,nDec)];
    end
    for jOver = 1:size(xOver,2)
        x0 = xOver(:,jOver);
        x1 = smooth(x0,smspan/numel(x0),'lowess');
        v1 = diff(x1)*trialSets.fs*-gainDiv; v1(end+1) = v1(end);
        x2 = -(x1-mean(x1((trialSets.imtime-1)*trialSets.fs+1:(trialSets.imtime-1)*trialSets.fs+trialSets.fs/2)));
        vOver(:,jOver,:) = [decimate(x1,nDec) decimate(v1,nDec)];
    end
    
    glmv5sim2(iFOV).spkMov = spkMov(:,roisv5,:);
    glmv5sim2(iFOV).posMov = vMov(:,:,1);
    glmv5sim2(iFOV).velMov = vMov(:,:,2);
    glmv5sim2(iFOV).spkCor = spkRew(:,roisv5,:);
    glmv5sim2(iFOV).posCor = vRew(:,:,1);
    glmv5sim2(iFOV).velCor = vRew(:,:,2);
    glmv5sim2(iFOV).spkUnder = spkUnder(:,roisv5,:);
    glmv5sim2(iFOV).posUnder = vUnder(:,:,1);
    glmv5sim2(iFOV).velUnder = vUnder(:,:,2);
    glmv5sim2(iFOV).spkOver = spkOver(:,roisv5,:);
    glmv5sim2(iFOV).posOver = vOver(:,:,1);
    glmv5sim2(iFOV).velOver = vOver(:,:,2);
end

clear dataPath0 gainDiv glmv5lobv0 glmv5sim20 iFOV iMovSort jMov jRew roisv4 roisv5 smspan msSmthing TS trialSets 
clear spkMov spkRew spkUnder spkOver v1 vMov vRew vUnder vOver x0 x1 x2 xMov xRew xUnder xOver

%% Get the data:
% % Setting the ranges manually:
% glmv5lobv(1).ranges = {[-29:0],[1:15],[46:75],[101:130],[1:15]}; %dk103
% glmv5lobv(2).ranges = {[-29:0],[1:15],[31:60],[131:160],[1:15]}; %dk105
% glmv5lobv(3).ranges = {[-29:0],[1:15],[76:105],[121:150],[1:15]}; %dk169
% glmv5lobv(4).ranges = {[-29:0],[1:15],[61:90],[106:135],[1:15]}; %dk171
% glmv5lobv(5).ranges = {[-29:0],[1:15],[61:90],[106:135],[1:15]}; %dk199
% 
% glmv5sim2(1).ranges = {[-29:0],[1:15],[46:75],[86:105],[1:15]}; %dk103
% glmv5sim2(2).ranges = {[-29:0],[1:15],[31:60],[86:115],[1:15]}; %dk105
% glmv5sim2(3).ranges = {[-29:0],[1:15],[76:105],[121:150],[1:15]}; %dk169
% glmv5sim2(4).ranges = {[-29:0],[1:15],[61:90],[106:135],[1:15]}; %dk171
% glmv5sim2(5).ranges = {[-29:0],[1:15],[61:90],[106:135],[1:15]}; %dk199

% Setting the ranges programatically:
behavAdapt = load(fullfile(paperFold,'Figures\Fig1\Sub1_v5Adaptation','vAdaptationCurves.mat'));
for iBehavAll = 1:numel(behavAdapt.sessionsv5)
    for jSeshLobv = 1:numel(seshv5lobv)
        if strcmpi(behavAdapt.sessionsv5(iBehavAll).name,seshv5lobv(jSeshLobv).name) && ...
                strcmpi(behavAdapt.sessionsv5(iBehavAll).date,seshv5lobv(jSeshLobv).date)
%             iBehavLobv(jSeshLobv) = iBehavAll;
            glmv5lobv(jSeshLobv).ranges{1,1} = behavAdapt.sessionsv5(iBehavAll).endPos.iTrialsBase - ...
                glmv5lobv(jSeshLobv).gainUp;
            glmv5lobv(jSeshLobv).ranges{2,1} = behavAdapt.sessionsv5(iBehavAll).endPos.iTrialsEarly(1:end/2) - ...
                glmv5lobv(jSeshLobv).gainUp;
            glmv5lobv(jSeshLobv).ranges{3,1} = behavAdapt.sessionsv5(iBehavAll).endPos.iTrialsMid - ...
                glmv5lobv(jSeshLobv).gainUp;
            glmv5lobv(jSeshLobv).ranges{4,1} = behavAdapt.sessionsv5(iBehavAll).endPos.iTrialsLate - ...
                glmv5lobv(jSeshLobv).gainUp;
            glmv5lobv(jSeshLobv).ranges{5,1} = behavAdapt.sessionsv5(iBehavAll).endPos.iTrialsWash - ...
                glmv5lobv(jSeshLobv).gainDown;
        end
    end
    for jSeshSim2 = 1:numel(seshv5sim2)
        if strcmpi(behavAdapt.sessionsv5(iBehavAll).name,seshv5sim2(jSeshSim2).name) && ...
                strcmpi(behavAdapt.sessionsv5(iBehavAll).date,seshv5sim2(jSeshSim2).date)
            glmv5sim2(jSeshSim2).ranges{1,1} = behavAdapt.sessionsv5(iBehavAll).endPos.iTrialsBase - ...
                glmv5sim2(jSeshSim2).gainUp;
            glmv5sim2(jSeshSim2).ranges{2,1} = behavAdapt.sessionsv5(iBehavAll).endPos.iTrialsEarly(1:end/2) - ...
                glmv5sim2(jSeshSim2).gainUp;
            glmv5sim2(jSeshSim2).ranges{3,1} = behavAdapt.sessionsv5(iBehavAll).endPos.iTrialsMid - ...
                glmv5sim2(jSeshSim2).gainUp;
            glmv5sim2(jSeshSim2).ranges{4,1} = behavAdapt.sessionsv5(iBehavAll).endPos.iTrialsLate - ...
                glmv5sim2(jSeshSim2).gainUp;
            glmv5sim2(jSeshSim2).ranges{5,1} = behavAdapt.sessionsv5(iBehavAll).endPos.iTrialsWash - ...
                glmv5sim2(jSeshSim2).gainDown;
        end
    end
end

clear behavAdapt iBehavAll jSeshLobv jSeshSim2 

%% Plot 15 trial rolling average for mov trials and reward trials:
glmFold = fullfile(paperFold, 'Figures\Fig5\Sub1_FOVexamples');
if ~exist(glmFold,'dir'); mkdir(glmFold); end
plotSets.trialGrp = 30;
plotSets.imfs = imfs;
plotSets.vrfs = vrfsDec;
for iPlot = 1:numel(glmv5lobv)    
    plotSets.name = seshv5lobv(iPlot).name;
    plotSets.date = seshv5lobv(iPlot).date;
    plotSets.fov = 'lobv';
    [figOut,dataOut] = fig5_adaptPlotFOV(glmv5lobv(iPlot),plotSets);
    if savebool
        figname = [sprintf('%s_lobv_movAdapt',seshv5lobv(iPlot).name)];
        savefig(figOut(1), fullfile(glmFold,[figname,'.fig']));
        saveas(figOut(1),fullfile(glmFold,[figname,'.png']));
        print(figOut(1),fullfile(glmFold,[figname,'.eps']), '-depsc', '-vector');
        figname = [sprintf('%s_lobv_rewAdapt',seshv5lobv(iPlot).name)];
        savefig(figOut(2), fullfile(glmFold,[figname,'.fig']));
        saveas(figOut(2),fullfile(glmFold,[figname,'.png']));
        print(figOut(2),fullfile(glmFold,[figname,'.eps']), '-depsc', '-vector');
        figname = [sprintf('%s_lobv_underAdapt_rewPCs',seshv5lobv(iPlot).name)];
        savefig(figOut(3), fullfile(glmFold,[figname,'.fig']));
        saveas(figOut(3),fullfile(glmFold,[figname,'.png']));
        print(figOut(3),fullfile(glmFold,[figname,'.eps']), '-depsc', '-vector');
        figname = [sprintf('%s_lobv_overAdapt_rewPCs',seshv5lobv(iPlot).name)];
        savefig(figOut(4), fullfile(glmFold,[figname,'.fig']));
        saveas(figOut(4),fullfile(glmFold,[figname,'.png']));
        print(figOut(4),fullfile(glmFold,[figname,'.eps']), '-depsc', '-vector');
        figname = [sprintf('%s_lobv_underAdapt_movPCs',seshv5lobv(iPlot).name)];
        savefig(figOut(5), fullfile(glmFold,[figname,'.fig']));
        saveas(figOut(5),fullfile(glmFold,[figname,'.png']));
        print(figOut(5),fullfile(glmFold,[figname,'.eps']), '-depsc', '-vector');
        figname = [sprintf('%s_lobv_overAdapt_movPCs',seshv5lobv(iPlot).name)];
        savefig(figOut(6), fullfile(glmFold,[figname,'.fig']));
        saveas(figOut(6),fullfile(glmFold,[figname,'.png']));
        print(figOut(6),fullfile(glmFold,[figname,'.eps']), '-depsc', '-vector');
        
    end    
end

for iPlot = 1:numel(glmv5sim2)
    plotSets.name = seshv5sim2(iPlot).name;
    plotSets.date = seshv5sim2(iPlot).date;
    plotSets.fov = 'sim2';
    [figOut,dataOut] = fig5_adaptPlotFOV(glmv5sim2(iPlot),plotSets);
    if savebool
        figname = [sprintf('%s_sim2_movAdapt',seshv5sim2(iPlot).name)];
        savefig(figOut(1), fullfile(glmFold,[figname,'.fig']));
        saveas(figOut(1),fullfile(glmFold,[figname,'.png']));
        print(figOut(1),fullfile(glmFold,[figname,'.eps']), '-depsc', '-vector');
        figname = [sprintf('%s_sim2_rewAdapt',seshv5sim2(iPlot).name)];
        savefig(figOut(2), fullfile(glmFold,[figname,'.fig']));
        saveas(figOut(2),fullfile(glmFold,[figname,'.png']));
        print(figOut(2),fullfile(glmFold,[figname,'.eps']), '-depsc', '-vector');
        figname = [sprintf('%s_sim2_underAdapt_rewPCs',seshv5lobv(iPlot).name)];
        savefig(figOut(3), fullfile(glmFold,[figname,'.fig']));
        saveas(figOut(3),fullfile(glmFold,[figname,'.png']));
        print(figOut(3),fullfile(glmFold,[figname,'.eps']), '-depsc', '-vector');
        figname = [sprintf('%s_sim2_overAdapt_rewPCs',seshv5lobv(iPlot).name)];
        savefig(figOut(4), fullfile(glmFold,[figname,'.fig']));
        saveas(figOut(4),fullfile(glmFold,[figname,'.png']));
        print(figOut(4),fullfile(glmFold,[figname,'.eps']), '-depsc', '-vector');
        figname = [sprintf('%s_sim2_underAdapt_movPCs',seshv5lobv(iPlot).name)];
        savefig(figOut(5), fullfile(glmFold,[figname,'.fig']));
        saveas(figOut(5),fullfile(glmFold,[figname,'.png']));
        print(figOut(5),fullfile(glmFold,[figname,'.eps']), '-depsc', '-vector');
        figname = [sprintf('%s_sim2_overAdapt_movPCs',seshv5lobv(iPlot).name)];
        savefig(figOut(6), fullfile(glmFold,[figname,'.fig']));
        saveas(figOut(6),fullfile(glmFold,[figname,'.png']));
        print(figOut(6),fullfile(glmFold,[figname,'.eps']), '-depsc', '-vector');
    end    
end

clear dataOut figOut iPlot plotSets figname

%%
% rangeBase = [-29:0]; rangeEarly = [1:15]; rangeMid = [46:75]; rangeLate = [96:125]; rangeWash = [1:15];


lobvdata = struct;
for iFOV = 1:numel(glmv5lobv)
    lobvdata(iFOV).name = seshv5lobv(iFOV).name;
    lobvdata(iFOV).iROIv4 = glmv5lobv(iFOV).iROIv4;
    lobvdata(iFOV).iROIv5 = glmv5lobv(iFOV).iROIv5;
    lobvdata(iFOV).sigBoolMov = glmv5lobv(iFOV).sigBoolMov;
    lobvdata(iFOV).sigBoolRew = glmv5lobv(iFOV).sigBoolRew;
    % Get trials for lobv, movement
    iTrBase = ismember(glmv5lobv(iFOV).iTrialsMov,glmv5lobv(iFOV).gainUp + glmv5lobv(iFOV).ranges{1});
    iTrEarly = ismember(glmv5lobv(iFOV).iTrialsMov,glmv5lobv(iFOV).gainUp + glmv5lobv(iFOV).ranges{2});
    iTrMid = ismember(glmv5lobv(iFOV).iTrialsMov,glmv5lobv(iFOV).gainUp + glmv5lobv(iFOV).ranges{3});
    iTrLate = ismember(glmv5lobv(iFOV).iTrialsMov,glmv5lobv(iFOV).gainUp + glmv5lobv(iFOV).ranges{4});
    iTrWash = ismember(glmv5lobv(iFOV).iTrialsMov,glmv5lobv(iFOV).gainDown + glmv5lobv(iFOV).ranges{5});    
    nTrialsMov = [sum(iTrBase) sum(iTrEarly) sum(iTrMid) sum(iTrLate) sum(iTrWash)];
    fprintf('Movement trials across adaptation: %s \n',num2str(nTrialsMov))
    
    % Get averages
    lobvdata(iFOV).spkMovBase = mean(glmv5lobv(iFOV).spkMov(end/2-imfs+1:end/2+imfs,glmv5lobv(iFOV).sigBoolMov,iTrBase),3);
    lobvdata(iFOV).spkMovBaseFR = mean(lobvdata(iFOV).spkMovBase(imfs-.3*imfs+1:imfs,:))';
    lobvdata(iFOV).spkMovBaseMu = mean(lobvdata(iFOV).spkMovBase,2);
    lobvdata(iFOV).velMovBase = mean(glmv5lobv(iFOV).velMov(end/2-vrfsDec+1:end/2+vrfsDec,iTrBase),2);
    
    lobvdata(iFOV).spkMovEarly = mean(glmv5lobv(iFOV).spkMov(end/2-imfs+1:end/2+imfs,glmv5lobv(iFOV).sigBoolMov,iTrEarly),3);
    lobvdata(iFOV).spkMovEarlyFR = mean(lobvdata(iFOV).spkMovEarly(imfs-.3*imfs+1:imfs,:))';
    lobvdata(iFOV).spkMovEarlyMu = mean(lobvdata(iFOV).spkMovEarly,2);
    lobvdata(iFOV).velMovEarly = mean(glmv5lobv(iFOV).velMov(end/2-vrfsDec+1:end/2+vrfsDec,iTrEarly),2);
    
    lobvdata(iFOV).spkMovMid = mean(glmv5lobv(iFOV).spkMov(end/2-imfs+1:end/2+imfs,glmv5lobv(iFOV).sigBoolMov,iTrMid),3);
    lobvdata(iFOV).spkMovMidFR = mean(lobvdata(iFOV).spkMovMid(imfs-.3*imfs+1:imfs,:))';
    lobvdata(iFOV).spkMovMidMu = mean(lobvdata(iFOV).spkMovMid,2);
    lobvdata(iFOV).velMovMid = mean(glmv5lobv(iFOV).velMov(end/2-vrfsDec+1:end/2+vrfsDec,iTrMid),2);
    
    lobvdata(iFOV).spkMovLate = mean(glmv5lobv(iFOV).spkMov(end/2-imfs+1:end/2+imfs,glmv5lobv(iFOV).sigBoolMov,iTrLate),3);
    lobvdata(iFOV).spkMovLateFR = mean(lobvdata(iFOV).spkMovLate(imfs-.3*imfs+1:imfs,:))';
    lobvdata(iFOV).spkMovLateMu = mean(lobvdata(iFOV).spkMovLate,2);
    lobvdata(iFOV).velMovLate = mean(glmv5lobv(iFOV).velMov(end/2-vrfsDec+1:end/2+vrfsDec,iTrLate),2);
    
    lobvdata(iFOV).spkMovWash = mean(glmv5lobv(iFOV).spkMov(end/2-imfs+1:end/2+imfs,glmv5lobv(iFOV).sigBoolMov,iTrWash),3);
    lobvdata(iFOV).spkMovWashFR = mean(lobvdata(iFOV).spkMovWash(imfs-.3*imfs+1:imfs,:))';
    lobvdata(iFOV).spkMovWashMu = mean(lobvdata(iFOV).spkMovWash,2);
    lobvdata(iFOV).velMovWash = mean(glmv5lobv(iFOV).velMov(end/2-vrfsDec+1:end/2+vrfsDec,iTrWash),2);
        
    
    % Get trials for lobv, reward
    iTrBase = ismember(glmv5lobv(iFOV).iTrialsRew,glmv5lobv(iFOV).gainUp + glmv5lobv(iFOV).ranges{1});
    iTrEarly = ismember(glmv5lobv(iFOV).iTrialsRew,glmv5lobv(iFOV).gainUp + glmv5lobv(iFOV).ranges{2});
    iTrMid = ismember(glmv5lobv(iFOV).iTrialsRew,glmv5lobv(iFOV).gainUp + glmv5lobv(iFOV).ranges{3});
    iTrLate = ismember(glmv5lobv(iFOV).iTrialsRew,glmv5lobv(iFOV).gainUp + glmv5lobv(iFOV).ranges{4});
    iTrWash = ismember(glmv5lobv(iFOV).iTrialsRew,glmv5lobv(iFOV).gainDown + glmv5lobv(iFOV).ranges{5});    
    nTrialsRew = [sum(iTrBase) sum(iTrEarly) sum(iTrMid) sum(iTrLate) sum(iTrWash)];
    fprintf('Reward trials across adaptation: %s \n',num2str(nTrialsRew))    
    
    % Get averages
    lobvdata(iFOV).spkRewBase = mean(glmv5lobv(iFOV).spkCor(end/2-imfs+1:end/2+imfs,glmv5lobv(iFOV).sigBoolRew,iTrBase),3);
    lobvdata(iFOV).spkRewBaseFR = mean(lobvdata(iFOV).spkRewBase(imfs+1:imfs+.2*imfs,:))';
    lobvdata(iFOV).spkRewBaseMu = mean(lobvdata(iFOV).spkRewBase,2);
    lobvdata(iFOV).velRewBase = mean(glmv5lobv(iFOV).velCor(end/2-vrfsDec+1:end/2+vrfsDec,iTrBase),2);
    
    lobvdata(iFOV).spkRewEarly = mean(glmv5lobv(iFOV).spkCor(end/2-imfs+1:end/2+imfs,glmv5lobv(iFOV).sigBoolRew,iTrEarly),3);
    lobvdata(iFOV).spkRewEarlyFR = mean(lobvdata(iFOV).spkRewEarly(imfs+1:imfs+.2*imfs,:))';
    lobvdata(iFOV).spkRewEarlyMu = mean(lobvdata(iFOV).spkRewEarly,2);
    lobvdata(iFOV).velRewEarly = mean(glmv5lobv(iFOV).velCor(end/2-vrfsDec+1:end/2+vrfsDec,iTrEarly),2);
    
    lobvdata(iFOV).spkRewMid = mean(glmv5lobv(iFOV).spkCor(end/2-imfs+1:end/2+imfs,glmv5lobv(iFOV).sigBoolRew,iTrMid),3);
    lobvdata(iFOV).spkRewMidFR = mean(lobvdata(iFOV).spkRewMid(imfs+1:imfs+.2*imfs,:))';
    lobvdata(iFOV).spkRewMidMu = mean(lobvdata(iFOV).spkRewMid,2);
    lobvdata(iFOV).velRewMid = mean(glmv5lobv(iFOV).velCor(end/2-vrfsDec+1:end/2+vrfsDec,iTrMid),2);
    
    lobvdata(iFOV).spkRewLate = mean(glmv5lobv(iFOV).spkCor(end/2-imfs+1:end/2+imfs,glmv5lobv(iFOV).sigBoolRew,iTrLate),3);
    lobvdata(iFOV).spkRewLateFR = mean(lobvdata(iFOV).spkRewLate(imfs+1:imfs+.2*imfs,:))';
    lobvdata(iFOV).spkRewLateMu = mean(lobvdata(iFOV).spkRewLate,2);
    lobvdata(iFOV).velRewLate = mean(glmv5lobv(iFOV).velCor(end/2-vrfsDec+1:end/2+vrfsDec,iTrLate),2);
    
    lobvdata(iFOV).spkRewWash = mean(glmv5lobv(iFOV).spkCor(end/2-imfs+1:end/2+imfs,glmv5lobv(iFOV).sigBoolRew,iTrWash),3); 
    lobvdata(iFOV).spkRewWashFR = mean(lobvdata(iFOV).spkRewWash(imfs+1:imfs+.2*imfs,:))';
    lobvdata(iFOV).spkRewWashMu = mean(lobvdata(iFOV).spkRewWash,2);
    lobvdata(iFOV).velRewWash = mean(glmv5lobv(iFOV).velCor(end/2-vrfsDec+1:end/2+vrfsDec,iTrWash),2);

    % Get trials for lobv, undershoots
    iTrBase = ismember(glmv5lobv(iFOV).iTrialsUnder,glmv5lobv(iFOV).gainUp + glmv5lobv(iFOV).ranges{1});
    iTrEarly = ismember(glmv5lobv(iFOV).iTrialsUnder,glmv5lobv(iFOV).gainUp + glmv5lobv(iFOV).ranges{2});
    iTrMid = ismember(glmv5lobv(iFOV).iTrialsUnder,glmv5lobv(iFOV).gainUp + glmv5lobv(iFOV).ranges{3});
    iTrLate = ismember(glmv5lobv(iFOV).iTrialsUnder,glmv5lobv(iFOV).gainUp + glmv5lobv(iFOV).ranges{4});
    iTrWash = ismember(glmv5lobv(iFOV).iTrialsUnder,glmv5lobv(iFOV).gainDown + glmv5lobv(iFOV).ranges{5});    
    nTrialsRew = [sum(iTrBase) sum(iTrEarly) sum(iTrMid) sum(iTrLate) sum(iTrWash)];
    fprintf('Undershoot trials across adaptation: %s \n',num2str(nTrialsRew))    
    
    % Get averages - reward PCs
    lobvdata(iFOV).spkUnderBase = mean(glmv5lobv(iFOV).spkUnder(end/2-imfs+1:end/2+imfs,glmv5lobv(iFOV).sigBoolRew,iTrBase),3);
    lobvdata(iFOV).spkUnderBaseFR = mean(lobvdata(iFOV).spkUnderBase(imfs+1:imfs+.2*imfs,:))';
    lobvdata(iFOV).spkUnderBaseMu = mean(lobvdata(iFOV).spkUnderBase,2);
    lobvdata(iFOV).velUnderBase = mean(glmv5lobv(iFOV).velUnder(end/2-vrfsDec+1:end/2+vrfsDec,iTrBase),2);
    
    lobvdata(iFOV).spkUnderEarly = mean(glmv5lobv(iFOV).spkUnder(end/2-imfs+1:end/2+imfs,glmv5lobv(iFOV).sigBoolRew,iTrEarly),3);
    lobvdata(iFOV).spkUnderEarlyFR = mean(lobvdata(iFOV).spkUnderEarly(imfs+1:imfs+.2*imfs,:))';
    lobvdata(iFOV).spkUnderEarlyMu = mean(lobvdata(iFOV).spkUnderEarly,2);
    lobvdata(iFOV).velUnderEarly = mean(glmv5lobv(iFOV).velUnder(end/2-vrfsDec+1:end/2+vrfsDec,iTrEarly),2);
    
    lobvdata(iFOV).spkUnderMid = mean(glmv5lobv(iFOV).spkUnder(end/2-imfs+1:end/2+imfs,glmv5lobv(iFOV).sigBoolRew,iTrMid),3);
    lobvdata(iFOV).spkUnderMidFR = mean(lobvdata(iFOV).spkUnderMid(imfs+1:imfs+.2*imfs,:))';
    lobvdata(iFOV).spkUnderMidMu = mean(lobvdata(iFOV).spkUnderMid,2);
    lobvdata(iFOV).velUnderMid = mean(glmv5lobv(iFOV).velUnder(end/2-vrfsDec+1:end/2+vrfsDec,iTrMid),2);
    
    lobvdata(iFOV).spkUnderLate = mean(glmv5lobv(iFOV).spkUnder(end/2-imfs+1:end/2+imfs,glmv5lobv(iFOV).sigBoolRew,iTrLate),3);
    lobvdata(iFOV).spkUnderLateFR = mean(lobvdata(iFOV).spkUnderLate(imfs+1:imfs+.2*imfs,:))';
    lobvdata(iFOV).spkUnderLateMu = mean(lobvdata(iFOV).spkUnderLate,2);
    lobvdata(iFOV).velUnderLate = mean(glmv5lobv(iFOV).velUnder(end/2-vrfsDec+1:end/2+vrfsDec,iTrLate),2);
    
    lobvdata(iFOV).spkUnderWash = mean(glmv5lobv(iFOV).spkUnder(end/2-imfs+1:end/2+imfs,glmv5lobv(iFOV).sigBoolRew,iTrWash),3); 
    lobvdata(iFOV).spkUnderWashFR = mean(lobvdata(iFOV).spkUnderWash(imfs+1:imfs+.2*imfs,:))';
    lobvdata(iFOV).spkUnderWashMu = mean(lobvdata(iFOV).spkUnderWash,2);
    lobvdata(iFOV).velUnderWash = mean(glmv5lobv(iFOV).velUnder(end/2-vrfsDec+1:end/2+vrfsDec,iTrWash),2);

    % Get averages - movement PCs
    lobvdata(iFOV).spkUnder_movBase = mean(glmv5lobv(iFOV).spkUnder(end/2-imfs+1:end/2+imfs,glmv5lobv(iFOV).sigBoolMov,iTrBase),3);
    lobvdata(iFOV).spkUnder_movBaseFR = mean(lobvdata(iFOV).spkUnder_movBase(imfs+1:imfs+.2*imfs,:))';
    lobvdata(iFOV).spkUnder_movBaseMu = mean(lobvdata(iFOV).spkUnder_movBase,2);
    lobvdata(iFOV).velUnderBase = mean(glmv5lobv(iFOV).velUnder(end/2-vrfsDec+1:end/2+vrfsDec,iTrBase),2);
    
    lobvdata(iFOV).spkUnder_movEarly = mean(glmv5lobv(iFOV).spkUnder(end/2-imfs+1:end/2+imfs,glmv5lobv(iFOV).sigBoolMov,iTrEarly),3);
    lobvdata(iFOV).spkUnder_movEarlyFR = mean(lobvdata(iFOV).spkUnder_movEarly(imfs+1:imfs+.2*imfs,:))';
    lobvdata(iFOV).spkUnder_movEarlyMu = mean(lobvdata(iFOV).spkUnder_movEarly,2);
    lobvdata(iFOV).velUnderEarly = mean(glmv5lobv(iFOV).velUnder(end/2-vrfsDec+1:end/2+vrfsDec,iTrEarly),2);
    
    lobvdata(iFOV).spkUnder_movMid = mean(glmv5lobv(iFOV).spkUnder(end/2-imfs+1:end/2+imfs,glmv5lobv(iFOV).sigBoolMov,iTrMid),3);
    lobvdata(iFOV).spkUnder_movMidFR = mean(lobvdata(iFOV).spkUnder_movMid(imfs+1:imfs+.2*imfs,:))';
    lobvdata(iFOV).spkUnder_movMidMu = mean(lobvdata(iFOV).spkUnder_movMid,2);
    lobvdata(iFOV).velUnderMid = mean(glmv5lobv(iFOV).velUnder(end/2-vrfsDec+1:end/2+vrfsDec,iTrMid),2);
    
    lobvdata(iFOV).spkUnder_movLate = mean(glmv5lobv(iFOV).spkUnder(end/2-imfs+1:end/2+imfs,glmv5lobv(iFOV).sigBoolMov,iTrLate),3);
    lobvdata(iFOV).spkUnder_movLateFR = mean(lobvdata(iFOV).spkUnder_movLate(imfs+1:imfs+.2*imfs,:))';
    lobvdata(iFOV).spkUnder_movLateMu = mean(lobvdata(iFOV).spkUnder_movLate,2);
    lobvdata(iFOV).velUnderLate = mean(glmv5lobv(iFOV).velUnder(end/2-vrfsDec+1:end/2+vrfsDec,iTrLate),2);
    
    lobvdata(iFOV).spkUnder_movWash = mean(glmv5lobv(iFOV).spkUnder(end/2-imfs+1:end/2+imfs,glmv5lobv(iFOV).sigBoolMov,iTrWash),3); 
    lobvdata(iFOV).spkUnder_movWashFR = mean(lobvdata(iFOV).spkUnder_movWash(imfs+1:imfs+.2*imfs,:))';
    lobvdata(iFOV).spkUnder_movWashMu = mean(lobvdata(iFOV).spkUnder_movWash,2);
    lobvdata(iFOV).velUnderWash = mean(glmv5lobv(iFOV).velUnder(end/2-vrfsDec+1:end/2+vrfsDec,iTrWash),2);

    % Get trials for lobv, overshoot
    iTrBase = ismember(glmv5lobv(iFOV).iTrialsOver,glmv5lobv(iFOV).gainUp + glmv5lobv(iFOV).ranges{1});
    iTrEarly = ismember(glmv5lobv(iFOV).iTrialsOver,glmv5lobv(iFOV).gainUp + glmv5lobv(iFOV).ranges{2});
    iTrMid = ismember(glmv5lobv(iFOV).iTrialsOver,glmv5lobv(iFOV).gainUp + glmv5lobv(iFOV).ranges{3});
    iTrLate = ismember(glmv5lobv(iFOV).iTrialsOver,glmv5lobv(iFOV).gainUp + glmv5lobv(iFOV).ranges{4});
    iTrWash = ismember(glmv5lobv(iFOV).iTrialsOver,glmv5lobv(iFOV).gainDown + glmv5lobv(iFOV).ranges{5});    
    nTrialsRew = [sum(iTrBase) sum(iTrEarly) sum(iTrMid) sum(iTrLate) sum(iTrWash)];
    fprintf('Overshoot trials across adaptation: %s \n',num2str(nTrialsRew))    
    
    % Get averages - rewPCs
    lobvdata(iFOV).spkOverBase = mean(glmv5lobv(iFOV).spkOver(end/2-imfs+1:end/2+imfs,glmv5lobv(iFOV).sigBoolRew,iTrBase),3);
    lobvdata(iFOV).spkOverBaseFR = mean(lobvdata(iFOV).spkOverBase(imfs+1:imfs+.2*imfs,:))';
    lobvdata(iFOV).spkOverBaseMu = mean(lobvdata(iFOV).spkOverBase,2);
    lobvdata(iFOV).velOverBase = mean(glmv5lobv(iFOV).velOver(end/2-vrfsDec+1:end/2+vrfsDec,iTrBase),2);
    
    lobvdata(iFOV).spkOverEarly = mean(glmv5lobv(iFOV).spkOver(end/2-imfs+1:end/2+imfs,glmv5lobv(iFOV).sigBoolRew,iTrEarly),3);
    lobvdata(iFOV).spkOverEarlyFR = mean(lobvdata(iFOV).spkOverEarly(imfs+1:imfs+.2*imfs,:))';
    lobvdata(iFOV).spkOverEarlyMu = mean(lobvdata(iFOV).spkOverEarly,2);
    lobvdata(iFOV).velOverEarly = mean(glmv5lobv(iFOV).velOver(end/2-vrfsDec+1:end/2+vrfsDec,iTrEarly),2);
    
    lobvdata(iFOV).spkOverMid = mean(glmv5lobv(iFOV).spkOver(end/2-imfs+1:end/2+imfs,glmv5lobv(iFOV).sigBoolRew,iTrMid),3);
    lobvdata(iFOV).spkOverMidFR = mean(lobvdata(iFOV).spkOverMid(imfs+1:imfs+.2*imfs,:))';
    lobvdata(iFOV).spkOverMidMu = mean(lobvdata(iFOV).spkOverMid,2);
    lobvdata(iFOV).velOverMid = mean(glmv5lobv(iFOV).velOver(end/2-vrfsDec+1:end/2+vrfsDec,iTrMid),2);
    
    lobvdata(iFOV).spkOverLate = mean(glmv5lobv(iFOV).spkOver(end/2-imfs+1:end/2+imfs,glmv5lobv(iFOV).sigBoolRew,iTrLate),3);
    lobvdata(iFOV).spkOverLateFR = mean(lobvdata(iFOV).spkOverLate(imfs+1:imfs+.2*imfs,:))';
    lobvdata(iFOV).spkOverLateMu = mean(lobvdata(iFOV).spkOverLate,2);
    lobvdata(iFOV).velOverLate = mean(glmv5lobv(iFOV).velOver(end/2-vrfsDec+1:end/2+vrfsDec,iTrLate),2);
    
    lobvdata(iFOV).spkOverWash = mean(glmv5lobv(iFOV).spkOver(end/2-imfs+1:end/2+imfs,glmv5lobv(iFOV).sigBoolRew,iTrWash),3); 
    lobvdata(iFOV).spkOverWashFR = mean(lobvdata(iFOV).spkOverWash(imfs+1:imfs+.2*imfs,:))';
    lobvdata(iFOV).spkOverWashMu = mean(lobvdata(iFOV).spkOverWash,2);
    lobvdata(iFOV).velOverWash = mean(glmv5lobv(iFOV).velOver(end/2-vrfsDec+1:end/2+vrfsDec,iTrWash),2);

    % Get averages - mov PCs
    lobvdata(iFOV).spkOver_movBase = mean(glmv5lobv(iFOV).spkOver(end/2-imfs+1:end/2+imfs,glmv5lobv(iFOV).sigBoolMov,iTrBase),3);
    lobvdata(iFOV).spkOver_movBaseFR = mean(lobvdata(iFOV).spkOver_movBase(imfs+1:imfs+.2*imfs,:))';
    lobvdata(iFOV).spkOver_movBaseMu = mean(lobvdata(iFOV).spkOver_movBase,2);
    lobvdata(iFOV).velOverBase = mean(glmv5lobv(iFOV).velOver(end/2-vrfsDec+1:end/2+vrfsDec,iTrBase),2);
    
    lobvdata(iFOV).spkOver_movEarly = mean(glmv5lobv(iFOV).spkOver(end/2-imfs+1:end/2+imfs,glmv5lobv(iFOV).sigBoolMov,iTrEarly),3);
    lobvdata(iFOV).spkOver_movEarlyFR = mean(lobvdata(iFOV).spkOver_movEarly(imfs+1:imfs+.2*imfs,:))';
    lobvdata(iFOV).spkOver_movEarlyMu = mean(lobvdata(iFOV).spkOver_movEarly,2);
    lobvdata(iFOV).velOverEarly = mean(glmv5lobv(iFOV).velOver(end/2-vrfsDec+1:end/2+vrfsDec,iTrEarly),2);
    
    lobvdata(iFOV).spkOver_movMid = mean(glmv5lobv(iFOV).spkOver(end/2-imfs+1:end/2+imfs,glmv5lobv(iFOV).sigBoolMov,iTrMid),3);
    lobvdata(iFOV).spkOver_movMidFR = mean(lobvdata(iFOV).spkOver_movMid(imfs+1:imfs+.2*imfs,:))';
    lobvdata(iFOV).spkOver_movMidMu = mean(lobvdata(iFOV).spkOver_movMid,2);
    lobvdata(iFOV).velOverMid = mean(glmv5lobv(iFOV).velOver(end/2-vrfsDec+1:end/2+vrfsDec,iTrMid),2);
    
    lobvdata(iFOV).spkOver_movLate = mean(glmv5lobv(iFOV).spkOver(end/2-imfs+1:end/2+imfs,glmv5lobv(iFOV).sigBoolMov,iTrLate),3);
    lobvdata(iFOV).spkOver_movLateFR = mean(lobvdata(iFOV).spkOver_movLate(imfs+1:imfs+.2*imfs,:))';
    lobvdata(iFOV).spkOver_movLateMu = mean(lobvdata(iFOV).spkOver_movLate,2);
    lobvdata(iFOV).velOverLate = mean(glmv5lobv(iFOV).velOver(end/2-vrfsDec+1:end/2+vrfsDec,iTrLate),2);
    
    lobvdata(iFOV).spkOver_movWash = mean(glmv5lobv(iFOV).spkOver(end/2-imfs+1:end/2+imfs,glmv5lobv(iFOV).sigBoolMov,iTrWash),3); 
    lobvdata(iFOV).spkOver_movWashFR = mean(lobvdata(iFOV).spkOver_movWash(imfs+1:imfs+.2*imfs,:))';
    lobvdata(iFOV).spkOver_movWashMu = mean(lobvdata(iFOV).spkOver_movWash,2);
    lobvdata(iFOV).velOverWash = mean(glmv5lobv(iFOV).velOver(end/2-vrfsDec+1:end/2+vrfsDec,iTrWash),2);
end

sim2data = struct;
for iFOV = 1:numel(glmv5sim2)
    sim2data(iFOV).name = seshv5sim2(iFOV).name;
    sim2data(iFOV).iROIv4 = glmv5sim2(iFOV).iROIv4;
    sim2data(iFOV).iROIv5 = glmv5sim2(iFOV).iROIv5;
    sim2data(iFOV).sigBoolMov = glmv5sim2(iFOV).sigBoolMov;
    sim2data(iFOV).sigBoolRew = glmv5sim2(iFOV).sigBoolRew;
    % Get trials for sim2, movement
    iTrBase = ismember(glmv5sim2(iFOV).iTrialsMov,glmv5sim2(iFOV).gainUp + glmv5sim2(iFOV).ranges{1});
    iTrEarly = ismember(glmv5sim2(iFOV).iTrialsMov,glmv5sim2(iFOV).gainUp + glmv5sim2(iFOV).ranges{2});
    iTrMid = ismember(glmv5sim2(iFOV).iTrialsMov,glmv5sim2(iFOV).gainUp + glmv5sim2(iFOV).ranges{3});
    iTrLate = ismember(glmv5sim2(iFOV).iTrialsMov,glmv5sim2(iFOV).gainUp + glmv5sim2(iFOV).ranges{4});
    iTrWash = ismember(glmv5sim2(iFOV).iTrialsMov,glmv5sim2(iFOV).gainDown + glmv5sim2(iFOV).ranges{5});    
    nTrialsMov = [sum(iTrBase) sum(iTrEarly) sum(iTrMid) sum(iTrLate) sum(iTrWash)];
    fprintf('Movement trials across adaptation: %s \n',num2str(nTrialsMov))
    
    % Get averages
    sim2data(iFOV).spkMovBase = mean(glmv5sim2(iFOV).spkMov(end/2-imfs+1:end/2+imfs,glmv5sim2(iFOV).sigBoolMov,iTrBase),3);
    sim2data(iFOV).spkMovBaseFR = mean(sim2data(iFOV).spkMovBase(imfs-.3*imfs+1:imfs,:))';
    sim2data(iFOV).spkMovBaseMu = mean(sim2data(iFOV).spkMovBase,2);
    sim2data(iFOV).velMovBase = mean(glmv5sim2(iFOV).velMov(end/2-vrfsDec+1:end/2+vrfsDec,iTrBase),2);
    
    sim2data(iFOV).spkMovEarly = mean(glmv5sim2(iFOV).spkMov(end/2-imfs+1:end/2+imfs,glmv5sim2(iFOV).sigBoolMov,iTrEarly),3);
    sim2data(iFOV).spkMovEarlyFR = mean(sim2data(iFOV).spkMovEarly(imfs-.3*imfs+1:imfs,:))';
    sim2data(iFOV).spkMovEarlyMu = mean(sim2data(iFOV).spkMovEarly,2);
    sim2data(iFOV).velMovEarly = mean(glmv5sim2(iFOV).velMov(end/2-vrfsDec+1:end/2+vrfsDec,iTrEarly),2);
    
    sim2data(iFOV).spkMovMid = mean(glmv5sim2(iFOV).spkMov(end/2-imfs+1:end/2+imfs,glmv5sim2(iFOV).sigBoolMov,iTrMid),3);
    sim2data(iFOV).spkMovMidFR = mean(sim2data(iFOV).spkMovMid(imfs-.3*imfs+1:imfs,:))';
    sim2data(iFOV).spkMovMidMu = mean(sim2data(iFOV).spkMovMid,2);
    sim2data(iFOV).velMovMid = mean(glmv5sim2(iFOV).velMov(end/2-vrfsDec+1:end/2+vrfsDec,iTrMid),2);
    
    sim2data(iFOV).spkMovLate = mean(glmv5sim2(iFOV).spkMov(end/2-imfs+1:end/2+imfs,glmv5sim2(iFOV).sigBoolMov,iTrLate),3);
    sim2data(iFOV).spkMovLateFR = mean(sim2data(iFOV).spkMovLate(imfs-.3*imfs+1:imfs,:))';
    sim2data(iFOV).spkMovLateMu = mean(sim2data(iFOV).spkMovLate,2);
    sim2data(iFOV).velMovLate = mean(glmv5sim2(iFOV).velMov(end/2-vrfsDec+1:end/2+vrfsDec,iTrLate),2);
    
    sim2data(iFOV).spkMovWash = mean(glmv5sim2(iFOV).spkMov(end/2-imfs+1:end/2+imfs,glmv5sim2(iFOV).sigBoolMov,iTrWash),3);
    sim2data(iFOV).spkMovWashFR = mean(sim2data(iFOV).spkMovWash(imfs-.3*imfs+1:imfs,:))';
    sim2data(iFOV).spkMovWashMu = mean(sim2data(iFOV).spkMovWash,2);
    sim2data(iFOV).velMovWash = mean(glmv5sim2(iFOV).velMov(end/2-vrfsDec+1:end/2+vrfsDec,iTrWash),2);
        
    
    % Get trials for sim2, reward
    iTrBase = ismember(glmv5sim2(iFOV).iTrialsRew,glmv5sim2(iFOV).gainUp + glmv5sim2(iFOV).ranges{1});
    iTrEarly = ismember(glmv5sim2(iFOV).iTrialsRew,glmv5sim2(iFOV).gainUp + glmv5sim2(iFOV).ranges{2});
    iTrMid = ismember(glmv5sim2(iFOV).iTrialsRew,glmv5sim2(iFOV).gainUp + glmv5sim2(iFOV).ranges{3});
    iTrLate = ismember(glmv5sim2(iFOV).iTrialsRew,glmv5sim2(iFOV).gainUp + glmv5sim2(iFOV).ranges{4});
    iTrWash = ismember(glmv5sim2(iFOV).iTrialsRew,glmv5sim2(iFOV).gainDown + glmv5sim2(iFOV).ranges{5});    
    nTrialsRew = [sum(iTrBase) sum(iTrEarly) sum(iTrMid) sum(iTrLate) sum(iTrWash)];
    fprintf('Reward trials across adaptation: %s \n',num2str(nTrialsRew))    
    
    % Get averages
    sim2data(iFOV).spkRewBase = mean(glmv5sim2(iFOV).spkCor(end/2-imfs+1:end/2+imfs,glmv5sim2(iFOV).sigBoolRew,iTrBase),3);
    sim2data(iFOV).spkRewBaseFR = mean(sim2data(iFOV).spkRewBase(imfs+1:imfs+.2*imfs,:))';
    sim2data(iFOV).spkRewBaseMu = mean(sim2data(iFOV).spkRewBase,2);
    sim2data(iFOV).velRewBase = mean(glmv5sim2(iFOV).velCor(end/2-vrfsDec+1:end/2+vrfsDec,iTrBase),2);
    
    sim2data(iFOV).spkRewEarly = mean(glmv5sim2(iFOV).spkCor(end/2-imfs+1:end/2+imfs,glmv5sim2(iFOV).sigBoolRew,iTrEarly),3);
    sim2data(iFOV).spkRewEarlyFR = mean(sim2data(iFOV).spkRewEarly(imfs+1:imfs+.2*imfs,:))';
    sim2data(iFOV).spkRewEarlyMu = mean(sim2data(iFOV).spkRewEarly,2);
    sim2data(iFOV).velRewEarly = mean(glmv5sim2(iFOV).velCor(end/2-vrfsDec+1:end/2+vrfsDec,iTrEarly),2);
    
    sim2data(iFOV).spkRewMid = mean(glmv5sim2(iFOV).spkCor(end/2-imfs+1:end/2+imfs,glmv5sim2(iFOV).sigBoolRew,iTrMid),3);
    sim2data(iFOV).spkRewMidFR = mean(sim2data(iFOV).spkRewMid(imfs+1:imfs+.2*imfs,:))';
    sim2data(iFOV).spkRewMidMu = mean(sim2data(iFOV).spkRewMid,2);
    sim2data(iFOV).velRewMid = mean(glmv5sim2(iFOV).velCor(end/2-vrfsDec+1:end/2+vrfsDec,iTrMid),2);
    
    sim2data(iFOV).spkRewLate = mean(glmv5sim2(iFOV).spkCor(end/2-imfs+1:end/2+imfs,glmv5sim2(iFOV).sigBoolRew,iTrLate),3);
    sim2data(iFOV).spkRewLateFR = mean(sim2data(iFOV).spkRewLate(imfs+1:imfs+.2*imfs,:))';
    sim2data(iFOV).spkRewLateMu = mean(sim2data(iFOV).spkRewLate,2);
    sim2data(iFOV).velRewLate = mean(glmv5sim2(iFOV).velCor(end/2-vrfsDec+1:end/2+vrfsDec,iTrLate),2);
    
    sim2data(iFOV).spkRewWash = mean(glmv5sim2(iFOV).spkCor(end/2-imfs+1:end/2+imfs,glmv5sim2(iFOV).sigBoolRew,iTrWash),3); 
    sim2data(iFOV).spkRewWashFR = mean(sim2data(iFOV).spkRewWash(imfs+1:imfs+.2*imfs,:))';
    sim2data(iFOV).spkRewWashMu = mean(sim2data(iFOV).spkRewWash,2);
    sim2data(iFOV).velRewWash = mean(glmv5sim2(iFOV).velCor(end/2-vrfsDec+1:end/2+vrfsDec,iTrWash),2);        
    

    % Get trials for lobv, undershoots
    iTrBase = ismember(glmv5sim2(iFOV).iTrialsUnder,glmv5sim2(iFOV).gainUp + glmv5sim2(iFOV).ranges{1});
    iTrEarly = ismember(glmv5sim2(iFOV).iTrialsUnder,glmv5sim2(iFOV).gainUp + glmv5sim2(iFOV).ranges{2});
    iTrMid = ismember(glmv5sim2(iFOV).iTrialsUnder,glmv5sim2(iFOV).gainUp + glmv5sim2(iFOV).ranges{3});
    iTrLate = ismember(glmv5sim2(iFOV).iTrialsUnder,glmv5sim2(iFOV).gainUp + glmv5sim2(iFOV).ranges{4});
    iTrWash = ismember(glmv5sim2(iFOV).iTrialsUnder,glmv5sim2(iFOV).gainDown + glmv5sim2(iFOV).ranges{5});    
    nTrialsRew = [sum(iTrBase) sum(iTrEarly) sum(iTrMid) sum(iTrLate) sum(iTrWash)];
    fprintf('Undershoot trials across adaptation: %s \n',num2str(nTrialsRew))    
    
    % Get averages - reward PCs
    sim2data(iFOV).spkUnderBase = mean(glmv5sim2(iFOV).spkUnder(end/2-imfs+1:end/2+imfs,glmv5sim2(iFOV).sigBoolRew,iTrBase),3);
    sim2data(iFOV).spkUnderBaseFR = mean(sim2data(iFOV).spkUnderBase(imfs+1:imfs+.2*imfs,:))';
    sim2data(iFOV).spkUnderBaseMu = mean(sim2data(iFOV).spkUnderBase,2);
    sim2data(iFOV).velUnderBase = mean(glmv5sim2(iFOV).velUnder(end/2-vrfsDec+1:end/2+vrfsDec,iTrBase),2);
    
    sim2data(iFOV).spkUnderEarly = mean(glmv5sim2(iFOV).spkUnder(end/2-imfs+1:end/2+imfs,glmv5sim2(iFOV).sigBoolRew,iTrEarly),3);
    sim2data(iFOV).spkUnderEarlyFR = mean(sim2data(iFOV).spkUnderEarly(imfs+1:imfs+.2*imfs,:))';
    sim2data(iFOV).spkUnderEarlyMu = mean(sim2data(iFOV).spkUnderEarly,2);
    sim2data(iFOV).velUnderEarly = mean(glmv5sim2(iFOV).velUnder(end/2-vrfsDec+1:end/2+vrfsDec,iTrEarly),2);
    
    sim2data(iFOV).spkUnderMid = mean(glmv5sim2(iFOV).spkUnder(end/2-imfs+1:end/2+imfs,glmv5sim2(iFOV).sigBoolRew,iTrMid),3);
    sim2data(iFOV).spkUnderMidFR = mean(sim2data(iFOV).spkUnderMid(imfs+1:imfs+.2*imfs,:))';
    sim2data(iFOV).spkUnderMidMu = mean(sim2data(iFOV).spkUnderMid,2);
    sim2data(iFOV).velUnderMid = mean(glmv5sim2(iFOV).velUnder(end/2-vrfsDec+1:end/2+vrfsDec,iTrMid),2);
    
    sim2data(iFOV).spkUnderLate = mean(glmv5sim2(iFOV).spkUnder(end/2-imfs+1:end/2+imfs,glmv5sim2(iFOV).sigBoolRew,iTrLate),3);
    sim2data(iFOV).spkUnderLateFR = mean(sim2data(iFOV).spkUnderLate(imfs+1:imfs+.2*imfs,:))';
    sim2data(iFOV).spkUnderLateMu = mean(sim2data(iFOV).spkUnderLate,2);
    sim2data(iFOV).velUnderLate = mean(glmv5sim2(iFOV).velUnder(end/2-vrfsDec+1:end/2+vrfsDec,iTrLate),2);
    
    sim2data(iFOV).spkUnderWash = mean(glmv5sim2(iFOV).spkUnder(end/2-imfs+1:end/2+imfs,glmv5sim2(iFOV).sigBoolRew,iTrWash),3); 
    sim2data(iFOV).spkUnderWashFR = mean(sim2data(iFOV).spkUnderWash(imfs+1:imfs+.2*imfs,:))';
    sim2data(iFOV).spkUnderWashMu = mean(sim2data(iFOV).spkUnderWash,2);
    sim2data(iFOV).velUnderWash = mean(glmv5sim2(iFOV).velUnder(end/2-vrfsDec+1:end/2+vrfsDec,iTrWash),2);

    % Get averages - movement PCs
    sim2data(iFOV).spkUnder_movBase = mean(glmv5sim2(iFOV).spkUnder(end/2-imfs+1:end/2+imfs,glmv5sim2(iFOV).sigBoolMov,iTrBase),3);
    sim2data(iFOV).spkUnder_movBaseFR = mean(sim2data(iFOV).spkUnder_movBase(imfs+1:imfs+.2*imfs,:))';
    sim2data(iFOV).spkUnder_movBaseMu = mean(sim2data(iFOV).spkUnder_movBase,2);
    sim2data(iFOV).velUnderBase = mean(glmv5sim2(iFOV).velUnder(end/2-vrfsDec+1:end/2+vrfsDec,iTrBase),2);
    
    sim2data(iFOV).spkUnder_movEarly = mean(glmv5sim2(iFOV).spkUnder(end/2-imfs+1:end/2+imfs,glmv5sim2(iFOV).sigBoolMov,iTrEarly),3);
    sim2data(iFOV).spkUnder_movEarlyFR = mean(sim2data(iFOV).spkUnder_movEarly(imfs+1:imfs+.2*imfs,:))';
    sim2data(iFOV).spkUnder_movEarlyMu = mean(sim2data(iFOV).spkUnder_movEarly,2);
    sim2data(iFOV).velUnderEarly = mean(glmv5sim2(iFOV).velUnder(end/2-vrfsDec+1:end/2+vrfsDec,iTrEarly),2);
    
    sim2data(iFOV).spkUnder_movMid = mean(glmv5sim2(iFOV).spkUnder(end/2-imfs+1:end/2+imfs,glmv5sim2(iFOV).sigBoolMov,iTrMid),3);
    sim2data(iFOV).spkUnder_movMidFR = mean(sim2data(iFOV).spkUnder_movMid(imfs+1:imfs+.2*imfs,:))';
    sim2data(iFOV).spkUnder_movMidMu = mean(sim2data(iFOV).spkUnder_movMid,2);
    sim2data(iFOV).velUnderMid = mean(glmv5sim2(iFOV).velUnder(end/2-vrfsDec+1:end/2+vrfsDec,iTrMid),2);
    
    sim2data(iFOV).spkUnder_movLate = mean(glmv5sim2(iFOV).spkUnder(end/2-imfs+1:end/2+imfs,glmv5sim2(iFOV).sigBoolMov,iTrLate),3);
    sim2data(iFOV).spkUnder_movLateFR = mean(sim2data(iFOV).spkUnder_movLate(imfs+1:imfs+.2*imfs,:))';
    sim2data(iFOV).spkUnder_movLateMu = mean(sim2data(iFOV).spkUnder_movLate,2);
    sim2data(iFOV).velUnderLate = mean(glmv5sim2(iFOV).velUnder(end/2-vrfsDec+1:end/2+vrfsDec,iTrLate),2);
    
    sim2data(iFOV).spkUnder_movWash = mean(glmv5sim2(iFOV).spkUnder(end/2-imfs+1:end/2+imfs,glmv5sim2(iFOV).sigBoolMov,iTrWash),3); 
    sim2data(iFOV).spkUnder_movWashFR = mean(sim2data(iFOV).spkUnder_movWash(imfs+1:imfs+.2*imfs,:))';
    sim2data(iFOV).spkUnder_movWashMu = mean(sim2data(iFOV).spkUnder_movWash,2);
    sim2data(iFOV).velUnderWash = mean(glmv5sim2(iFOV).velUnder(end/2-vrfsDec+1:end/2+vrfsDec,iTrWash),2);

    % Get trials for sim2, overshoot
    iTrBase = ismember(glmv5sim2(iFOV).iTrialsOver,glmv5sim2(iFOV).gainUp + glmv5sim2(iFOV).ranges{1});
    iTrEarly = ismember(glmv5sim2(iFOV).iTrialsOver,glmv5sim2(iFOV).gainUp + glmv5sim2(iFOV).ranges{2});
    iTrMid = ismember(glmv5sim2(iFOV).iTrialsOver,glmv5sim2(iFOV).gainUp + glmv5sim2(iFOV).ranges{3});
    iTrLate = ismember(glmv5sim2(iFOV).iTrialsOver,glmv5sim2(iFOV).gainUp + glmv5sim2(iFOV).ranges{4});
    iTrWash = ismember(glmv5sim2(iFOV).iTrialsOver,glmv5sim2(iFOV).gainDown + glmv5sim2(iFOV).ranges{5});    
    nTrialsRew = [sum(iTrBase) sum(iTrEarly) sum(iTrMid) sum(iTrLate) sum(iTrWash)];
    fprintf('Overshoot trials across adaptation: %s \n',num2str(nTrialsRew))    
    
    % Get averages - rewPCs
    sim2data(iFOV).spkOverBase = mean(glmv5sim2(iFOV).spkOver(end/2-imfs+1:end/2+imfs,glmv5sim2(iFOV).sigBoolRew,iTrBase),3);
    sim2data(iFOV).spkOverBaseFR = mean(sim2data(iFOV).spkOverBase(imfs+1:imfs+.2*imfs,:))';
    sim2data(iFOV).spkOverBaseMu = mean(sim2data(iFOV).spkOverBase,2);
    sim2data(iFOV).velOverBase = mean(glmv5sim2(iFOV).velOver(end/2-vrfsDec+1:end/2+vrfsDec,iTrBase),2);
    
    sim2data(iFOV).spkOverEarly = mean(glmv5sim2(iFOV).spkOver(end/2-imfs+1:end/2+imfs,glmv5sim2(iFOV).sigBoolRew,iTrEarly),3);
    sim2data(iFOV).spkOverEarlyFR = mean(sim2data(iFOV).spkOverEarly(imfs+1:imfs+.2*imfs,:))';
    sim2data(iFOV).spkOverEarlyMu = mean(sim2data(iFOV).spkOverEarly,2);
    sim2data(iFOV).velOverEarly = mean(glmv5sim2(iFOV).velOver(end/2-vrfsDec+1:end/2+vrfsDec,iTrEarly),2);
    
    sim2data(iFOV).spkOverMid = mean(glmv5sim2(iFOV).spkOver(end/2-imfs+1:end/2+imfs,glmv5sim2(iFOV).sigBoolRew,iTrMid),3);
    sim2data(iFOV).spkOverMidFR = mean(sim2data(iFOV).spkOverMid(imfs+1:imfs+.2*imfs,:))';
    sim2data(iFOV).spkOverMidMu = mean(sim2data(iFOV).spkOverMid,2);
    sim2data(iFOV).velOverMid = mean(glmv5sim2(iFOV).velOver(end/2-vrfsDec+1:end/2+vrfsDec,iTrMid),2);
    
    sim2data(iFOV).spkOverLate = mean(glmv5sim2(iFOV).spkOver(end/2-imfs+1:end/2+imfs,glmv5sim2(iFOV).sigBoolRew,iTrLate),3);
    sim2data(iFOV).spkOverLateFR = mean(sim2data(iFOV).spkOverLate(imfs+1:imfs+.2*imfs,:))';
    sim2data(iFOV).spkOverLateMu = mean(sim2data(iFOV).spkOverLate,2);
    sim2data(iFOV).velOverLate = mean(glmv5sim2(iFOV).velOver(end/2-vrfsDec+1:end/2+vrfsDec,iTrLate),2);
    
    sim2data(iFOV).spkOverWash = mean(glmv5sim2(iFOV).spkOver(end/2-imfs+1:end/2+imfs,glmv5sim2(iFOV).sigBoolRew,iTrWash),3); 
    sim2data(iFOV).spkOverWashFR = mean(sim2data(iFOV).spkOverWash(imfs+1:imfs+.2*imfs,:))';
    sim2data(iFOV).spkOverWashMu = mean(sim2data(iFOV).spkOverWash,2);
    sim2data(iFOV).velOverWash = mean(glmv5sim2(iFOV).velOver(end/2-vrfsDec+1:end/2+vrfsDec,iTrWash),2);

    % Get averages - mov PCs
    sim2data(iFOV).spkOver_movBase = mean(glmv5sim2(iFOV).spkOver(end/2-imfs+1:end/2+imfs,glmv5sim2(iFOV).sigBoolMov,iTrBase),3);
    sim2data(iFOV).spkOver_movBaseFR = mean(sim2data(iFOV).spkOver_movBase(imfs+1:imfs+.2*imfs,:))';
    sim2data(iFOV).spkOver_movBaseMu = mean(sim2data(iFOV).spkOver_movBase,2);
    sim2data(iFOV).velOverBase = mean(glmv5sim2(iFOV).velOver(end/2-vrfsDec+1:end/2+vrfsDec,iTrBase),2);
    
    sim2data(iFOV).spkOver_movEarly = mean(glmv5sim2(iFOV).spkOver(end/2-imfs+1:end/2+imfs,glmv5sim2(iFOV).sigBoolMov,iTrEarly),3);
    sim2data(iFOV).spkOver_movEarlyFR = mean(sim2data(iFOV).spkOver_movEarly(imfs+1:imfs+.2*imfs,:))';
    sim2data(iFOV).spkOver_movEarlyMu = mean(sim2data(iFOV).spkOver_movEarly,2);
    sim2data(iFOV).velOverEarly = mean(glmv5sim2(iFOV).velOver(end/2-vrfsDec+1:end/2+vrfsDec,iTrEarly),2);
    
    sim2data(iFOV).spkOver_movMid = mean(glmv5sim2(iFOV).spkOver(end/2-imfs+1:end/2+imfs,glmv5sim2(iFOV).sigBoolMov,iTrMid),3);
    sim2data(iFOV).spkOver_movMidFR = mean(sim2data(iFOV).spkOver_movMid(imfs+1:imfs+.2*imfs,:))';
    sim2data(iFOV).spkOver_movMidMu = mean(sim2data(iFOV).spkOver_movMid,2);
    sim2data(iFOV).velOverMid = mean(glmv5sim2(iFOV).velOver(end/2-vrfsDec+1:end/2+vrfsDec,iTrMid),2);
    
    sim2data(iFOV).spkOver_movLate = mean(glmv5sim2(iFOV).spkOver(end/2-imfs+1:end/2+imfs,glmv5sim2(iFOV).sigBoolMov,iTrLate),3);
    sim2data(iFOV).spkOver_movLateFR = mean(sim2data(iFOV).spkOver_movLate(imfs+1:imfs+.2*imfs,:))';
    sim2data(iFOV).spkOver_movLateMu = mean(sim2data(iFOV).spkOver_movLate,2);
    sim2data(iFOV).velOverLate = mean(glmv5sim2(iFOV).velOver(end/2-vrfsDec+1:end/2+vrfsDec,iTrLate),2);
    
    sim2data(iFOV).spkOver_movWash = mean(glmv5sim2(iFOV).spkOver(end/2-imfs+1:end/2+imfs,glmv5sim2(iFOV).sigBoolMov,iTrWash),3); 
    sim2data(iFOV).spkOver_movWashFR = mean(sim2data(iFOV).spkOver_movWash(imfs+1:imfs+.2*imfs,:))';
    sim2data(iFOV).spkOver_movWashMu = mean(sim2data(iFOV).spkOver_movWash,2);
    sim2data(iFOV).velOverWash = mean(glmv5sim2(iFOV).velOver(end/2-vrfsDec+1:end/2+vrfsDec,iTrWash),2);
end

%% Make averages and plot:
% Lobv mov
lobvmov = imfs*cat(1,mean(cat(2,lobvdata.spkMovBaseMu),2),mean(cat(2,lobvdata.spkMovEarlyMu),2),...
    mean(cat(2,lobvdata.spkMovMidMu),2),mean(cat(2,lobvdata.spkMovLateMu),2),...
    mean(cat(2,lobvdata.spkMovWashMu),2));
lobvmovSEM = imfs*cat(1,std(cat(2,lobvdata.spkMovBaseMu),[],2),std(cat(2,lobvdata.spkMovEarlyMu),[],2),...
    std(cat(2,lobvdata.spkMovMidMu),[],2),std(cat(2,lobvdata.spkMovLateMu),[],2),...
    std(cat(2,lobvdata.spkMovWashMu),[],2))/sqrt(numel(lobvdata));
f1 = figure; 
tplot = (0.5:1:imfs*10)';
h1 = boundedline(tplot,lobvmov,lobvmovSEM,'cmap',[0 0 0],'alpha'); set(h1,'linewidth',1.5);
xlim([0 imfs*10]); xticks([0:imfs:imfs*10]);
ylim([0 12]); yticks([0:2:15]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
title('lobv Movement - mov PCs');

% sim2 mov
sim2mov = imfs*cat(1,mean(cat(2,sim2data.spkMovBaseMu),2),mean(cat(2,sim2data.spkMovEarlyMu),2),...
    mean(cat(2,sim2data.spkMovMidMu),2),mean(cat(2,sim2data.spkMovLateMu),2),...
    mean(cat(2,sim2data.spkMovWashMu),2));
sim2movSEM = imfs*cat(1,std(cat(2,sim2data.spkMovBaseMu),[],2),std(cat(2,sim2data.spkMovEarlyMu),[],2),...
    std(cat(2,sim2data.spkMovMidMu),[],2),std(cat(2,sim2data.spkMovLateMu),[],2),...
    std(cat(2,sim2data.spkMovWashMu),[],2))/sqrt(numel(sim2data));
f2 = figure; 
tplot = (0.5:1:imfs*10)';
h1 = boundedline(tplot,sim2mov,sim2movSEM,'cmap',[0 0 0],'alpha'); set(h1,'linewidth',1.5);
xlim([0 imfs*10]); xticks([0:imfs:imfs*10]);
ylim([0 12]); yticks([0:2:15]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
title('sim2 Movement - mov PCs')

% lobv reward
lobvrew = imfs*cat(1,mean(cat(2,lobvdata.spkRewBaseMu),2),mean(cat(2,lobvdata.spkRewEarlyMu),2),...
    mean(cat(2,lobvdata.spkRewMidMu),2),mean(cat(2,lobvdata.spkRewLateMu),2),...
    mean(cat(2,lobvdata.spkRewWashMu),2));
lobvrewSEM = imfs*cat(1,std(cat(2,lobvdata.spkRewBaseMu),[],2),std(cat(2,lobvdata.spkRewEarlyMu),[],2),...
    std(cat(2,lobvdata.spkRewMidMu),[],2),std(cat(2,lobvdata.spkRewLateMu),[],2),...
    std(cat(2,lobvdata.spkRewWashMu),[],2))/sqrt(numel(lobvdata));
f3 = figure; 
tplot = (0.5:1:imfs*10)';
h1 = boundedline(tplot,lobvrew,lobvrewSEM,'cmap',[0 0 0],'alpha'); set(h1,'linewidth',1.5);
xlim([0 imfs*10]); xticks([0:imfs:imfs*10]);
ylim([0 12]); yticks([0:2:15]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
title('lobv Reward - rew PCs')

% sim2 reward
sim2rew = imfs*cat(1,mean(cat(2,sim2data.spkRewBaseMu),2),mean(cat(2,sim2data.spkRewEarlyMu),2),...
    mean(cat(2,sim2data.spkRewMidMu),2),mean(cat(2,sim2data.spkRewLateMu),2),...
    mean(cat(2,sim2data.spkRewWashMu),2));
sim2rewSEM = imfs*cat(1,std(cat(2,sim2data.spkRewBaseMu),[],2),std(cat(2,sim2data.spkRewEarlyMu),[],2),...
    std(cat(2,sim2data.spkRewMidMu),[],2),std(cat(2,sim2data.spkRewLateMu),[],2),...
    std(cat(2,sim2data.spkRewWashMu),[],2))/sqrt(numel(sim2data));
f4 = figure; 
tplot = (0.5:1:imfs*10)';
h1 = boundedline(tplot,sim2rew,sim2rewSEM,'cmap',[0 0 0],'alpha'); set(h1,'linewidth',1.5);
xlim([0 imfs*10]); xticks([0:30:imfs*10]);
ylim([0 12]); yticks([0:2:15]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
title('sim2 Reward - rew PCs');

% lobv under - rew cells
lobvUnderMu_rew = imfs*cat(1,mean(cat(2,lobvdata.spkUnderBaseMu),2, 'omitnan'),mean(cat(2,lobvdata.spkUnderEarlyMu),2, 'omitnan'),...
    mean(cat(2,lobvdata.spkUnderMidMu),2, 'omitnan'),mean(cat(2,lobvdata.spkUnderLateMu),2, 'omitnan'),...
    mean(cat(2,lobvdata.spkUnderWashMu),2, 'omitnan'));
lobvUnderSEM_rew = imfs*cat(1,std(cat(2,lobvdata.spkUnderBaseMu),[],2, 'omitnan'),std(cat(2,lobvdata.spkUnderEarlyMu),[],2, 'omitnan'),...
    std(cat(2,lobvdata.spkUnderMidMu),[],2, 'omitnan'),std(cat(2,lobvdata.spkUnderLateMu),[],2, 'omitnan'),...
    std(cat(2,lobvdata.spkUnderWashMu),[],2, 'omitnan'))/sqrt(numel(lobvdata));
f5 = figure; 
tplot = (0.5:1:imfs*10)';
h1 = boundedline(tplot,lobvUnderMu_rew,lobvUnderSEM_rew,'cmap',[0 0 0],'alpha'); set(h1,'linewidth',1.5);
xlim([0 imfs*10]); xticks([0:imfs:imfs*10]);
ylim([0 12]); yticks([0:2:15]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
title('lobv Undershoots - rew PCs')

% lobv under - mov cells
lobvUnderMu_mov = imfs*cat(1,mean(cat(2,lobvdata.spkUnder_movBaseMu),2, 'omitnan'),mean(cat(2,lobvdata.spkUnder_movEarlyMu),2, 'omitnan'),...
    mean(cat(2,lobvdata.spkUnder_movMidMu),2, 'omitnan'),mean(cat(2,lobvdata.spkUnder_movLateMu),2, 'omitnan'),...
    mean(cat(2,lobvdata.spkUnder_movWashMu),2, 'omitnan'));
lobvUnderSEM_mov = imfs*cat(1,std(cat(2,lobvdata.spkUnder_movBaseMu),[],2, 'omitnan'),std(cat(2,lobvdata.spkUnder_movEarlyMu),[],2, 'omitnan'),...
    std(cat(2,lobvdata.spkUnder_movMidMu),[],2, 'omitnan'),std(cat(2,lobvdata.spkUnder_movLateMu),[],2, 'omitnan'),...
    std(cat(2,lobvdata.spkUnder_movWashMu),[],2, 'omitnan'))/sqrt(numel(lobvdata));
f6 = figure; 
tplot = (0.5:1:imfs*10)';
h1 = boundedline(tplot,lobvUnderMu_mov,lobvUnderSEM_mov,'cmap',[0 0 0],'alpha'); set(h1,'linewidth',1.5);
xlim([0 imfs*10]); xticks([0:imfs:imfs*10]);
ylim([0 12]); yticks([0:2:15]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
title('lobv Undershoots - mov PCs')

% lobv over - rew cells
lobvOverMu_rew = imfs*cat(1,mean(cat(2,lobvdata.spkOverBaseMu),2, 'omitnan'),mean(cat(2,lobvdata.spkOverEarlyMu),2, 'omitnan'),...
    mean(cat(2,lobvdata.spkOverMidMu),2, 'omitnan'),mean(cat(2,lobvdata.spkOverLateMu),2, 'omitnan'),...
    mean(cat(2,lobvdata.spkOverWashMu),2, 'omitnan'));
lobvOverSEM_rew = imfs*cat(1,std(cat(2,lobvdata.spkOverBaseMu),[],2, 'omitnan'),std(cat(2,lobvdata.spkOverEarlyMu),[],2, 'omitnan'),...
    std(cat(2,lobvdata.spkOverMidMu),[],2, 'omitnan'),std(cat(2,lobvdata.spkOverLateMu),[],2, 'omitnan'),...
    std(cat(2,lobvdata.spkOverWashMu),[],2, 'omitnan'))/sqrt(numel(lobvdata));
f7 = figure; 
tplot = (0.5:1:imfs*10)';
h1 = boundedline(tplot,lobvOverMu_rew,lobvOverSEM_rew,'cmap',[0 0 0],'alpha'); set(h1,'linewidth',1.5);
xlim([0 imfs*10]); xticks([0:imfs:imfs*10]);
ylim([0 12]); yticks([0:2:15]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
title('lobv Overshoots - rew PCs')

% lobv over - mov cells
lobvOverMu_mov = imfs*cat(1,mean(cat(2,lobvdata.spkOver_movBaseMu),2, 'omitnan'),mean(cat(2,lobvdata.spkOver_movEarlyMu),2, 'omitnan'),...
    mean(cat(2,lobvdata.spkOver_movMidMu),2, 'omitnan'),mean(cat(2,lobvdata.spkOver_movLateMu),2, 'omitnan'),...
    mean(cat(2,lobvdata.spkOver_movWashMu),2, 'omitnan'));
lobvOverSEM_mov = imfs*cat(1,std(cat(2,lobvdata.spkOver_movBaseMu),[],2, 'omitnan'),std(cat(2,lobvdata.spkOver_movEarlyMu),[],2, 'omitnan'),...
    std(cat(2,lobvdata.spkOver_movMidMu),[],2, 'omitnan'),std(cat(2,lobvdata.spkOver_movLateMu),[],2, 'omitnan'),...
    std(cat(2,lobvdata.spkOver_movWashMu),[],2, 'omitnan'))/sqrt(numel(lobvdata));
f8 = figure; 
tplot = (0.5:1:imfs*10)';
h1 = boundedline(tplot,lobvOverMu_mov,lobvOverSEM_mov,'cmap',[0 0 0],'alpha'); set(h1,'linewidth',1.5);
xlim([0 imfs*10]); xticks([0:imfs:imfs*10]);
ylim([0 12]); yticks([0:2:15]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
title('lobv Overshoots - mov PCs')

% sim2 under - rew cells
sim2UnderMu_rew = imfs*cat(1,mean(cat(2,sim2data.spkUnderBaseMu),2, 'omitnan'),mean(cat(2,sim2data.spkUnderEarlyMu),2, 'omitnan'),...
    mean(cat(2,sim2data.spkUnderMidMu),2, 'omitnan'),mean(cat(2,sim2data.spkUnderLateMu),2, 'omitnan'),...
    mean(cat(2,sim2data.spkUnderWashMu),2, 'omitnan'));
sim2UnderSEM_rew = imfs*cat(1,std(cat(2,sim2data.spkUnderBaseMu),[],2, 'omitnan'),std(cat(2,sim2data.spkUnderEarlyMu),[],2, 'omitnan'),...
    std(cat(2,sim2data.spkUnderMidMu),[],2, 'omitnan'),std(cat(2,sim2data.spkUnderLateMu),[],2, 'omitnan'),...
    std(cat(2,sim2data.spkUnderWashMu),[],2, 'omitnan'))/sqrt(numel(sim2data));
f9 = figure; 
tplot = (0.5:1:imfs*10)';
h1 = boundedline(tplot,sim2UnderMu_rew,sim2UnderSEM_rew,'cmap',[0 0 0],'alpha'); set(h1,'linewidth',1.5);
xlim([0 imfs*10]); xticks([0:imfs:imfs*10]);
ylim([0 12]); yticks([0:2:15]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
title('sim2 Undershoots - rew PCs')

% sim2 under - mov cells
sim2UnderMu_mov = imfs*cat(1,mean(cat(2,sim2data.spkUnder_movBaseMu),2, 'omitnan'),mean(cat(2,sim2data.spkUnder_movEarlyMu),2, 'omitnan'),...
    mean(cat(2,sim2data.spkUnder_movMidMu),2, 'omitnan'),mean(cat(2,sim2data.spkUnder_movLateMu),2, 'omitnan'),...
    mean(cat(2,sim2data.spkUnder_movWashMu),2, 'omitnan'));
sim2UnderSEM_mov = imfs*cat(1,std(cat(2,sim2data.spkUnder_movBaseMu),[],2, 'omitnan'),std(cat(2,sim2data.spkUnder_movEarlyMu),[],2, 'omitnan'),...
    std(cat(2,sim2data.spkUnder_movMidMu),[],2, 'omitnan'),std(cat(2,sim2data.spkUnder_movLateMu),[],2, 'omitnan'),...
    std(cat(2,sim2data.spkUnder_movWashMu),[],2, 'omitnan'))/sqrt(numel(sim2data));
f10 = figure; 
tplot = (0.5:1:imfs*10)';
h1 = boundedline(tplot,sim2UnderMu_mov,sim2UnderSEM_mov,'cmap',[0 0 0],'alpha'); set(h1,'linewidth',1.5);
xlim([0 imfs*10]); xticks([0:imfs:imfs*10]);
ylim([0 12]); yticks([0:2:15]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
title('sim2 Undershoots - mov PCs')

% sim2 over - rew cells
sim2OverMu_rew = imfs*cat(1,mean(cat(2,sim2data.spkOverBaseMu),2, 'omitnan'),mean(cat(2,sim2data.spkOverEarlyMu),2, 'omitnan'),...
    mean(cat(2,sim2data.spkOverMidMu),2, 'omitnan'),mean(cat(2,sim2data.spkOverLateMu),2, 'omitnan'),...
    mean(cat(2,sim2data.spkOverWashMu),2, 'omitnan'));
sim2OverSEM_rew = imfs*cat(1,std(cat(2,sim2data.spkOverBaseMu),[],2, 'omitnan'),std(cat(2,sim2data.spkOverEarlyMu),[],2, 'omitnan'),...
    std(cat(2,sim2data.spkOverMidMu),[],2, 'omitnan'),std(cat(2,sim2data.spkOverLateMu),[],2, 'omitnan'),...
    std(cat(2,sim2data.spkOverWashMu),[],2, 'omitnan'))/sqrt(numel(sim2data));
f11 = figure; 
tplot = (0.5:1:imfs*10)';
h1 = boundedline(tplot,sim2OverMu_rew,sim2OverSEM_rew,'cmap',[0 0 0],'alpha'); set(h1,'linewidth',1.5);
xlim([0 imfs*10]); xticks([0:imfs:imfs*10]);
ylim([0 12]); yticks([0:2:15]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
title('sim2 Overshoots - rew PCs')

% sim2 over - mov cells
sim2OverMu_mov = imfs*cat(1,mean(cat(2,sim2data.spkOver_movBaseMu),2, 'omitnan'),mean(cat(2,sim2data.spkOver_movEarlyMu),2, 'omitnan'),...
    mean(cat(2,sim2data.spkOver_movMidMu),2, 'omitnan'),mean(cat(2,sim2data.spkOver_movLateMu),2, 'omitnan'),...
    mean(cat(2,sim2data.spkOver_movWashMu),2, 'omitnan'));
sim2OverSEM_mov = imfs*cat(1,std(cat(2,sim2data.spkOver_movBaseMu),[],2, 'omitnan'),std(cat(2,sim2data.spkOver_movEarlyMu),[],2, 'omitnan'),...
    std(cat(2,sim2data.spkOver_movMidMu),[],2, 'omitnan'),std(cat(2,sim2data.spkOver_movLateMu),[],2, 'omitnan'),...
    std(cat(2,sim2data.spkOver_movWashMu),[],2, 'omitnan'))/sqrt(numel(sim2data));
f12 = figure; 
tplot = (0.5:1:imfs*10)';
h1 = boundedline(tplot,sim2OverMu_mov,sim2OverSEM_mov,'cmap',[0 0 0],'alpha'); set(h1,'linewidth',1.5);
xlim([0 imfs*10]); xticks([0:imfs:imfs*10]);
ylim([0 12]); yticks([0:2:15]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
title('sim2 Overshoots - mov PCs')


if savebool
    savFold = fullfile(paperFold,'Figures','Fig5','Sub2_behavSummary');    
    figname = 'lobvmov';
    savefig(f1, fullfile(savFold,[figname,'.fig']));
    saveas(f1,fullfile(savFold,[figname,'.png']));
    print(f1,fullfile(savFold,[figname,'.eps']), '-depsc', '-vector');
    figname = 'sim2mov';
    savefig(f2, fullfile(savFold,[figname,'.fig']));
    saveas(f2,fullfile(savFold,[figname,'.png']));
    print(f2,fullfile(savFold,[figname,'.eps']), '-depsc', '-vector');
    figname = 'lobvrew';
    savefig(f3, fullfile(savFold,[figname,'.fig']));
    saveas(f3,fullfile(savFold,[figname,'.png']));
    print(f3,fullfile(savFold,[figname,'.eps']), '-depsc', '-vector');
    figname = 'sim2rew';
    savefig(f4, fullfile(savFold,[figname,'.fig']));
    saveas(f4,fullfile(savFold,[figname,'.png']));
    print(f4,fullfile(savFold,[figname,'.eps']), '-depsc', '-vector');
    figname = 'lobvUnder-rew';
    savefig(f5, fullfile(savFold,[figname,'.fig']));
    saveas(f5,fullfile(savFold,[figname,'.png']));
    print(f5,fullfile(savFold,[figname,'.eps']), '-depsc', '-vector');
    figname = 'lobvUnder-mov';
    savefig(f6, fullfile(savFold,[figname,'.fig']));
    saveas(f6,fullfile(savFold,[figname,'.png']));
    print(f6,fullfile(savFold,[figname,'.eps']), '-depsc', '-vector');
    figname = 'lobvOver-rew';
    savefig(f7, fullfile(savFold,[figname,'.fig']));
    saveas(f7,fullfile(savFold,[figname,'.png']));
    print(f7,fullfile(savFold,[figname,'.eps']), '-depsc', '-vector');
    figname = 'lobvOver-mov';
    savefig(f8, fullfile(savFold,[figname,'.fig']));
    saveas(f8,fullfile(savFold,[figname,'.png']));
    print(f8,fullfile(savFold,[figname,'.eps']), '-depsc', '-vector');
    figname = 'sim2Under-rew';
    savefig(f9, fullfile(savFold,[figname,'.fig']));
    saveas(f9,fullfile(savFold,[figname,'.png']));
    print(f9,fullfile(savFold,[figname,'.eps']), '-depsc', '-vector');
    figname = 'sim2Under-mov';
    savefig(f10, fullfile(savFold,[figname,'.fig']));
    saveas(f10,fullfile(savFold,[figname,'.png']));
    print(f10,fullfile(savFold,[figname,'.eps']), '-depsc', '-vector');
    figname = 'sim2Over-rew';
    savefig(f11, fullfile(savFold,[figname,'.fig']));
    saveas(f11,fullfile(savFold,[figname,'.png']));
    print(f11,fullfile(savFold,[figname,'.eps']), '-depsc', '-vector');
    figname = 'sim2Over-mov';
    savefig(f12, fullfile(savFold,[figname,'.fig']));
    saveas(f12,fullfile(savFold,[figname,'.png']));
    print(f12,fullfile(savFold,[figname,'.eps']), '-depsc', '-vector');
        
    save(fullfile(savFold,'v5_adaptSummary.mat'),'lobvdata','sim2data');
end

%% Scatters:
% lv - mov
lvMov_base = cat(1,lobvdata.spkMovBaseFR)*imfs;
lvMov_early = cat(1,lobvdata.spkMovEarlyFR)*imfs;
lvMov_mid = cat(1,lobvdata.spkMovMidFR)*imfs;
lvMov_late = cat(1,lobvdata.spkMovLateFR)*imfs;
lvMov_wash = cat(1,lobvdata.spkMovWashFR)*imfs;
% ls - mov
lsMov_base = cat(1,sim2data.spkMovBaseFR)*imfs;
lsMov_early = cat(1,sim2data.spkMovEarlyFR)*imfs;
lsMov_mid = cat(1,sim2data.spkMovMidFR)*imfs;
lsMov_late = cat(1,sim2data.spkMovLateFR)*imfs;
lsMov_wash = cat(1,sim2data.spkMovWashFR)*imfs;
% lv - rew
lvRew_base = cat(1,lobvdata.spkRewBaseFR)*imfs;
lvRew_early = cat(1,lobvdata.spkRewEarlyFR)*imfs;
lvRew_mid = cat(1,lobvdata.spkRewMidFR)*imfs;
lvRew_late = cat(1,lobvdata.spkRewLateFR)*imfs;
lvRew_wash = cat(1,lobvdata.spkRewWashFR)*imfs;
% ls - rew
lsRew_base = cat(1,sim2data.spkRewBaseFR)*imfs;
lsRew_early = cat(1,sim2data.spkRewEarlyFR)*imfs;
lsRew_mid = cat(1,sim2data.spkRewMidFR)*imfs;
lsRew_late = cat(1,sim2data.spkRewLateFR)*imfs;
lsRew_wash = cat(1,sim2data.spkRewWashFR)*imfs;

% lv - under-rew
lvUnderRew_base = cat(1,lobvdata.spkUnderBaseFR)*imfs;
lvUnderRew_early = cat(1,lobvdata.spkUnderEarlyFR)*imfs;
lvUnderRew_mid = cat(1,lobvdata.spkUnderMidFR)*imfs;
lvUnderRew_late = cat(1,lobvdata.spkUnderLateFR)*imfs;
lvUnderRew_wash = cat(1,lobvdata.spkUnderWashFR)*imfs;
% lv - under-mov
lvUnderMov_base = cat(1,lobvdata.spkUnder_movBaseFR)*imfs;
lvUnderMov_early = cat(1,lobvdata.spkUnder_movEarlyFR)*imfs;
lvUnderMov_mid = cat(1,lobvdata.spkUnder_movMidFR)*imfs;
lvUnderMov_late = cat(1,lobvdata.spkUnder_movLateFR)*imfs;
lvUnderMov_wash = cat(1,lobvdata.spkUnder_movWashFR)*imfs;
% lv - overRew
lvOverRew_base = cat(1,lobvdata.spkOverBaseFR)*imfs;
lvOverRew_early = cat(1,lobvdata.spkOverEarlyFR)*imfs;
lvOverRew_mid = cat(1,lobvdata.spkOverMidFR)*imfs;
lvOverRew_late = cat(1,lobvdata.spkOverLateFR)*imfs;
lvOverRew_wash = cat(1,lobvdata.spkOverWashFR)*imfs;
% lv - over-mov
lvOverMov_base = cat(1,lobvdata.spkOver_movBaseFR)*imfs;
lvOverMov_early = cat(1,lobvdata.spkOver_movEarlyFR)*imfs;
lvOverMov_mid = cat(1,lobvdata.spkOver_movMidFR)*imfs;
lvOverMov_late = cat(1,lobvdata.spkOver_movLateFR)*imfs;
lvOverMov_wash = cat(1,lobvdata.spkOver_movWashFR)*imfs;

% ls - under-rew
lsUnderRew_base = cat(1,sim2data.spkUnderBaseFR)*imfs;
lsUnderRew_early = cat(1,sim2data.spkUnderEarlyFR)*imfs;
lsUnderRew_mid = cat(1,sim2data.spkUnderMidFR)*imfs;
lsUnderRew_late = cat(1,sim2data.spkUnderLateFR)*imfs;
lsUnderRew_wash = cat(1,sim2data.spkUnderWashFR)*imfs;
% ls - under-mov
lsUnderMov_base = cat(1,sim2data.spkUnder_movBaseFR)*imfs;
lsUnderMov_early = cat(1,sim2data.spkUnder_movEarlyFR)*imfs;
lsUnderMov_mid = cat(1,sim2data.spkUnder_movMidFR)*imfs;
lsUnderMov_late = cat(1,sim2data.spkUnder_movLateFR)*imfs;
lsUnderMov_wash = cat(1,sim2data.spkUnder_movWashFR)*imfs;
% ls - overRew
lsOverRew_base = cat(1,sim2data.spkOverBaseFR)*imfs;
lsOverRew_early = cat(1,sim2data.spkOverEarlyFR)*imfs;
lsOverRew_mid = cat(1,sim2data.spkOverMidFR)*imfs;
lsOverRew_late = cat(1,sim2data.spkOverLateFR)*imfs;
lsOverRew_wash = cat(1,sim2data.spkOverWashFR)*imfs;
% ls - over-mov
lsOverMov_base = cat(1,sim2data.spkOver_movBaseFR)*imfs;
lsOverMov_early = cat(1,sim2data.spkOver_movEarlyFR)*imfs;
lsOverMov_mid = cat(1,sim2data.spkOver_movMidFR)*imfs;
lsOverMov_late = cat(1,sim2data.spkOver_movLateFR)*imfs;
lsOverMov_wash = cat(1,sim2data.spkOver_movWashFR)*imfs;

% Movement - plot
fscat1 = figure;
subplot(2,4,1); hold on;
scatter(lvMov_base,lvMov_early,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel('Lobule V - movement')
title('Base vs Early')

subplot(2,4,2); hold on;
scatter(lvMov_base,lvMov_mid,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Mid')

subplot(2,4,3); hold on;
scatter(lvMov_base,lvMov_late,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Late')

subplot(2,4,4); hold on;
scatter(lvMov_base,lvMov_wash,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Wash')

subplot(2,4,5); hold on;
scatter(lsMov_base,lsMov_early,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Early')
ylabel('Lobule Simplex - Movement')

subplot(2,4,6); hold on;
scatter(lsMov_base,lsMov_mid,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Mid')

subplot(2,4,7); hold on;
scatter(lsMov_base,lsMov_late,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Late')

subplot(2,4,8); hold on;
scatter(lsMov_base,lsMov_wash,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Wash')

% Reward - plot
fscat2 = figure;
subplot(2,4,1); hold on;
scatter(lvRew_base,lvRew_early,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel('Lobule V - reward')
title('Base vs Early')

subplot(2,4,2); hold on;
scatter(lvRew_base,lvRew_mid,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Mid')

subplot(2,4,3); hold on;
scatter(lvRew_base,lvRew_late,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Late')

subplot(2,4,4); hold on;
scatter(lvRew_base,lvRew_wash,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Wash')

subplot(2,4,5); hold on;
scatter(lsRew_base,lsRew_early,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Early')
ylabel('Lobule Simplex - Reward')

subplot(2,4,6); hold on;
scatter(lsRew_base,lsRew_mid,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Mid')

subplot(2,4,7); hold on;
scatter(lsRew_base,lsRew_late,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Late')

subplot(2,4,8); hold on;
scatter(lsRew_base,lsRew_wash,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Wash')

% Undershoot-rew - plot
fscat3 = figure;
subplot(2,4,1); hold on;
scatter(lvUnderRew_base,lvUnderRew_early,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel('Lobule V - Undershoots (rew PCs)')
title('Base vs Early')

subplot(2,4,2); hold on;
scatter(lvUnderRew_base,lvUnderRew_mid,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Mid')

subplot(2,4,3); hold on;
scatter(lvUnderRew_base,lvUnderRew_late,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Late')

subplot(2,4,4); hold on;
scatter(lvUnderRew_base,lvUnderRew_wash,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Wash')

subplot(2,4,5); hold on;
scatter(lsUnderRew_base,lsUnderRew_early,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Early')
ylabel('Lobule simplex - Undershoots (rew PCs)')

subplot(2,4,6); hold on;
scatter(lsUnderRew_base,lsUnderRew_mid,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Mid')

subplot(2,4,7); hold on;
scatter(lsUnderRew_base,lsUnderRew_late,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Late')

subplot(2,4,8); hold on;
scatter(lsUnderRew_base,lsUnderRew_wash,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Wash')

% Undershoot-mov - plot
fscat4 = figure;
subplot(2,4,1); hold on;
scatter(lvUnderMov_base,lvUnderMov_early,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel('Lobule V - Undershoots (mov PCs)')
title('Base vs Early')

subplot(2,4,2); hold on;
scatter(lvUnderMov_base,lvUnderMov_mid,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Mid')

subplot(2,4,3); hold on;
scatter(lvUnderMov_base,lvUnderMov_late,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Late')

subplot(2,4,4); hold on;
scatter(lvUnderMov_base,lvUnderMov_wash,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Wash')

subplot(2,4,5); hold on;
scatter(lsUnderMov_base,lsUnderMov_early,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Early')
ylabel('Lobule simplex - Undershoots (mov PCs)')

subplot(2,4,6); hold on;
scatter(lsUnderMov_base,lsUnderMov_mid,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Mid')

subplot(2,4,7); hold on;
scatter(lsUnderMov_base,lsUnderMov_late,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Late')

subplot(2,4,8); hold on;
scatter(lsUnderMov_base,lsUnderMov_wash,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Wash')

% Overshoot-rew - plot
fscat5 = figure;
subplot(2,4,1); hold on;
scatter(lvOverRew_base,lvOverRew_early,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel('Lobule V - Overshoots (rew PCs)')
title('Base vs Early')

subplot(2,4,2); hold on;
scatter(lvOverRew_base,lvOverRew_mid,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Mid')

subplot(2,4,3); hold on;
scatter(lvOverRew_base,lvOverRew_late,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Late')

subplot(2,4,4); hold on;
scatter(lvOverRew_base,lvOverRew_wash,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Wash')

subplot(2,4,5); hold on;
scatter(lsOverRew_base,lsOverRew_early,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Early')
ylabel('Lobule simplex - Overshoots (rew PCs)')

subplot(2,4,6); hold on;
scatter(lsOverRew_base,lsOverRew_mid,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Mid')

subplot(2,4,7); hold on;
scatter(lsOverRew_base,lsOverRew_late,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Late')

subplot(2,4,8); hold on;
scatter(lsOverRew_base,lsOverRew_wash,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Wash')

% Overshoot-mov - plot
fscat6 = figure;
subplot(2,4,1); hold on;
scatter(lvOverMov_base,lvOverMov_early,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel('Lobule V - Overshoots (mov PCs)')
title('Base vs Early')

subplot(2,4,2); hold on;
scatter(lvOverMov_base,lvOverMov_mid,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Mid')

subplot(2,4,3); hold on;
scatter(lvOverMov_base,lvOverMov_late,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Late')

subplot(2,4,4); hold on;
scatter(lvOverMov_base,lvOverMov_wash,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Wash')

subplot(2,4,5); hold on;
scatter(lsOverMov_base,lsOverMov_early,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Early')
ylabel('Lobule simplex - Overshoots (mov PCs)')

subplot(2,4,6); hold on;
scatter(lsOverMov_base,lsOverMov_mid,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Mid')

subplot(2,4,7); hold on;
scatter(lsOverMov_base,lsOverMov_late,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Late')

subplot(2,4,8); hold on;
scatter(lsOverMov_base,lsOverMov_wash,'filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
xlim([0,15]); ylim([0,15]); xticks([0:3:15]); yticks([0:3:15]); axis square;
hline = refline(1,0);
hline.Color = 'k'; hline.LineStyle = '--'; 
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Base vs Wash')

if savebool    
    savFold = fullfile(paperFold,'Figures','Fig5','Sub2_behavSummary');
    figname = 'movScatter';
    savefig(fscat1, fullfile(savFold,[figname,'.fig']));
    saveas(fscat1,fullfile(savFold,[figname,'.png']));
    print(fscat1,fullfile(savFold,[figname,'.eps']), '-depsc', '-vector');
    figname = 'rewScatter';
    savefig(fscat2, fullfile(savFold,[figname,'.fig']));
    saveas(fscat2,fullfile(savFold,[figname,'.png']));
    print(fscat2,fullfile(savFold,[figname,'.eps']), '-depsc', '-vector'); 
    figname = 'underScatter-rew';
    savefig(fscat3, fullfile(savFold,[figname,'.fig']));
    saveas(fscat3,fullfile(savFold,[figname,'.png']));
    print(fscat3,fullfile(savFold,[figname,'.eps']), '-depsc', '-vector'); 
    figname = 'underScatter-mov';
    savefig(fscat4, fullfile(savFold,[figname,'.fig']));
    saveas(fscat4,fullfile(savFold,[figname,'.png']));
    print(fscat4,fullfile(savFold,[figname,'.eps']), '-depsc', '-vector'); 
    figname = 'overScatter-rew';
    savefig(fscat5, fullfile(savFold,[figname,'.fig']));
    saveas(fscat5,fullfile(savFold,[figname,'.png']));
    print(fscat5,fullfile(savFold,[figname,'.eps']), '-depsc', '-vector'); 
    figname = 'overScatter-mov';
    savefig(fscat6, fullfile(savFold,[figname,'.fig']));
    saveas(fscat6,fullfile(savFold,[figname,'.png']));
    print(fscat6,fullfile(savFold,[figname,'.eps']), '-depsc', '-vector'); 
end