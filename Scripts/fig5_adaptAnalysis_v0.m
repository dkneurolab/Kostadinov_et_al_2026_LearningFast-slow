%% Settings:
% Decimation:
nDec = 50;

%% Load in stuff we need:
path45 = 'E:\Analysis_current\paperFigs\Fig1\Sub3_behavFigs';
cd(path45);
load('matchDataBehav_old.mat')
matchData = rmfield(matchData,{'on1','on2','rewC1','rewC2','rewR1','rewR2'});
    
pathSig = 'C:\Users\Dimitar2\Dropbox\DK papers\2022_Kostadinov_FastSlowLearning\Figures\Fig3\Sub3_SummaryBars\spk_gaussian_alpha50';
cd(pathSig);
load('sigSummary.mat');

pathv5 = 'C:\Users\Dimitar2\Dropbox\DK papers\2022_Kostadinov_FastSlowLearning\Figures\Fig1\Sub1_percMovCorr';
cd(pathv5);
load('percMovCorrect_v5.mat');

%% Trim

seshv5lobv = sessionsv5([10,13,15,17,19]);
seshv5sim2 = sessionsv5([11,12,14,16,18]);

glmv5lobv0 = glmSummary([1,4,6,8,12]);
glmv5sim20 = glmSummary([2,3,5,7,11]);

% matchV5lobv = matchData([3,8,14,17,27]);
% matchV5sim2 = matchData([4,7,11,20,30]);
matchV5lobv = matchData([4,10,18,23,37]);
matchV5sim2 = matchData([5,9,15,26,40]);
clear glmSummary sessionsv5 pathv5 pathSig path145

%% Go through and dig out correct trials
dataPath0 = 'E:\Analysis_current\01_Behav+imaging\Version_5';

% First LobV
glmv5lobv = struct;
for iFOV = 1:numel(glmv5lobv0)
    TS = load(fullfile(dataPath0,seshv5lobv(iFOV).name,[seshv5lobv(iFOV).date,'_',seshv5lobv(iFOV).fov],'trialsonly.mat'));
    roisv4 = matchV5lobv(iFOV).iROIs1;
    roisv5 = matchV5lobv(iFOV).iROIs2;
    
    % Take only good ROIs
    glmv5lobv(iFOV).sigBool = glmv5lobv0(iFOV).sigBool(roisv4);
    glmv5lobv(iFOV).sigBoolVis = glmv5lobv0(iFOV).sigBoolVis(roisv4);
    glmv5lobv(iFOV).sigBoolMov = glmv5lobv0(iFOV).sigBoolMov(roisv4);
    glmv5lobv(iFOV).sigBoolRew = glmv5lobv0(iFOV).sigBoolRew(roisv4);
    glmv5lobv(iFOV).sigBoolLick = glmv5lobv0(iFOV).sigBoolLick(roisv4);
    glmv5lobv(iFOV).gainUp = seshv5lobv(iFOV).sets.changetrial;
    glmv5lobv(iFOV).gainDown = seshv5lobv(iFOV).sets.changebacktrial;
    
    glmv5lobv(iFOV).iTrials = cat(1,TS.trialstructs.LCstruct.trialind); 
    spkMov = cat(3,TS.trialstructs.LCstruct.spkLmov);
    xMov = cat(3,TS.trialstructs.LCstruct.VRLmov); xMov = squeeze(xMov(:,3,:));
    spkRew = cat(3,TS.trialstructs.LCstruct.spkLoff);
    xRew = cat(3,TS.trialstructs.LCstruct.VRLoff); xRew = squeeze(xRew(:,3,:));
    
    % Initialize settings and arrays for behaviour:
    trialSets = TS.sets;
    nDec = 50;
    gainDiv = 37.5/360*4.5*pi/trialSets.vrgain;
    smspan = trialSets.fs*.06; % 60 ms smoothing span
    vrfsDec = trialSets.fs/nDec;
    vMov = zeros(size(xMov,1)/50,size(xMov,2),2);
    vRew = vMov;
    
    % Iterate through movements - should be the same for mov and rew
    for jMov = 1:size(xMov,2)
        x0 = xMov(:,jMov);
        x1 = smooth(x0,smspan/numel(x0),'lowess');
        v1 = diff(x1)*trialSets.fs*-gainDiv; v1(end+1) = v1(end);
        x2 = -(x1-mean(x1((trialSets.imtime-1)*trialSets.fs+1:(trialSets.imtime-1)*trialSets.fs+trialSets.fs/2)));
        vMov(:,jMov,:) = [decimate(x1,nDec) decimate(v1,nDec)];
        
        x0 = xRew(:,jMov);
        x1 = smooth(x0,smspan/numel(x0),'lowess');
        v1 = diff(x1)*trialSets.fs*-gainDiv; v1(end+1) = v1(end);
        x2 = -(x1-mean(x1((trialSets.imtime-1)*trialSets.fs+1:(trialSets.imtime-1)*trialSets.fs+trialSets.fs/2)));
        vRew(:,jMov,:) = [decimate(x1,nDec) decimate(v1,nDec)];
    end
    
    glmv5lobv(iFOV).spkMov = spkMov(:,roisv5,:);
    glmv5lobv(iFOV).posMov = vMov(:,:,1);
    glmv5lobv(iFOV).velMov = vMov(:,:,2);
    glmv5lobv(iFOV).spkCor = spkRew(:,roisv5,:);
    glmv5lobv(iFOV).posCor = vRew(:,:,1);
    glmv5lobv(iFOV).velCor = vRew(:,:,2);
end

% Now do sim2
glmv5sim2 = struct;
for iFOV = 1:numel(glmv5sim20)
    TS = load(fullfile(dataPath0,seshv5sim2(iFOV).name,[seshv5sim2(iFOV).date,'_',seshv5sim2(iFOV).fov],'trialsonly.mat'));
    roisv4 = matchV5sim2(iFOV).iROIs1;
    roisv5 = matchV5sim2(iFOV).iROIs2;
    
    % Take only good ROIs
    glmv5sim2(iFOV).sigBool = glmv5sim20(iFOV).sigBool(roisv4);
    glmv5sim2(iFOV).sigBoolVis = glmv5sim20(iFOV).sigBoolVis(roisv4);
    glmv5sim2(iFOV).sigBoolMov = glmv5sim20(iFOV).sigBoolMov(roisv4);
    glmv5sim2(iFOV).sigBoolRew = glmv5sim20(iFOV).sigBoolRew(roisv4);
    glmv5sim2(iFOV).sigBoolLick = glmv5sim20(iFOV).sigBoolLick(roisv4);
    glmv5sim2(iFOV).gainUp = seshv5sim2(iFOV).sets.changetrial;
    glmv5sim2(iFOV).gainDown = seshv5sim2(iFOV).sets.changebacktrial;
    
    glmv5sim2(iFOV).iTrials = cat(1,TS.trialstructs.LCstruct.trialind); 
    spkMov = cat(3,TS.trialstructs.LCstruct.spkLmov);
    xMov = cat(3,TS.trialstructs.LCstruct.VRLmov); xMov = squeeze(xMov(:,3,:));
    spkRew = cat(3,TS.trialstructs.LCstruct.spkLoff);
    xRew = cat(3,TS.trialstructs.LCstruct.VRLoff); xRew = squeeze(xRew(:,3,:));
    
    % Initialize settings and arrays for behaviour:
    trialSets = TS.sets;
    nDec = 50;
    gainDiv = 37.5/360*4.5*pi/trialSets.vrgain;
    smspan = trialSets.fs*.06; % 60 ms smoothing span
    vrfsDec = trialSets.fs/nDec;
    vMov = zeros(size(xMov,1)/nDec,size(xMov,2),2);
    vRew = vMov;
    
    % Iterate through movements - should be the same for mov and rew
    for jMov = 1:size(xMov,2)
        x0 = xMov(:,jMov);
        x1 = smooth(x0,smspan/numel(x0),'lowess');
        v1 = diff(x1)*trialSets.fs*-gainDiv; v1(end+1) = v1(end);
        x2 = -(x1-mean(x1((trialSets.imtime-1)*trialSets.fs+1:(trialSets.imtime-1)*trialSets.fs+trialSets.fs/2)));
        vMov(:,jMov,:) = [decimate(x1,nDec) decimate(v1,nDec)];
        
        x0 = xRew(:,jMov);
        x1 = smooth(x0,smspan/numel(x0),'lowess');
        v1 = diff(x1)*trialSets.fs*-gainDiv; v1(end+1) = v1(end);
        x2 = -(x1-mean(x1((trialSets.imtime-1)*trialSets.fs+1:(trialSets.imtime-1)*trialSets.fs+trialSets.fs/2)));
        vRew(:,jMov,:) = [decimate(x1,nDec) decimate(v1,nDec)];
    end
    
    glmv5sim2(iFOV).spkMov = spkMov(:,roisv5,:);
    glmv5sim2(iFOV).posMov = vMov(:,:,1);
    glmv5sim2(iFOV).velMov = vMov(:,:,2);
    glmv5sim2(iFOV).spkCor = spkRew(:,roisv5,:);
    glmv5sim2(iFOV).posCor = vRew(:,:,1);
    glmv5sim2(iFOV).velCor = vRew(:,:,2);
end

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
behavAdapt = load(fullfile('C:\Users\Dimitar2\Dropbox\DK papers\2022_Kostadinov_FastSlowLearning\Figures\Fig1\Sub1_v5Adaptation',...
    'vAdaptationCurves.mat'));
iBehavLobv = [10,13,15,18,21];
iBehavSim2 = [11,12,14,17,20];
for iBehav = 1:numel(iBehavLobv)
    glmv5lobv(iBehav).ranges{1,1} = behavAdapt.sessionsv5(iBehavLobv(iBehav)).endPos.iTrialsBase - ...
        glmv5lobv(iBehav).gainUp;
    glmv5lobv(iBehav).ranges{2,1} = behavAdapt.sessionsv5(iBehavLobv(iBehav)).endPos.iTrialsEarly(1:end/2) - ...
        glmv5lobv(iBehav).gainUp;
    glmv5lobv(iBehav).ranges{3,1} = behavAdapt.sessionsv5(iBehavLobv(iBehav)).endPos.iTrialsMid - ...
        glmv5lobv(iBehav).gainUp;
    glmv5lobv(iBehav).ranges{4,1} = behavAdapt.sessionsv5(iBehavLobv(iBehav)).endPos.iTrialsLate - ...
        glmv5lobv(iBehav).gainUp;
    glmv5lobv(iBehav).ranges{5,1} = behavAdapt.sessionsv5(iBehavLobv(iBehav)).endPos.iTrialsWash - ...
        glmv5lobv(iBehav).gainDown;
end

for iBehav = 1:numel(iBehavSim2)
    glmv5sim2(iBehav).ranges{1,1} = behavAdapt.sessionsv5(iBehavSim2(iBehav)).endPos.iTrialsBase - ...
        glmv5sim2(iBehav).gainUp;
    glmv5sim2(iBehav).ranges{2,1} = behavAdapt.sessionsv5(iBehavSim2(iBehav)).endPos.iTrialsEarly(1:end/2) - ...
        glmv5sim2(iBehav).gainUp;
    glmv5sim2(iBehav).ranges{3,1} = behavAdapt.sessionsv5(iBehavSim2(iBehav)).endPos.iTrialsMid - ...
        glmv5sim2(iBehav).gainUp;
    glmv5sim2(iBehav).ranges{4,1} = behavAdapt.sessionsv5(iBehavSim2(iBehav)).endPos.iTrialsLate - ...
        glmv5sim2(iBehav).gainUp;
    glmv5sim2(iBehav).ranges{5,1} = behavAdapt.sessionsv5(iBehavSim2(iBehav)).endPos.iTrialsWash - ...
        glmv5sim2(iBehav).gainDown;
end

%% Plot 15 trial rolling average for mov trials and reward trials:
figFold = 'C:\Users\Dimitar2\Dropbox\DK papers\2022_Kostadinov_FastSlowLearning\Figures\Fig5\Sub1_FOVexamples';
if ~exist(figFold,'dir'); mkdir(figFold); end
trialGrp = 30;
for iPlot = 1:numel(glmv5lobv)
    [figOut,dataOut] = fig5_adaptPlotFOV(glmv5lobv(iPlot),trialGrp);
    if savebool
        figname = [sprintf('%s_lobv_movAdapt',lobvdata(iPlot).name)];
        savefig(figOut(1), fullfile(figFold,[figname,'.fig']));
        saveas(figOut(1),fullfile(figFold,[figname,'.png']));
        print(figOut(1),fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
        figname = [sprintf('%s_lobv_rewAdapt',lobvdata(iPlot).name)];
        savefig(figOut(2), fullfile(figFold,[figname,'.fig']));
        saveas(figOut(2),fullfile(figFold,[figname,'.png']));
        print(figOut(2),fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    end    
end

for iPlot = 1:numel(glmv5sim2)
    [figOut,dataOut] = fig5_adaptPlotFOV(glmv5sim2(iPlot),trialGrp);
    if savebool
        figname = [sprintf('%s_sim2_movAdapt',sim2data(iPlot).name)];
        savefig(figOut(1), fullfile(figFold,[figname,'.fig']));
        saveas(figOut(1),fullfile(figFold,[figname,'.png']));
        print(figOut(1),fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
        figname = [sprintf('%s_sim2_rewAdapt',sim2data(iPlot).name)];
        savefig(figOut(2), fullfile(figFold,[figname,'.fig']));
        saveas(figOut(2),fullfile(figFold,[figname,'.png']));
        print(figOut(2),fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    end    
end


%%
% rangeBase = [-29:0];
% rangeEarly = [1:15];
% rangeMid = [46:75];
% rangeLate = [96:125];
% rangeWash = [1:15];

% rangeBase = [-59:0];
% rangeEarly = [1:20];
% rangeMid = [51:70];
% rangeLate = [101:120];
% rangeWash = [1:10];


lobvdata = struct;
for iFOV = 1:numel(glmv5lobv)
    lobvdata(iFOV).name = seshv5lobv(iFOV).name;
    iTrBase = ismember(glmv5lobv(iFOV).iTrials,glmv5lobv(iFOV).gainUp + glmv5lobv(iFOV).ranges{1});
    iTrEarly = ismember(glmv5lobv(iFOV).iTrials,glmv5lobv(iFOV).gainUp + glmv5lobv(iFOV).ranges{2});
    iTrMid = ismember(glmv5lobv(iFOV).iTrials,glmv5lobv(iFOV).gainUp + glmv5lobv(iFOV).ranges{3});
    iTrLate = ismember(glmv5lobv(iFOV).iTrials,glmv5lobv(iFOV).gainUp + glmv5lobv(iFOV).ranges{4});
    iTrWash = ismember(glmv5lobv(iFOV).iTrials,glmv5lobv(iFOV).gainDown + glmv5lobv(iFOV).ranges{5});
    
    nTrials = [sum(iTrBase) sum(iTrEarly) sum(iTrMid) sum(iTrLate) sum(iTrWash)]
    
    lobvdata(iFOV).spkMovBase = mean(glmv5lobv(iFOV).spkMov(31:90,glmv5lobv(iFOV).sigBoolMov,iTrBase),3);
    lobvdata(iFOV).spkMovBaseMu = mean(lobvdata(iFOV).spkMovBase,2);
    lobvdata(iFOV).velMovBase = mean(glmv5lobv(iFOV).velMov(101:300,iTrBase),2);
    lobvdata(iFOV).spkRewBase = mean(glmv5lobv(iFOV).spkCor(31:90,glmv5lobv(iFOV).sigBoolRew,iTrBase),3);
    lobvdata(iFOV).spkRewBaseMu = mean(lobvdata(iFOV).spkRewBase,2);
    lobvdata(iFOV).velRewBase = mean(glmv5lobv(iFOV).velCor(101:300,iTrBase),2);
    
    lobvdata(iFOV).spkMovEarly = mean(glmv5lobv(iFOV).spkMov(31:90,glmv5lobv(iFOV).sigBoolMov,iTrEarly),3);
    lobvdata(iFOV).spkMovEarlyMu = mean(lobvdata(iFOV).spkMovEarly,2);
    lobvdata(iFOV).velMovEarly = mean(glmv5lobv(iFOV).velMov(101:300,iTrEarly),2);
    lobvdata(iFOV).spkRewEarly = mean(glmv5lobv(iFOV).spkCor(31:90,glmv5lobv(iFOV).sigBoolRew,iTrEarly),3);
    lobvdata(iFOV).spkRewEarlyMu = mean(lobvdata(iFOV).spkRewEarly,2);
    lobvdata(iFOV).velRewEarly = mean(glmv5lobv(iFOV).velCor(101:300,iTrEarly),2);
    
    lobvdata(iFOV).spkMovMid = mean(glmv5lobv(iFOV).spkMov(31:90,glmv5lobv(iFOV).sigBoolMov,iTrMid),3);
    lobvdata(iFOV).spkMovMidMu = mean(lobvdata(iFOV).spkMovMid,2);
    lobvdata(iFOV).velMovMid = mean(glmv5lobv(iFOV).velMov(101:300,iTrMid),2);
    lobvdata(iFOV).spkRewMid = mean(glmv5lobv(iFOV).spkCor(31:90,glmv5lobv(iFOV).sigBoolRew,iTrMid),3);
    lobvdata(iFOV).spkRewMidMu = mean(lobvdata(iFOV).spkRewMid,2);
    lobvdata(iFOV).velRewMid = mean(glmv5lobv(iFOV).velCor(101:300,iTrMid),2);
    
    lobvdata(iFOV).spkMovLate = mean(glmv5lobv(iFOV).spkMov(31:90,glmv5lobv(iFOV).sigBoolMov,iTrLate),3);
    lobvdata(iFOV).spkMovLateMu = mean(lobvdata(iFOV).spkMovLate,2);
    lobvdata(iFOV).velMovLate = mean(glmv5lobv(iFOV).velMov(101:300,iTrLate),2);
    lobvdata(iFOV).spkRewLate = mean(glmv5lobv(iFOV).spkCor(31:90,glmv5lobv(iFOV).sigBoolRew,iTrLate),3);
    lobvdata(iFOV).spkRewLateMu = mean(lobvdata(iFOV).spkRewLate,2);
    lobvdata(iFOV).velRewLate = mean(glmv5lobv(iFOV).velCor(101:300,iTrLate),2);
    
    lobvdata(iFOV).spkMovWash = mean(glmv5lobv(iFOV).spkMov(31:90,glmv5lobv(iFOV).sigBoolMov,iTrWash),3);
    lobvdata(iFOV).spkMovWashMu = mean(lobvdata(iFOV).spkMovWash,2);
    lobvdata(iFOV).velMovWash = mean(glmv5lobv(iFOV).velMov(101:300,iTrWash),2);
    lobvdata(iFOV).spkRewWash = mean(glmv5lobv(iFOV).spkCor(31:90,glmv5lobv(iFOV).sigBoolRew,iTrWash),3); 
    lobvdata(iFOV).spkRewWashMu = mean(lobvdata(iFOV).spkRewWash,2);
    lobvdata(iFOV).velRewWash = mean(glmv5lobv(iFOV).velCor(101:300,iTrWash),2);
    
end

sim2data = struct;
for iFOV = 1:numel(glmv5sim2)
    sim2data(iFOV).name = seshv5sim2(iFOV).name;
    iTrBase = ismember(glmv5sim2(iFOV).iTrials,glmv5sim2(iFOV).gainUp + glmv5sim2(iFOV).ranges{1});
    iTrEarly = ismember(glmv5sim2(iFOV).iTrials,glmv5sim2(iFOV).gainUp + glmv5sim2(iFOV).ranges{2});
    iTrMid = ismember(glmv5sim2(iFOV).iTrials,glmv5sim2(iFOV).gainUp + glmv5sim2(iFOV).ranges{3});
    iTrLate = ismember(glmv5sim2(iFOV).iTrials,glmv5sim2(iFOV).gainUp + glmv5sim2(iFOV).ranges{4});
    iTrWash = ismember(glmv5sim2(iFOV).iTrials,glmv5sim2(iFOV).gainDown + glmv5sim2(iFOV).ranges{5});
    
    nTrialss2 = ([sum(iTrBase) sum(iTrEarly) sum(iTrMid) sum(iTrLate) sum(iTrWash)])
    
    sim2data(iFOV).spkMovBase = mean(glmv5sim2(iFOV).spkMov(31:90,glmv5sim2(iFOV).sigBoolMov,iTrBase),3);
    sim2data(iFOV).spkMovBaseMu = mean(sim2data(iFOV).spkMovBase,2);
    sim2data(iFOV).velMovBase = mean(glmv5sim2(iFOV).velMov(101:300,iTrBase),2);
    sim2data(iFOV).spkRewBase = mean(glmv5sim2(iFOV).spkCor(31:90,glmv5sim2(iFOV).sigBoolRew,iTrBase),3);
    sim2data(iFOV).spkRewBaseMu = mean(sim2data(iFOV).spkRewBase,2);
    sim2data(iFOV).velRewBase = mean(glmv5sim2(iFOV).velCor(101:300,iTrBase),2);
    
    sim2data(iFOV).spkMovEarly = mean(glmv5sim2(iFOV).spkMov(31:90,glmv5sim2(iFOV).sigBoolMov,iTrEarly),3);
    sim2data(iFOV).spkMovEarlyMu = mean(sim2data(iFOV).spkMovEarly,2);
    sim2data(iFOV).velMovEarly = mean(glmv5sim2(iFOV).velMov(101:300,iTrEarly),2);
    sim2data(iFOV).spkRewEarly = mean(glmv5sim2(iFOV).spkCor(31:90,glmv5sim2(iFOV).sigBoolRew,iTrEarly),3);
    sim2data(iFOV).spkRewEarlyMu = mean(sim2data(iFOV).spkRewEarly,2);
    sim2data(iFOV).velRewEarly = mean(glmv5sim2(iFOV).velCor(101:300,iTrEarly),2);
    
    sim2data(iFOV).spkMovMid = mean(glmv5sim2(iFOV).spkMov(31:90,glmv5sim2(iFOV).sigBoolMov,iTrMid),3);
    sim2data(iFOV).spkMovMidMu = mean(sim2data(iFOV).spkMovMid,2);
    sim2data(iFOV).velMovMid = mean(glmv5sim2(iFOV).velMov(101:300,iTrMid),2);
    sim2data(iFOV).spkRewMid = mean(glmv5sim2(iFOV).spkCor(31:90,glmv5sim2(iFOV).sigBoolRew,iTrMid),3);
    sim2data(iFOV).spkRewMidMu = mean(sim2data(iFOV).spkRewMid,2);
    sim2data(iFOV).velRewMid = mean(glmv5sim2(iFOV).velCor(101:300,iTrMid),2);
    
    sim2data(iFOV).spkMovLate = mean(glmv5sim2(iFOV).spkMov(31:90,glmv5sim2(iFOV).sigBoolMov,iTrLate),3);
    sim2data(iFOV).spkMovLateMu = mean(sim2data(iFOV).spkMovLate,2);
    sim2data(iFOV).velMovLate = mean(glmv5sim2(iFOV).velMov(101:300,iTrLate),2);
    sim2data(iFOV).spkRewLate = mean(glmv5sim2(iFOV).spkCor(31:90,glmv5sim2(iFOV).sigBoolRew,iTrLate),3);
    sim2data(iFOV).spkRewLateMu = mean(sim2data(iFOV).spkRewLate,2);
    sim2data(iFOV).velRewLate = mean(glmv5sim2(iFOV).velCor(101:300,iTrLate),2);
    
    sim2data(iFOV).spkMovWash = mean(glmv5sim2(iFOV).spkMov(31:90,glmv5sim2(iFOV).sigBoolMov,iTrWash),3);
    sim2data(iFOV).spkMovWashMu = mean(sim2data(iFOV).spkMovWash,2);
    sim2data(iFOV).velMovWash = mean(glmv5sim2(iFOV).velMov(101:300,iTrWash),2);
    sim2data(iFOV).spkRewWash = mean(glmv5sim2(iFOV).spkCor(31:90,glmv5sim2(iFOV).sigBoolRew,iTrWash),3); 
    sim2data(iFOV).spkRewWashMu = mean(sim2data(iFOV).spkRewWash,2);
    sim2data(iFOV).velRewWash = mean(glmv5sim2(iFOV).velCor(101:300,iTrWash),2);
    
end

%
lobvmov = 30*cat(1,mean(cat(2,lobvdata.spkMovBaseMu),2),mean(cat(2,lobvdata.spkMovEarlyMu),2),...
    mean(cat(2,lobvdata.spkMovMidMu),2),mean(cat(2,lobvdata.spkMovLateMu),2),...
    mean(cat(2,lobvdata.spkMovWashMu),2));
lobvmovSEM = 30*cat(1,std(cat(2,lobvdata.spkMovBaseMu),[],2),std(cat(2,lobvdata.spkMovEarlyMu),[],2),...
    std(cat(2,lobvdata.spkMovMidMu),[],2),std(cat(2,lobvdata.spkMovLateMu),[],2),...
    std(cat(2,lobvdata.spkMovWashMu),[],2))/sqrt(numel(lobvdata));
f1 = figure; 
tplot = (0.5:1:300)';
h1 = boundedline(tplot,lobvmov,lobvmovSEM,'cmap',[0 0 0],'alpha'); set(h1,'linewidth',1.5);
xlim([0 300]); xticks([0:30:300]);
ylim([0 12]); yticks([0:2:15]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
title('lobv mov');

sim2mov = 30*cat(1,mean(cat(2,sim2data.spkMovBaseMu),2),mean(cat(2,sim2data.spkMovEarlyMu),2),...
    mean(cat(2,sim2data.spkMovMidMu),2),mean(cat(2,sim2data.spkMovLateMu),2),...
    mean(cat(2,sim2data.spkMovWashMu),2));
sim2movSEM = 30*cat(1,std(cat(2,sim2data.spkMovBaseMu),[],2),std(cat(2,sim2data.spkMovEarlyMu),[],2),...
    std(cat(2,sim2data.spkMovMidMu),[],2),std(cat(2,sim2data.spkMovLateMu),[],2),...
    std(cat(2,sim2data.spkMovWashMu),[],2))/sqrt(numel(sim2data));
f2 = figure; 
tplot = (0.5:1:300)';
h1 = boundedline(tplot,sim2mov,sim2movSEM,'cmap',[0 0 0],'alpha'); set(h1,'linewidth',1.5);
xlim([0 300]); xticks([0:30:300]);
ylim([0 12]); yticks([0:2:15]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
title('sim2 mov')

lobvrew = 30*cat(1,mean(cat(2,lobvdata.spkRewBaseMu),2),mean(cat(2,lobvdata.spkRewEarlyMu),2),...
    mean(cat(2,lobvdata.spkRewMidMu),2),mean(cat(2,lobvdata.spkRewLateMu),2),...
    mean(cat(2,lobvdata.spkRewWashMu),2));
lobvrewSEM = 30*cat(1,std(cat(2,lobvdata.spkRewBaseMu),[],2),std(cat(2,lobvdata.spkRewEarlyMu),[],2),...
    std(cat(2,lobvdata.spkRewMidMu),[],2),std(cat(2,lobvdata.spkRewLateMu),[],2),...
    std(cat(2,lobvdata.spkRewWashMu),[],2))/sqrt(numel(lobvdata));
f3 = figure; 
tplot = (0.5:1:300)';
h1 = boundedline(tplot,lobvrew,lobvrewSEM,'cmap',[0 0 0],'alpha'); set(h1,'linewidth',1.5);
xlim([0 300]); xticks([0:30:300]);
ylim([0 12]); yticks([0:2:15]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
title('lobv rew')

sim2rew = 30*cat(1,mean(cat(2,sim2data.spkRewBaseMu),2),mean(cat(2,sim2data.spkRewEarlyMu),2),...
    mean(cat(2,sim2data.spkRewMidMu),2),mean(cat(2,sim2data.spkRewLateMu),2),...
    mean(cat(2,sim2data.spkRewWashMu),2));
sim2rewSEM = 30*cat(1,std(cat(2,sim2data.spkRewBaseMu),[],2),std(cat(2,sim2data.spkRewEarlyMu),[],2),...
    std(cat(2,sim2data.spkRewMidMu),[],2),std(cat(2,sim2data.spkRewLateMu),[],2),...
    std(cat(2,sim2data.spkRewWashMu),[],2))/sqrt(numel(sim2data));

f4 = figure; 
tplot = (0.5:1:300)';
h1 = boundedline(tplot,sim2rew,sim2rewSEM,'cmap',[0 0 0],'alpha'); set(h1,'linewidth',1.5);
xlim([0 300]); xticks([0:30:300]);
ylim([0 12]); yticks([0:2:15]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
title('sim2 rew');

if savebool
savFold = 'C:\Users\Dimitar2\Dropbox\DK papers\2022_Kostadinov_FastSlowLearning\Figures\Fig5\Sub2_behavSummary';
figname = 'lobvmov';
savefig(f1, fullfile(savFold,[figname,'.fig']));
saveas(f1,fullfile(savFold,[figname,'.png']));
print(f1,fullfile(savFold,[figname,'.eps']), '-depsc', '-painters');
figname = 'sim2mov';
savefig(f2, fullfile(savFold,[figname,'.fig']));
saveas(f2,fullfile(savFold,[figname,'.png']));
print(f2,fullfile(savFold,[figname,'.eps']), '-depsc', '-painters');
figname = 'lobvrew';
savefig(f3, fullfile(savFold,[figname,'.fig']));
saveas(f3,fullfile(savFold,[figname,'.png']));
print(f3,fullfile(savFold,[figname,'.eps']), '-depsc', '-painters');
figname = 'sim2rew';
savefig(f4, fullfile(savFold,[figname,'.fig']));
saveas(f4,fullfile(savFold,[figname,'.png']));
print(f4,fullfile(savFold,[figname,'.eps']), '-depsc', '-painters');
end

