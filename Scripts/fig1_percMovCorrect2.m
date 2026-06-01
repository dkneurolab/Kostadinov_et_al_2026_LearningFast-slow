function fig1_percMovCorrect2(datasets,dataFold,paperFold,adaptTrialGrp,matchColors,closebool,savebool)
%% Concatenate sessions and
mouseStructs = fieldnames(datasets);
sessions = [];
for iSesh = 1:numel(mouseStructs)
    sessions = [sessions datasets.(mouseStructs{iSesh})]; %#ok<AGROW>
end 
seshV1 = sessions(cat(1,sessions.version) == 1);
seshV4 = sessions(cat(1,sessions.version) == 4);
sessionsv5 = sessions(cat(1,sessions.version) == 5);

behavFold = fullfile(dataFold,'04_Behav-processed');
figFold = fullfile(paperFold,'Figures','Fig1','Sub1_percMovCorr');
figFold2 = fullfile(dataFold,'paperFigs','Fig1','Sub1_percMovCorr');

%% Version 1 section
if ~exist(fullfile(figFold,'percMovCorrect_v1.mat'),'file')
    % If this has not already been done, iterate over session and extract 
    % percent correct and fraction of trials with movement

for iSesh = 1:numel(seshV1)
    
    % Load in data
    behavData = load(fullfile(dataFold,'04_Behav-processed',sesh(iSesh).name,sesh(iSesh).date,'Behaviour_4s',sprintf('%s_%s.mat',seshV1(iSesh).name,seshV1(iSesh).date)));
    behavData = behavData.datastruct;
    trials = behavData.paqdata.behav.trials;
    inits = behavData.inits;
    sets = behavData.sets;
    gaindiv = 37.5/360*4.5*pi/sets.rotgain;
    if sets.vrfs == 60
        smspan = inits.paqfs*.04; % 40 ms smoothing span
    elseif sets. vrfs == 30
        smspan = inits.paqfs*.06; % 60 ms smoothing span
    end
    
    seshV1(iSesh).nTrials = numel(trials);
    seshV1(iSesh).iTrialsReal = false(seshV1(iSesh).nTrials,1);
    seshV1(iSesh).movBool = false(seshV1(iSesh).nTrials,1);
    seshV1(iSesh).corrBool = false(seshV1(iSesh).nTrials,1);
    
    % Determine if movement happened and if it resulted in correct trial
    for j = 1:numel(trials)
        x0 = trials(j).data(inits.paqfs*inits.mintime+1:end-inits.paqfs*inits.mintime,10);
        x1 = smooth(x0,smspan/numel(x0),'lowess');
        v1 = diff(x1)*inits.paqfs*-gaindiv; v1(end+1) = v1(end);
        if numel(v1) > 15*inits.paqfs; v1 = v1(1:15*inits.paqfs); end
        if numel(v1)/inits.paqfs <= sets.responsetime
            seshV1(iSesh).iTrialsReal(j) = true;
        end
        if any(abs(v1) > 1.5) % movement threshold is 1.5 cm/s mm of wheel translation
%         if any(posTrace > sets.rotgain*360/45/pi()) % movement threshold is 1 mm of wheel translation
            seshV1(iSesh).movBool(j) = true;
        end
        if seshV1(iSesh).iTrialsReal(j) & trials(j).correct
            seshV1(iSesh).corrBool(j) = true;
        end
    end
    seshV1(iSesh).movPerc = sum(seshV1(iSesh).movBool)/seshV1(iSesh).nTrials;
    seshV1(iSesh).nTrialsReal = sum(seshV1(iSesh).iTrialsReal);
    if seshV1(iSesh).nTrialsReal > 0
        seshV1(iSesh).corrPerc = sum(seshV1(iSesh).corrBool)/seshV1(iSesh).nTrialsReal;
    else
        seshV1(iSesh).corrPerc = 0;
    end
    seshV1(iSesh).sets = sets;
    
    switch lower(seshV1(iSesh).name)
        case 'dk103'
            seshV1(iSesh).matchColor = matchColors(1,:);
        case 'dk105'
            seshV1(iSesh).matchColor = matchColors(2,:);
        case 'dk169'
            seshV1(iSesh).matchColor = matchColors(3,:);
        case 'dk171'
            seshV1(iSesh).matchColor = matchColors(4,:);
        case 'dk194'
            seshV1(iSesh).matchColor = matchColors(5,:);
        case 'dk199'
            seshV1(iSesh).matchColor = matchColors(6,:);
        case 'dk052'
            seshV1(iSesh).matchColor = matchColors(7,:);
        case 'dk063'
            seshV1(iSesh).matchColor = matchColors(8,:);
        case 'dk070'
            seshV1(iSesh).matchColor = matchColors(9,:);
        case 'dk138'
            seshV1(iSesh).matchColor = matchColors(10,:);
    end

end

if savebool
    if ~exist(figFold,'dir'); mkdir(figFold); end
    save(fullfile(figFold,'percMovCorrect_v1.mat'),'seshV1');
    if ~exist(figFold2,'dir'); mkdir(figFold2); end
    save(fullfile(figFold2,'percMovCorrect_v1.mat'),'seshV1');
end

else
   % If this has already been done, just load the file to save time
    seshV1 = load(fullfile(figFold,'percMovCorrect_v1.mat'));
    seshV1 = seshV1.sessionsv1;
end 

%% Version 4 section
if ~exist(fullfile(figFold,'percMovCorrect_v4.mat'),'file')
% If this has not already been done, iterate over session and extract 
% percent correct and fraction of trials with movement
for iSesh = 1:numel(seshV4)
    
    % Load in data
    behavData = load(fullfile(dataFold,'04_Behav-processed',sesh(iSesh).name,sesh(iSesh).date,'Behaviour_4s',sprintf('%s_%s.mat',seshV1(iSesh).name,seshV1(iSesh).date)));
    behavData = behavData.datastruct;
    trials = behavData.paqdata.behav.trials;
    inits = behavData.inits;
    sets = behavData.sets;
    gaindiv = 37.5/360*4.5*pi/sets.rotgain;
    if sets.vrfs == 60
        smspan = inits.paqfs*.04; % 40 ms smoothing span
    elseif sets. vrfs == 30
        smspan = inits.paqfs*.06; % 60 ms smoothing span
    end
    if strcmpi(seshV4(iSesh).name,'dk194') & strcmp(seshV4(iSesh).date,'2020-10-02')
        trials = trials(1:180); % Account for lack of motivation at end of this session
    end
    
    % Remove non-responsive trials at end
    trials = trials(1:find(cat(1,trials.realtrial) > 0,1, 'last'));
    behavstability = Rot2p_behavstability(trials,inits,sets,0);
    
    % Initiate 
    seshV4(iSesh).nTrials = numel(trials);
    seshV4(iSesh).nTrialsReal = sum(cat(1,trials.realtrial));
    seshV4(iSesh).iTrialsReal = false(seshV4(iSesh).nTrials,1);
    seshV4(iSesh).movBool = false(seshV4(iSesh).nTrials,1);
    seshV4(iSesh).corrBool = false(seshV4(iSesh).nTrials,1);
    
    % Determine if movement happened and if it resulted in correct trial
    for j = 1:numel(trials)
        x0 = trials(j).data(inits.paqfs*inits.mintime+1:end-inits.paqfs*inits.mintime,10);
        x1 = smooth(x0,smspan/numel(x0),'lowess');
        v1 = diff(x1)*inits.paqfs*-gaindiv; v1(end+1) = v1(end);
        if numel(v1)/inits.paqfs <= sets.responsetime
            seshV4(iSesh).iTrialsReal(j) = true;
        end
        if numel(v1) > 15*inits.paqfs; v1 = v1(1:15*inits.paqfs); end
        if any(abs(v1) > 1.5) % movement threshold is 1.5 cm/s of wheel translation
%         if any(posTrace > sets.rotgain*360/45/pi()) % movement threshold is 1 mm of wheel translation
            seshV4(iSesh).movBool(j) = true;
        end
        if seshV4(iSesh).iTrialsReal(j) & trials(j).correct
            seshV4(iSesh).corrBool(j) = true;
        end
    end
    seshV4(iSesh).movPerc = sum(seshV4(iSesh).movBool)/seshV4(iSesh).nTrials;
%     sessionsv4(i).nTrialsReal = sum(sessionsv4(i).iTrialsReal);
    if seshV4(iSesh).nTrialsReal > 0
        seshV4(iSesh).corrPerc = sum(seshV4(iSesh).corrBool)/seshV4(iSesh).nTrialsReal;
        seshV4(iSesh).xEndMu = mean(behavstability.params.xendall);
        seshV4(iSesh).xEndMed = median(behavstability.params.xendall);
    else
        seshV4(iSesh).corrPerc = 0;
    end
    
    seshV4(iSesh).sets = sets;
    seshV4(iSesh).params = rmfield(behavstability.params,{'trialnum','blocksize','rhoU','rhoC','rhoUmu','rhoUsd','rhoCmu','rhoCsd','rhoOmu','rhoOsd','numC'});
    seshV4(iSesh).offXall = behavstability.offxall;
    seshV4(iSesh).movVall = behavstability.movvall;
    


    % Decide whether session is v4 day 1 or v4 dayN (trained) - logic is
    % that if the mouse and fov names are the same as that of the
    % subsequent session, it is day 1, otherwise it is day N AND the last
    % session should always be day N because there are no day 1 sessions
    % without an day N session
%     if i < numel(sessionsv4)
%         if strcmp(sessionsv4(i).name,sessionsv4(i+1).name) & strcmp(sessionsv4(i).fov,sessionsv4(i+1).fov)
%             sessionsv4(i).day1 = true;
%             sessionsv4(i).dayN = false;
%         else
%             sessionsv4(i).day1 = false;
%             sessionsv4(i).dayN = true;
%         end
%     else
%         sessionsv4(end).day1 = false;
%         sessionsv4(end).dayN = true;
%     end
    
    % Alternative version - first session should always be day 1, others
    % should not mach if they are day 1
    if iSesh > 1 
        if strcmp(seshV4(iSesh).name,seshV4(iSesh-1).name) & strcmp(seshV4(iSesh).fov,seshV4(iSesh-1).fov)
            seshV4(iSesh).day1 = false;
            seshV4(iSesh).dayN = true;
        else
            seshV4(iSesh).day1 = true;
            seshV4(iSesh).dayN = false;
        end
    else
        seshV4(iSesh).day1 = true;
        seshV4(iSesh).dayN = false;
    end
    
    
    switch lower(seshV4(iSesh).name)
        case 'dk103'
            seshV4(iSesh).matchColor = matchColors(1,:);
        case 'dk105'
            seshV4(iSesh).matchColor = matchColors(2,:);
        case 'dk169'
            seshV4(iSesh).matchColor = matchColors(3,:);
        case 'dk171'
            seshV4(iSesh).matchColor = matchColors(4,:);
        case 'dk194'
            seshV4(iSesh).matchColor = matchColors(5,:);
        case 'dk199'
            seshV4(iSesh).matchColor = matchColors(6,:);
        case 'dk052'
            seshV4(iSesh).matchColor = matchColors(7,:);
        case 'dk063'
            seshV4(iSesh).matchColor = matchColors(8,:);
        case 'dk070'
            seshV4(iSesh).matchColor = matchColors(9,:);
        case 'dk138'
            seshV4(iSesh).matchColor = matchColors(10,:);
    end
end


if savebool
    if ~exist(figFold,'dir'); mkdir(figFold); end
    save(fullfile(figFold,'percMovCorrect_v4.mat'),'sessionsv4');
    if ~exist(figFold2,'dir'); mkdir(figFold2); end
    save(fullfile(figFold2,'percMovCorrect_v4.mat'),'sessionsv4');
end
    
else
    % If this has already been done, just load the file to save time
    seshV4 = load(fullfile(figFold,'percMovCorrect_v4.mat'));
    seshV4 = seshV4.sessionsv4;
end

%% Version 5 section
if ~exist(fullfile(figFold,'percMovCorrect_v5.mat'),'file')
% if ~exist(fullfile(figFold,['percMovCorrect_v5_',num2str(adaptTrialGrp),'trials.mat']),'file')
% If this has not already been done, iterate over session and extract 
% percent correct and fraction of trials with movement
for iSesh = 1:numel(sessionsv5)
    dataDir = dir(fullfile(behavFold,sessionsv5(iSesh).name,sessionsv5(iSesh).date,'Behav*',[sessionsv5(iSesh).name,'_',sessionsv5(iSesh).date,'.mat']));
    behavData = load(fullfile(dataDir.folder,dataDir.name)); behavData = behavData.datastruct;
    trials = behavData.paqdata.behav.trials;
    inits = behavData.inits;
    sets = Rot2p_changetrials(behavData.sets);
    gaindiv = 37.5/360*4.5*pi/sets.rotgain;
    if sets.vrfs == 60
        smspan = inits.paqfs*.04; % 40 ms smoothing span
    elseif sets. vrfs == 30
        smspan = inits.paqfs*.06; % 60 ms smoothing span
    end
    if strcmpi(sessionsv5(iSesh).name,'dk105') & strcmp(sessionsv5(iSesh).date,'2018-03-07')
        sets.changebacktrial = sets.changetrial+10*(sets.changebackmult-1);
    end
    trials = trials(1:find(cat(1,trials.realtrial) > 0,1, 'last'));
    
    behavstability = Rot2p_behavstability(trials,inits,sets,0);
    
    sessionsv5(iSesh).nTrials = numel(trials);
    sessionsv5(iSesh).nTrialsReal = sum(cat(1,trials.realtrial));
    sessionsv5(iSesh).nTrialsBase = sum(cat(1,trials(1:sets.changetrial).realtrial));
    sessionsv5(iSesh).nTrialsEarly = sum(cat(1,trials(sets.changetrial+1:sets.changetrial+adaptTrialGrp).realtrial));
    sessionsv5(iSesh).nTrialsLate = sum(cat(1,trials(sets.changetrial+121-adaptTrialGrp:sets.changetrial+120).realtrial));
    
    if ~sets.changeback
        sessionsv5(iSesh).nTrialsWash = 0;
    else
        nWashTrials = numel(trials)-sets.changebacktrial;
        if nWashTrials > adaptTrialGrp; nWashTrials = adaptTrialGrp; end
        sessionsv5(iSesh).nTrialsWash = sum(cat(1,trials(sets.changebacktrial+1:sets.changebacktrial+nWashTrials).realtrial));
    end
    
    sessionsv5(iSesh).iTrialsReal = false(sessionsv5(iSesh).nTrials,1);
    sessionsv5(iSesh).movBool = false(sessionsv5(iSesh).nTrials,1);
    sessionsv5(iSesh).corrBool = false(sessionsv5(iSesh).nTrials,1);

    % Determine if movement happened and if it resulted in correct trial
    for j = 1:numel(trials)
        x0 = trials(j).data(inits.paqfs*inits.mintime+1:end-inits.paqfs*inits.mintime,10);
        x1 = smooth(x0,smspan/numel(x0),'lowess');
        v1 = diff(x1)*inits.paqfs*-gaindiv; v1(end+1) = v1(end);
        if numel(v1) > 15*inits.paqfs; v1 = v1(1:15*inits.paqfs); end
        if numel(v1)/inits.paqfs <= sets.responsetime
            sessionsv5(iSesh).iTrialsReal(j) = true;
        end
        if any(abs(v1) > 1.5) % movement threshold is 1.5 cm/s of wheel translation
%         if any(posTrace > sets.rotgain*360/45/pi()) % movement threshold is 1 mm of wheel translation
            sessionsv5(iSesh).movBool(j) = true;
        end
        if sessionsv5(iSesh).iTrialsReal(j) & trials(j).correct
            sessionsv5(iSesh).corrBool(j) = true;
        end
    end
    sessionsv5(iSesh).movPerc = sum(sessionsv5(iSesh).movBool)/sessionsv5(iSesh).nTrials;
    if sessionsv5(iSesh).nTrialsReal > 0
        sessionsv5(iSesh).corrPerc = sum(sessionsv5(iSesh).corrBool)/sessionsv5(iSesh).nTrialsReal;
        sessionsv5(iSesh).corrPercBase = sum(sessionsv5(iSesh).corrBool(1:sets.changetrial))/sessionsv5(iSesh).nTrialsBase;
        sessionsv5(iSesh).corrPercEarly = sum(sessionsv5(iSesh).corrBool(sets.changetrial+1:sets.changetrial+adaptTrialGrp))/sessionsv5(iSesh).nTrialsEarly;
        sessionsv5(iSesh).corrPercLate = sum(sessionsv5(iSesh).corrBool(sets.changetrial+121-adaptTrialGrp:sets.changetrial+120))/sessionsv5(iSesh).nTrialsLate;
        if ~sets.changeback
            sessionsv5(iSesh).corrPercWash = nan;
        else
            sessionsv5(iSesh).corrPercWash = sum(sessionsv5(iSesh).corrBool(sets.changebacktrial+1:sets.changebacktrial+nWashTrials))/sessionsv5(iSesh).nTrialsWash;
        end
        xEndAll = nan(size(sessionsv5(iSesh).iTrialsReal)); xEndAll(behavstability.params.trialinds) = behavstability.params.xendall;
        
        sessionsv5(iSesh).xEndMu_Base = nanmean(xEndAll(1:sets.changetrial));
        sessionsv5(iSesh).xEndMed_Base = nanmedian(xEndAll(1:sets.changetrial));
        sessionsv5(iSesh).xEndMu_Early = nanmean(xEndAll(sets.changetrial+1:sets.changetrial+adaptTrialGrp));
        sessionsv5(iSesh).xEndMed_Early = nanmedian(xEndAll(sets.changetrial+1:sets.changetrial+adaptTrialGrp));
        sessionsv5(iSesh).xEndMu_Late = nanmean(xEndAll(sets.changetrial+121-adaptTrialGrp:sets.changetrial+120));
        sessionsv5(iSesh).xEndMed_Late = nanmedian(xEndAll(sets.changetrial+121-adaptTrialGrp:sets.changetrial+120));
        if ~sets.changeback
            sessionsv5(iSesh).xEndMu_Wash = nan;
            sessionsv5(iSesh).xEndMed_Wash = nan;
        else
            sessionsv5(iSesh).xEndMu_Wash = nanmean(xEndAll(sets.changebacktrial+1:sets.changebacktrial+nWashTrials));
            sessionsv5(iSesh).xEndMed_Wash = nanmedian(xEndAll(sets.changebacktrial+1:sets.changebacktrial+nWashTrials));
        end
    else
        sessionsv5(iSesh).corrPerc = 0;
    end
    sessionsv5(iSesh).sets = sets;
    sessionsv5(iSesh).params = rmfield(behavstability.params,{'trialnum','blocksize','rhoU','rhoC','rhoUmu','rhoUsd','rhoCmu','rhoCsd','rhoOmu','rhoOsd','numC'});
    sessionsv5(iSesh).offXall = behavstability.offxall;
    sessionsv5(iSesh).movVall = behavstability.movvall;
    
    switch lower(sessionsv5(iSesh).name)
        case 'dk103'
            sessionsv5(iSesh).matchColor = matchColors(1,:);
        case 'dk105'
            sessionsv5(iSesh).matchColor = matchColors(2,:);
        case 'dk169'
            sessionsv5(iSesh).matchColor = matchColors(3,:);
        case 'dk171'
            sessionsv5(iSesh).matchColor = matchColors(4,:);
        case 'dk194'
            sessionsv5(iSesh).matchColor = matchColors(5,:);
        case 'dk199'
            sessionsv5(iSesh).matchColor = matchColors(6,:);
        case 'dk052'
            sessionsv5(iSesh).matchColor = matchColors(7,:);
        case 'dk063'
            sessionsv5(iSesh).matchColor = matchColors(8,:);
        case 'dk070'
            sessionsv5(iSesh).matchColor = matchColors(9,:);
        case 'dk138'
            sessionsv5(iSesh).matchColor = matchColors(10,:);
    end
end


if savebool
    if ~exist(figFold,'dir'); mkdir(figFold); end
    save(fullfile(figFold,'percMovCorrect_v5.mat'),'sessionsv5');
%     save(fullfile(figFold,['percMovCorrect_v5_',num2str(adaptTrialGrp),'trials.mat']),'sessionsv5');
    if ~exist(figFold2,'dir'); mkdir(figFold2); end
    save(fullfile(figFold2,'percMovCorrect_v5.mat'),'sessionsv5');
%     save(fullfile(figFold2,['percMovCorrect_v5_',num2str(adaptTrialGrp),'trials.mat']),'sessionsv5');
end
    
else
    % If this has already been done, just load the file to save time
    sessionsv5 = load(fullfile(figFold,'percMovCorrect_v5.mat'),'sessionsv5');
%     sessionsv5 = load(fullfile(figFold,['percMovCorrect_v5_',num2str(adaptTrialGrp),'trials.mat']));
    sessionsv5 = sessionsv5.sessionsv5;
    
    for iSesh = 1:numel(sessionsv5)
        sets = sessionsv5(iSesh).sets;
        
        sessionsv5(iSesh).nTrialsBase = sum(sessionsv5(iSesh).iTrialsReal(1:sets.changetrial));
        sessionsv5(iSesh).nTrialsEarly = sum(sessionsv5(iSesh).iTrialsReal(sets.changetrial+1:sets.changetrial+adaptTrialGrp));
        sessionsv5(iSesh).nTrialsLate = sum(sessionsv5(iSesh).iTrialsReal(sets.changetrial+121-adaptTrialGrp:sets.changetrial+120));
%         sessionsv5(i).nTrialsBase = sum(cat(1,trials(1:sets.changetrial).realtrial));
%         sessionsv5(i).nTrialsEarly = sum(cat(1,trials(sets.changetrial+1:sets.changetrial+adaptTrialGrp).realtrial));
%         sessionsv5(i).nTrialsLate = sum(cat(1,trials(sets.changetrial+121-adaptTrialGrp:sets.changetrial+120).realtrial));
        
        if ~sets.changeback
            sessionsv5(iSesh).nTrialsWash = 0;
        else
            nWashTrials = sessionsv5(iSesh).nTrials-sets.changebacktrial;
            if nWashTrials > adaptTrialGrp; nWashTrials = adaptTrialGrp; end
            sessionsv5(iSesh).nTrialsWash = sum(sessionsv5(iSesh).iTrialsReal(sets.changebacktrial+1:sets.changebacktrial+nWashTrials));
%             sessionsv5(i).nTrialsWash = sum(cat(1,trials(sets.changebacktrial+1:sets.changebacktrial+nWashTrials).realtrial));
        end
        
        if sessionsv5(iSesh).nTrialsReal > 0
            sessionsv5(iSesh).corrPerc = sum(sessionsv5(iSesh).corrBool)/sessionsv5(iSesh).nTrialsReal;
            sessionsv5(iSesh).corrPercBase = sum(sessionsv5(iSesh).corrBool(1:sets.changetrial))/sessionsv5(iSesh).nTrialsBase;
            sessionsv5(iSesh).corrPercEarly = sum(sessionsv5(iSesh).corrBool(sets.changetrial+1:sets.changetrial+adaptTrialGrp))/sessionsv5(iSesh).nTrialsEarly;
            sessionsv5(iSesh).corrPercLate = sum(sessionsv5(iSesh).corrBool(sets.changetrial+121-adaptTrialGrp:sets.changetrial+120))/sessionsv5(iSesh).nTrialsLate;
            if ~sets.changeback
                sessionsv5(iSesh).corrPercWash = nan;
            else
                sessionsv5(iSesh).corrPercWash = sum(sessionsv5(iSesh).corrBool(sets.changebacktrial+1:sets.changebacktrial+nWashTrials))/sessionsv5(iSesh).nTrialsWash;
            end
            xEndAll = nan(size(sessionsv5(iSesh).iTrialsReal)); xEndAll(sessionsv5(iSesh).params.trialinds) = sessionsv5(iSesh).params.xendall;
            
            sessionsv5(iSesh).xEndMu_Base = nanmean(xEndAll(1:sets.changetrial));
            sessionsv5(iSesh).xEndMed_Base = nanmedian(xEndAll(1:sets.changetrial));
            sessionsv5(iSesh).xEndMu_Early = nanmean(xEndAll(sets.changetrial+1:sets.changetrial+adaptTrialGrp));
            sessionsv5(iSesh).xEndMed_Early = nanmedian(xEndAll(sets.changetrial+1:sets.changetrial+adaptTrialGrp));
            sessionsv5(iSesh).xEndMu_Late = nanmean(xEndAll(sets.changetrial+121-adaptTrialGrp:sets.changetrial+120));
            sessionsv5(iSesh).xEndMed_Late = nanmedian(xEndAll(sets.changetrial+121-adaptTrialGrp:sets.changetrial+120));
            
            if ~sets.changeback
                sessionsv5(iSesh).xEndMu_Wash = nan;
                sessionsv5(iSesh).xEndMed_Wash = nan;
            else
                sessionsv5(iSesh).xEndMu_Wash = nanmean(xEndAll(sets.changebacktrial+1:sets.changebacktrial+nWashTrials));
                sessionsv5(iSesh).xEndMed_Wash = nanmedian(xEndAll(sets.changebacktrial+1:sets.changebacktrial+nWashTrials));
            end
            
        else
            sessionsv5(iSesh).corrPerc = 0;
        end
    end
    
end

clear datastruct gaindiv i inits j nWashTrials sets smspan trials v1 x0 x1 xEndAll

%% Make plots

sessionsv4_1 = seshV4(cat(1,seshV4.day1));
sessionsv4_n = seshV4(cat(1,seshV4.dayN));

% First figure: fraction of trials with movement
figout(1) = figure;
subplot(2,2,1); hold on;
smov(1) = scatter(rand(size(seshV1))/5+0.9,cat(1,seshV1.movPerc),72,cat(1,seshV1.matchColor)/2,'filled');
smov(1).MarkerEdgeColor = 'none'; smov(1).MarkerFaceAlpha = 0.5;
smov(2) = scatter(rand(size(sessionsv4_1))/5+1.9,cat(1,sessionsv4_1.movPerc),72,cat(1,sessionsv4_1.matchColor)/2,'filled');
smov(2).MarkerEdgeColor = 'none'; smov(2).MarkerFaceAlpha = 0.5;
smov(3) = scatter(rand(size(sessionsv4_n))/5+2.9,cat(1,sessionsv4_n.movPerc),72,cat(1,sessionsv4_n.matchColor)/2,'filled');
smov(3).MarkerEdgeColor = 'none'; smov(3).MarkerFaceAlpha = 0.5;
smov(4) = scatter(rand(size(sessionsv5))/5+3.9,cat(1,sessionsv5.movPerc),72,cat(1,sessionsv5.matchColor)/2,'filled');
smov(4).MarkerEdgeColor = 'none'; smov(4).MarkerFaceAlpha = 0.5;
ylim([0 1]); yticks([0:0.1:1]);
xlim([0.5 7.5]); xticks([1:4]); xticklabels({'Naive','Trained day 1','Trained day N','Adaptation'});
xtickangle(45);
ylabel('Proportion of trials with movement');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
axis square;

subplot(2,2,2);
hold on;
scorr(1) = scatter(rand(size(sessionsv4_1))/5+.9,cat(1,sessionsv4_1.corrPerc),72,cat(1,sessionsv4_1.matchColor)/2,'filled');
scorr(1).MarkerEdgeColor = 'none'; scorr(1).MarkerFaceAlpha = 0.5;
scorr(2) = scatter(rand(size(sessionsv4_n))/5+2.4,cat(1,sessionsv4_n.corrPerc),72,cat(1,sessionsv4_n.matchColor)/2,'filled');
scorr(2).MarkerEdgeColor = 'none'; scorr(2).MarkerFaceAlpha = 0.5;

scorr(3) = scatter(rand(size(sessionsv5))/5+3.9,cat(1,sessionsv5.corrPercBase),72,cat(1,sessionsv5.matchColor)/2,'filled');
scorr(3).MarkerEdgeColor = 'none'; scorr(3).MarkerFaceAlpha = 0.5;
scorr(4) = scatter(rand(size(sessionsv5))/5+4.9,cat(1,sessionsv5.corrPercEarly),72,cat(1,sessionsv5.matchColor)/2,'filled');
scorr(4).MarkerEdgeColor = 'none'; scorr(4).MarkerFaceAlpha = 0.5;
scorr(5) = scatter(rand(size(sessionsv5))/5+5.9,cat(1,sessionsv5.corrPercLate),72,cat(1,sessionsv5.matchColor)/2,'filled');
scorr(5).MarkerEdgeColor = 'none'; scorr(5).MarkerFaceAlpha = 0.5;
scorr(6) = scatter(rand(size(sessionsv5))/5+6.9,cat(1,sessionsv5.corrPercWash),72,cat(1,sessionsv5.matchColor)/2,'filled');
scorr(6).MarkerEdgeColor = 'none'; scorr(6).MarkerFaceAlpha = 0.5;
ylim([0 1]); yticks([0:0.1:1]);
xlim([0.5 7.5]); xticks([1,2.5,4:7]); xticklabels({'Trained day 1','Trained day N','Adapt. baseline',sprintf('Adapt. early (%s)',num2str(adaptTrialGrp)),'Adapt. late','Adapt. washout'});
xtickangle(45);
ylabel('Proportion of trials correct');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
axis square;

subplot(2,2,3);
hold on;
scorr(1) = scatter(rand(size(sessionsv4_1))/5+.9,cat(1,sessionsv4_1.xEndMu),72,cat(1,sessionsv4_1.matchColor)/2,'filled');
scorr(1).MarkerEdgeColor = 'none'; scorr(1).MarkerFaceAlpha = 0.5;
scorr(2) = scatter(rand(size(sessionsv4_n))/5+2.4,cat(1,sessionsv4_n.xEndMu),72,cat(1,sessionsv4_n.matchColor)/2,'filled');
scorr(2).MarkerEdgeColor = 'none'; scorr(2).MarkerFaceAlpha = 0.5;
scorr(3) = scatter(rand(size(sessionsv5))/5+3.9,cat(1,sessionsv5.xEndMu_Base),72,cat(1,sessionsv5.matchColor)/2,'filled');
scorr(3).MarkerEdgeColor = 'none'; scorr(3).MarkerFaceAlpha = 0.5;
scorr(4) = scatter(rand(size(sessionsv5))/5+4.9,cat(1,sessionsv5.xEndMu_Early),72,cat(1,sessionsv5.matchColor)/2,'filled');
scorr(4).MarkerEdgeColor = 'none'; scorr(4).MarkerFaceAlpha = 0.5;
scorr(5) = scatter(rand(size(sessionsv5))/5+5.9,cat(1,sessionsv5.xEndMu_Late),72,cat(1,sessionsv5.matchColor)/2,'filled');
scorr(5).MarkerEdgeColor = 'none'; scorr(5).MarkerFaceAlpha = 0.5;
scorr(6) = scatter(rand(size(sessionsv5))/5+6.9,cat(1,sessionsv5.xEndMu_Wash),72,cat(1,sessionsv5.matchColor)/2,'filled');
scorr(6).MarkerEdgeColor = 'none'; scorr(6).MarkerFaceAlpha = 0.5;
ylim([-1 1]); yticks([-1:0.25:1]);
xlim([0.5 7.5]); xticks([1,2.5,4:7]); xticklabels({'Trained day 1','Trained day N','Adapt. baseline',sprintf('Adapt. early (%s)',num2str(adaptTrialGrp)),'Adapt. late','Adapt. washout'});
xtickangle(45);
ylabel('Mean end position');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
axis square;


% subplot(1,3,3);
% hold on;
% scorr(1) = scatter(rand(size(sessionsv4_1))/5+.9,cat(1,sessionsv4_1.xEndMed),72,cat(1,sessionsv4_1.matchColor)/2,'filled');
% scorr(1).MarkerEdgeColor = 'none'; scorr(1).MarkerFaceAlpha = 0.5;
% scorr(2) = scatter(rand(size(sessionsv4_n))/5+2.4,cat(1,sessionsv4_n.xEndMed),72,cat(1,sessionsv4_n.matchColor)/2,'filled');
% scorr(2).MarkerEdgeColor = 'none'; scorr(2).MarkerFaceAlpha = 0.5;
% scorr(3) = scatter(rand(size(sessionsv5))/5+3.9,cat(1,sessionsv5.xEndMed_Base),72,cat(1,sessionsv5.matchColor)/2,'filled');
% scorr(3).MarkerEdgeColor = 'none'; scorr(3).MarkerFaceAlpha = 0.5;
% scorr(4) = scatter(rand(size(sessionsv5))/5+4.9,cat(1,sessionsv5.xEndMed_Early),72,cat(1,sessionsv5.matchColor)/2,'filled');
% scorr(4).MarkerEdgeColor = 'none'; scorr(4).MarkerFaceAlpha = 0.5;
% scorr(5) = scatter(rand(size(sessionsv5))/5+5.9,cat(1,sessionsv5.xEndMed_Late),72,cat(1,sessionsv5.matchColor)/2,'filled');
% scorr(5).MarkerEdgeColor = 'none'; scorr(5).MarkerFaceAlpha = 0.5;
% scorr(6) = scatter(rand(size(sessionsv5))/5+6.9,cat(1,sessionsv5.xEndMed_Wash),72,cat(1,sessionsv5.matchColor)/2,'filled');
% scorr(6).MarkerEdgeColor = 'none'; scorr(6).MarkerFaceAlpha = 0.5;
% ylim([-1 1]); yticks([-1:0.25:1]);
% xlim([0.5 7.5]); xticks([1,2.5,4:7]); xticklabels({'Trained day 1','Trained day N','Adapt. baseline',sprintf('Adapt. early (%s)',num2str(adaptTrialGrp)),'Adapt. late','Adapt. washout'});
% xtickangle(45);
% ylabel('Median end position');
% set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
% axis square;


% First figure: fraction of trials with movement
figout(2) = figure;
subplot(2,2,1); hold on;
bar(1,mean(cat(1,seshV1.movPerc)),0.6,'EdgeColor','none','FaceColor',[.6 0 0],'FaceAlpha',0.5);
errorbar(1, mean(cat(1,seshV1.movPerc)), std(cat(1,seshV1.movPerc))/sqrt(numel(seshV1)),'k','CapSize',0,'LineWidth',2);
bar(2,mean(cat(1,sessionsv4_1.movPerc)),0.6,'EdgeColor','none','FaceColor',[.6 .6 .6],'FaceAlpha',0.5);
errorbar(2, mean(cat(1,sessionsv4_1.movPerc)), std(cat(1,sessionsv4_1.movPerc))/sqrt(numel(sessionsv4_1)),'k','CapSize',0,'LineWidth',2);
bar(3,mean(cat(1,sessionsv4_n.movPerc)),0.6,'EdgeColor','none','FaceColor',[.3 .3 .3],'FaceAlpha',0.5);
errorbar(3, mean(cat(1,sessionsv4_n.movPerc)), std(cat(1,sessionsv4_n.movPerc))/sqrt(numel(sessionsv4_n)),'k','CapSize',0,'LineWidth',2);
bar(4,mean(cat(1,sessionsv5.movPerc)),0.6,'EdgeColor','none','FaceColor',[.0 .2 1/3],'FaceAlpha',0.5);
errorbar(4, mean(cat(1,sessionsv5.movPerc)), std(cat(1,sessionsv5.movPerc))/sqrt(numel(sessionsv5)),'k','CapSize',0,'LineWidth',2);

ylim([0 1]); yticks([0:0.2:1]);
xlim([0.5 7.5]); xticks([1:4]); xticklabels({'Naive','Trained day 1','Trained day N','Adaptation'});
xtickangle(45);
ylabel('Proportion of trials with movement');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
axis square;

subplot(2,2,3);
hold on;
bar(1,mean(cat(1,sessionsv4_1.corrPerc)),0.6,'EdgeColor','none','FaceColor',[.6 .6 .6],'FaceAlpha',0.5);
errorbar(1, mean(cat(1,sessionsv4_1.corrPerc)), std(cat(1,sessionsv4_1.corrPerc))/sqrt(numel(sessionsv4_1)),'k','CapSize',0,'LineWidth',2);
bar(2.5,mean(cat(1,sessionsv4_n.corrPerc)),0.6,'EdgeColor','none','FaceColor',[.3 .3 .3],'FaceAlpha',0.5);
errorbar(2.5, mean(cat(1,sessionsv4_n.corrPerc)), std(cat(1,sessionsv4_n.corrPerc))/sqrt(numel(sessionsv4_n)),'k','CapSize',0,'LineWidth',2);
bar(4,mean(cat(1,sessionsv5.corrPercBase)),0.6,'EdgeColor','none','FaceColor',[.0 .2 1/3],'FaceAlpha',0.5);
errorbar(4, mean(cat(1,sessionsv5.corrPercBase)), std(cat(1,sessionsv5.corrPercBase))/sqrt(numel(sessionsv5)),'k','CapSize',0,'LineWidth',2);
bar(5,mean(cat(1,sessionsv5.corrPercEarly)),0.6,'EdgeColor','none','FaceColor',[.0 .4 2/3],'FaceAlpha',0.5);
errorbar(5, mean(cat(1,sessionsv5.corrPercEarly)), std(cat(1,sessionsv5.corrPercEarly))/sqrt(numel(sessionsv5)),'k','CapSize',0,'LineWidth',2);
bar(6,mean(cat(1,sessionsv5.corrPercLate)),0.6,'EdgeColor','none','FaceColor',[.0 .6 1],'FaceAlpha',0.5);
errorbar(6, mean(cat(1,sessionsv5.corrPercLate)), std(cat(1,sessionsv5.corrPercLate))/sqrt(numel(sessionsv5)),'k','CapSize',0,'LineWidth',2);
bar(7,nanmean(cat(1,sessionsv5.corrPercWash)),0.6,'EdgeColor','none','FaceColor',[.0 .9 .6],'FaceAlpha',0.5);
errorbar(7, nanmean(cat(1,sessionsv5.corrPercWash)), nanstd(cat(1,sessionsv5.corrPercWash))/sqrt(numel(cat(1,sessionsv5.corrPercWash))),'k','CapSize',0,'LineWidth',2);
ylim([0 1]); yticks([0:0.2:1]);
xlim([0.5 7.5]); xticks([1,2.5,4:7]); xticklabels({'Trained day 1','Trained day N','Adapt. baseline',sprintf('Adapt. early (%s)',num2str(adaptTrialGrp)),'Adapt. late','Adapt. washout'});
xtickangle(45);
ylabel('Proportion of trials correct');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
axis square;

subplot(2,2,2);
hold on;
bar(1,mean(cat(1,sessionsv4_1.xEndMu)),0.6,'EdgeColor','none','FaceColor',[.6 .6 .6],'FaceAlpha',0.5);
errorbar(1, mean(cat(1,sessionsv4_1.xEndMu)), std(cat(1,sessionsv4_1.xEndMu))/sqrt(numel(sessionsv4_1)),'k','CapSize',0,'LineWidth',2);
bar(2.5,mean(cat(1,sessionsv4_n.xEndMu)),0.6,'EdgeColor','none','FaceColor',[.3 .3 .3],'FaceAlpha',0.5);
errorbar(2.5, mean(cat(1,sessionsv4_n.xEndMu)), std(cat(1,sessionsv4_n.xEndMu))/sqrt(numel(sessionsv4_n)),'k','CapSize',0,'LineWidth',2);
bar(4,mean(cat(1,sessionsv5.xEndMu_Base)),0.6,'EdgeColor','none','FaceColor',[.0 .2 1/3],'FaceAlpha',0.5);
errorbar(4, mean(cat(1,sessionsv5.xEndMu_Base)), std(cat(1,sessionsv5.xEndMu_Base))/sqrt(numel(sessionsv5)),'k','CapSize',0,'LineWidth',2);
bar(5,mean(cat(1,sessionsv5.xEndMu_Early)),0.6,'EdgeColor','none','FaceColor',[.0 .4 2/3],'FaceAlpha',0.5);
errorbar(5, mean(cat(1,sessionsv5.xEndMu_Early)), std(cat(1,sessionsv5.xEndMu_Early))/sqrt(numel(sessionsv5)),'k','CapSize',0,'LineWidth',2);
bar(6,mean(cat(1,sessionsv5.xEndMu_Late)),0.6,'EdgeColor','none','FaceColor',[.0 .6 1],'FaceAlpha',0.5);
errorbar(6, mean(cat(1,sessionsv5.xEndMu_Late)), std(cat(1,sessionsv5.xEndMu_Late))/sqrt(numel(sessionsv5)),'k','CapSize',0,'LineWidth',2);
bar(7,nanmean(cat(1,sessionsv5.xEndMu_Wash)),0.6,'EdgeColor','none','FaceColor',[.0 .9 .6],'FaceAlpha',0.5);
errorbar(7, nanmean(cat(1,sessionsv5.xEndMu_Wash)), nanstd(cat(1,sessionsv5.xEndMu_Wash))/sqrt(numel(cat(1,sessionsv5.xEndMu_Wash))),'k','CapSize',0,'LineWidth',2);
ylim([-1 1]); yticks([-1:0.25:1]);
xlim([0.5 7.5]); xticks([1,2.5,4:7]); xticklabels({}); %xticklabels({'Trained day 1','Trained day N','Adapt. baseline',sprintf('Adapt. early (%s)',num2str(adaptTrialGrp)),'Adapt. late','Adapt. washout'});
xtickangle(45);
ylabel('Mean displacement');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
axis square;

subplot(2,2,4);
hold on;
rectangle('Position',[0 2/3 8 2/3],'FaceColor',[.9 .9 .9 .5],'EdgeColor','none')
yline(1,'LineWidth',1);
bar(1,1-mean(cat(1,sessionsv4_1.xEndMu)),0.6,'EdgeColor','none','FaceColor',[.6 .6 .6],'FaceAlpha',0.5);
errorbar(1, 1-mean(cat(1,sessionsv4_1.xEndMu)), std(cat(1,sessionsv4_1.xEndMu))/sqrt(numel(sessionsv4_1)),'k','CapSize',0,'LineWidth',2);
bar(2.5,1-mean(cat(1,sessionsv4_n.xEndMu)),0.6,'EdgeColor','none','FaceColor',[.3 .3 .3],'FaceAlpha',0.5);
errorbar(2.5, 1-mean(cat(1,sessionsv4_n.xEndMu)), std(cat(1,sessionsv4_n.xEndMu))/sqrt(numel(sessionsv4_n)),'k','CapSize',0,'LineWidth',2);
bar(4,1-mean(cat(1,sessionsv5.xEndMu_Base)),0.6,'EdgeColor','none','FaceColor',[.0 .2 1/3],'FaceAlpha',0.5);
errorbar(4, 1-mean(cat(1,sessionsv5.xEndMu_Base)), std(cat(1,sessionsv5.xEndMu_Base))/sqrt(numel(sessionsv5)),'k','CapSize',0,'LineWidth',2);
bar(5,1.6*(1-mean(cat(1,sessionsv5.xEndMu_Early))),0.6,'EdgeColor','none','FaceColor',[.0 .4 2/3],'FaceAlpha',0.5);
errorbar(5, 1.6*(1-mean(cat(1,sessionsv5.xEndMu_Early))), 1.6*std(cat(1,sessionsv5.xEndMu_Early))/sqrt(numel(sessionsv5)),'k','CapSize',0,'LineWidth',2);
bar(6,1.6*(1-mean(cat(1,sessionsv5.xEndMu_Late))),0.6,'EdgeColor','none','FaceColor',[.0 .6 1],'FaceAlpha',0.5);
errorbar(6, 1.6*(1-mean(cat(1,sessionsv5.xEndMu_Late))), 1.6*std(cat(1,sessionsv5.xEndMu_Late))/sqrt(numel(sessionsv5)),'k','CapSize',0,'LineWidth',2);
bar(7,1-nanmean(cat(1,sessionsv5.xEndMu_Wash)),0.6,'EdgeColor','none','FaceColor',[.0 .9 .6],'FaceAlpha',0.5);
errorbar(7, 1-nanmean(cat(1,sessionsv5.xEndMu_Wash)), nanstd(cat(1,sessionsv5.xEndMu_Wash))/sqrt(numel(cat(1,sessionsv5.xEndMu_Wash))),'k','CapSize',0,'LineWidth',2);
ylim([0 5/3]); yticks([0:0.25:1.5]);
xlim([0.5 7.5]); xticks([1,2.5,4:7]); xticklabels({'Trained day 1','Trained day N','Adapt. baseline',sprintf('Adapt. early (%s)',num2str(adaptTrialGrp)),sprintf('Adapt. late (%s)',num2str(adaptTrialGrp)),'Adapt. washout'});
xtickangle(45);
ylabel('Mean end position');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
axis square;

if savebool
    if ~exist(figFold,'dir'); mkdir(figFold); end
    figname = sprintf('percMovCorrect_%i_trials',adaptTrialGrp);
    savefig(figout(1), fullfile(figFold,[figname,'.fig']));
    saveas(figout(1),fullfile(figFold,[figname,'.png']));
    print(figout(1),fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    
    if ~exist(figFold2,'dir'); mkdir(figFold2); end
    savefig(figout(1), fullfile(figFold2,[figname,'.fig']));
    saveas(figout(1),fullfile(figFold2,[figname,'.png']));
    print(figout(1),fullfile(figFold2,[figname,'.eps']), '-depsc', '-painters');
    
    figname = sprintf('percMovCorrect_%i_trials_bar',adaptTrialGrp);
    savefig(figout(2), fullfile(figFold,[figname,'.fig']));
    saveas(figout(2),fullfile(figFold,[figname,'.png']));
    print(figout(2),fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    
    if ~exist(figFold2,'dir'); mkdir(figFold2); end
    savefig(figout(2), fullfile(figFold2,[figname,'.fig']));
    saveas(figout(2),fullfile(figFold2,[figname,'.png']));
    print(figout(2),fullfile(figFold2,[figname,'.eps']), '-depsc', '-painters');


end


%% Stats!

% corrv4_1 = cat(1,sessionsv4_1.corrPerc); corrv4_1 = [corrv4_1 ones(size(corrv4_1))];
% corrv4_n = cat(1,sessionsv4_n.corrPerc); corrv4_n = [corrv4_n 2*ones(size(corrv4_n))];
% corrv5base = cat(1,sessionsv5.corrPercBase); corrv5base = [corrv5base 3*ones(size(corrv5base))];
% corrv5early = cat(1,sessionsv5.corrPercEarly); corrv5early = [corrv5early 4*ones(size(corrv5early))];
% corrv5late = cat(1,sessionsv5.corrPercLate); corrv5late = [corrv5late 5*ones(size(corrv5late))];
% corrv5wash = cat(1,sessionsv5.corrPercWash); corrv5wash = [corrv5wash 6*ones(size(corrv5wash))];
% 
% % [corrP, corrTab, corrStats] = anova1([corrv4_1(:,1); corrv4_n(:,1); corrv5base(:,1)],[corrv4_1(:,2); corrv4_n(:,2); corrv5base(:,2)]);
% % corrC = multcompare(corrStats, 'Ctype','bonferroni');
% 
% [corrP, corrTab, corrStats] = anova1([corrv5base(:,1); corrv5early(:,1); corrv5late(:,1); corrv5wash(:,1)],[corrv5base(:,2); corrv5early(:,2); corrv5late(:,2); corrv5wash(:,2)]);
% corrC = multcompare(corrStats, 'Ctype','bonferroni');
% 
% 
% movv4_1 = cat(1,sessionsv4_1.xEndMu); movv4_1 = [movv4_1 ones(size(movv4_1))];
% movv4_n = cat(1,sessionsv4_n.xEndMu); movv4_n = [movv4_n 2*ones(size(movv4_n))];
% movv5base = cat(1,sessionsv5.xEndMu_Base); movv5base = [movv5base 3*ones(size(movv5base))];
% movv5early = cat(1,sessionsv5.xEndMu_Early); movv5early = [movv5early 4*ones(size(movv5early))];
% movv5late = cat(1,sessionsv5.xEndMu_Late); movv5late = [movv5late 5*ones(size(movv5late))];
% movv5wash = cat(1,sessionsv5.xEndMu_Wash); movv5wash = [movv5wash 6*ones(size(movv5wash))];
% 
% % [movP, movTab, movStats] = anova1([movv4_1(:,1); movv4_n(:,1); movv5base(:,1)],[movv4_1(:,2); movv4_n(:,2); movv5base(:,2)]);
% % movC = multcompare(movStats, 'Ctype','bonferroni');
% 
% [movP, movTab, movStats] = anova1([movv5base(:,1); movv5early(:,1); movv5late(:,1); movv5wash(:,1)],[movv5base(:,2); movv5early(:,2); movv5late(:,2); movv5wash(:,2)]);
% movC = multcompare(movStats, 'Ctype','bonferroni');
% movC2 = multcompare(movStats); %, 'Ctype','bonferroni');
% 
% 
% offv4_1 = 1-cat(1,sessionsv4_1.xEndMu); offv4_1 = [offv4_1 ones(size(offv4_1))];
% offv4_n = 1-cat(1,sessionsv4_n.xEndMu); offv4_n = [offv4_n 2*ones(size(offv4_n))];
% offv5base = 1-cat(1,sessionsv5.xEndMu_Base); offv5base = [offv5base 3*ones(size(offv5base))];
% offv5early = 1.6*(1-cat(1,sessionsv5.xEndMu_Early)); offv5early = [offv5early 4*ones(size(offv5early))];
% offv5late = 1.6*(1-cat(1,sessionsv5.xEndMu_Late)); offv5late = [offv5late 5*ones(size(offv5late))];
% offv5wash = 1-cat(1,sessionsv5.xEndMu_Wash); offv5wash = [offv5wash 6*ones(size(offv5wash))];
% 
% % [offP, offTab, offStats] = anova1([offv4_1(:,1); offv4_n(:,1); offv5base(:,1)],[offv4_1(:,2); offv4_n(:,2); offv5base(:,2)]);
% % offC = multcompare(offStats, 'Ctype','bonferroni');
% 
% [offP, offTab, offStats] = anova1([offv5base(:,1); offv5early(:,1); offv5late(:,1); offv5wash(:,1)],[offv5base(:,2); offv5early(:,2); offv5late(:,2); offv5wash(:,2)]);
% offC = multcompare(offStats, 'Ctype','bonferroni');
% offC2 = multcompare(offStats); %, 'Ctype','bonferroni');

if closebool
    close all
end
end