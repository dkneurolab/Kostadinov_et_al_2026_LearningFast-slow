function fig1_percMovCorrect(sessions,dataFold,paperFold,adaptTrialGrp,matchColors,closebool,savebool)
%% Concatenate sessions and split by version and set basic parameters
% Split sessions
sessionsv1 = sessions(cat(1,sessions.version) == 1);
sessionsv4 = sessions(cat(1,sessions.version) == 4);
sessionsv5 = sessions(cat(1,sessions.version) == 5);

% Set directories for behavioural data and output figures
behavFold = fullfile(dataFold,'04_Behav-processed');
figFold = fullfile(paperFold,'Figures','Fig1','Sub1_percMovCorr');

% Set behav parameters
movThreshold = 1.5; % Threshold for called a wheel turn a real movement (cm/s)

%% Version 1 section
if ~exist(fullfile(figFold,'percMovCorrect_v1.mat'),'file')
    % If this has not already been done, iterate over session and extract
    % percent correct and fraction of trials with movement

    for i = 1:numel(sessionsv1)
        % Load behavioural session
        dataDir = dir(fullfile(behavFold,sessionsv1(i).name,sessionsv1(i).date,'Behav*',[sessionsv1(i).name,'_',sessionsv1(i).date,'.mat']));
        datastruct = load(fullfile(dataDir.folder,dataDir.name)); datastruct = datastruct.datastruct;
        trials = datastruct.paqdata.behav.trials;
        inits = datastruct.inits;
        sets = datastruct.sets;
        gaindiv = sets.start1/360*4.5*pi/sets.rotgain; % Convert vr gain into cm (Lego wheel has D = 4.5 cm)
        if sets.vrfs == 60 % For sessions when VR ran at 60 hz
            smspan = inits.paqfs*.04; % 40 ms smoothing span
        elseif sets. vrfs == 30 % For sessions when VR ran at 30 hz
            smspan = inits.paqfs*.06; % 60 ms smoothing span
        end

        % Extract basic features from each trial
        sessionsv1(i).nTrials = numel(trials);
        sessionsv1(i).iTrialsReal = false(sessionsv1(i).nTrials,1);
        sessionsv1(i).movBool = false(sessionsv1(i).nTrials,1);
        sessionsv1(i).corrBool = false(sessionsv1(i).nTrials,1);
        % Determine if movement happened and if it resulted in correct trial
        for j = 1:numel(trials)
            x0 = trials(j).data(inits.paqfs*inits.mintime+1:end-inits.paqfs*inits.mintime,10);
            x1 = smooth(x0,smspan/numel(x0),'lowess');
            v1 = diff(x1)*inits.paqfs*-gaindiv; v1(end+1) = v1(end);
            if numel(v1) > 15*inits.paqfs; v1 = v1(1:15*inits.paqfs); end
            if numel(v1)/inits.paqfs <= sets.responsetime
                sessionsv1(i).iTrialsReal(j) = true;
            end
            % If wheel reached threshold cm/s on a trail --> movement happened
            if any(abs(v1) > movThreshold) % movement threshold is 1.5 cm/s mm of wheel translation
                sessionsv1(i).movBool(j) = true;
            end
            if sessionsv1(i).iTrialsReal(j) && trials(j).correct
                sessionsv1(i).corrBool(j) = true;
            end
        end
        % Calculate percent of trials with movement
        sessionsv1(i).movPerc = sum(sessionsv1(i).movBool)/sessionsv1(i).nTrials;
        sessionsv1(i).nTrialsReal = sum(sessionsv1(i).iTrialsReal);
        if sessionsv1(i).nTrialsReal > 0
            sessionsv1(i).corrPerc = sum(sessionsv1(i).corrBool)/sessionsv1(i).nTrialsReal;
        else
            sessionsv1(i).corrPerc = 0;
        end
        % Keep local settings file
        sessionsv1(i).sets = sets;
        
        switch lower(sessionsv1(i).name)
            case 'dk103'
                sessionsv1(i).matchColor = matchColors(1,:);
            case 'dk105'
                sessionsv1(i).matchColor = matchColors(2,:);
            case 'dk169'
                sessionsv1(i).matchColor = matchColors(3,:);
            case 'dk171'
                sessionsv1(i).matchColor = matchColors(4,:);
            case 'dk194'
                sessionsv1(i).matchColor = matchColors(5,:);
            case 'dk199'
                sessionsv1(i).matchColor = matchColors(6,:);
            case 'dk052'
                sessionsv1(i).matchColor = matchColors(7,:);
            case 'dk063'
                sessionsv1(i).matchColor = matchColors(8,:);
            case 'dk070'
                sessionsv1(i).matchColor = matchColors(9,:);
            case 'dk138'
                sessionsv1(i).matchColor = matchColors(10,:);
        end

    end

    if savebool
        if ~exist(figFold,'dir'); mkdir(figFold); end
        save(fullfile(figFold,'percMovCorrect_v1.mat'),'sessionsv1');
        if ~exist(figFold2,'dir'); mkdir(figFold2); end
        save(fullfile(figFold2,'percMovCorrect_v1.mat'),'sessionsv1');
    end

else
    % If this has already been done, just load the file to save time
    sessionsv1 = load(fullfile(figFold,'percMovCorrect_v1.mat'));
    sessionsv1 = sessionsv1.sessionsv1;
end

%% Version 4 section - same as above but for v4 (learning / expert task version)
if exist(fullfile(figFold,'percMovCorrect_v4.mat'),'file')
    % If this has not already been done, iterate over session and extract
    % percent correct and fraction of trials with movement
    for i = 1:numel(sessionsv4)
        dataDir = dir(fullfile(behavFold,sessionsv4(i).name,sessionsv4(i).date,'Behav*',[sessionsv4(i).name,'_',sessionsv4(i).date,'.mat']));
        datastruct = load(fullfile(dataDir.folder,dataDir.name)); datastruct = datastruct.datastruct;
        trials = datastruct.paqdata.behav.trials;
        inits = datastruct.inits;
        sets = datastruct.sets;
        gaindiv = 37.5/360*4.5*pi/sets.rotgain;
        if sets.vrfs == 60
            smspan = inits.paqfs*.04; % 40 ms smoothing span
        elseif sets. vrfs == 30
            smspan = inits.paqfs*.06; % 60 ms smoothing span
        end
        if strcmpi(sessionsv4(i).name,'dk194') & strcmp(sessionsv4(i).date,'2020-10-02')
            trials = trials(1:180); % Account for lack of motivation at end of this session
        end

        % Remove non-responsive trials at end
        trials = trials(1:find(cat(1,trials.realtrial) > 0,1, 'last'));

        behavstability = Rot2p_behavstability(trials,inits,sets,0);

        sessionsv4(i).nTrials = numel(trials);
        sessionsv4(i).nTrialsReal = sum(cat(1,trials.realtrial));
        sessionsv4(i).iTrialsReal = false(sessionsv4(i).nTrials,1);
        sessionsv4(i).movBool = false(sessionsv4(i).nTrials,1);
        sessionsv4(i).corrBool = false(sessionsv4(i).nTrials,1);
        % Determine if movement happened and if it resulted in correct trial
        for j = 1:numel(trials)
            x0 = trials(j).data(inits.paqfs*inits.mintime+1:end-inits.paqfs*inits.mintime,10);
            x1 = smooth(x0,smspan/numel(x0),'lowess');
            v1 = diff(x1)*inits.paqfs*-gaindiv; v1(end+1) = v1(end);
            if numel(v1)/inits.paqfs <= sets.responsetime
                sessionsv4(i).iTrialsReal(j) = true;
            end
            if numel(v1) > 15*inits.paqfs; v1 = v1(1:15*inits.paqfs); end
            if any(abs(v1) > 1.5) % movement threshold is 1.5 cm/s of wheel translation
                %         if any(posTrace > sets.rotgain*360/45/pi()) % movement threshold is 1 mm of wheel translation
                sessionsv4(i).movBool(j) = true;
            end
            if sessionsv4(i).iTrialsReal(j) & trials(j).correct
                sessionsv4(i).corrBool(j) = true;
            end
        end
        sessionsv4(i).movPerc = sum(sessionsv4(i).movBool)/sessionsv4(i).nTrials;
        %     sessionsv4(i).nTrialsReal = sum(sessionsv4(i).iTrialsReal);
        if sessionsv4(i).nTrialsReal > 0
            sessionsv4(i).corrPerc = sum(sessionsv4(i).corrBool)/sessionsv4(i).nTrialsReal;
            sessionsv4(i).xEndMu = mean(behavstability.params.xendall);
            sessionsv4(i).xEndMed = median(behavstability.params.xendall);
        else
            sessionsv4(i).corrPerc = 0;
        end

        sessionsv4(i).sets = sets;
        sessionsv4(i).params = rmfield(behavstability.params,{'trialnum','blocksize','rhoU','rhoC','rhoUmu','rhoUsd','rhoCmu','rhoCsd','rhoOmu','rhoOsd','numC'});
        sessionsv4(i).offXall = behavstability.offxall;
        sessionsv4(i).movVall = behavstability.movvall;

        % Alternative version - first session should always be day 1, others
        % should not mach if they are day 1
        if i > 1
            if strcmp(sessionsv4(i).name,sessionsv4(i-1).name) & strcmp(sessionsv4(i).fov,sessionsv4(i-1).fov)
                sessionsv4(i).day1 = false;
                sessionsv4(i).dayN = true;
            else
                sessionsv4(i).day1 = true;
                sessionsv4(i).dayN = false;
            end
        else
            sessionsv4(i).day1 = true;
            sessionsv4(i).dayN = false;
        end


        switch lower(sessionsv4(i).name)
            case 'dk103'
                sessionsv4(i).matchColor = matchColors(1,:);
            case 'dk105'
                sessionsv4(i).matchColor = matchColors(2,:);
            case 'dk169'
                sessionsv4(i).matchColor = matchColors(3,:);
            case 'dk171'
                sessionsv4(i).matchColor = matchColors(4,:);
            case 'dk194'
                sessionsv4(i).matchColor = matchColors(5,:);
            case 'dk199'
                sessionsv4(i).matchColor = matchColors(6,:);
            case 'dk052'
                sessionsv4(i).matchColor = matchColors(7,:);
            case 'dk063'
                sessionsv4(i).matchColor = matchColors(8,:);
            case 'dk070'
                sessionsv4(i).matchColor = matchColors(9,:);
            case 'dk138'
                sessionsv4(i).matchColor = matchColors(10,:);
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
    sessionsv4 = load(fullfile(figFold,'percMovCorrect_v4.mat'));
    sessionsv4 = sessionsv4.sessionsv4;
end

%% Version 5 section
if ~exist(fullfile(figFold,'percMovCorrect_v5.mat'),'file')
    % if ~exist(fullfile(figFold,['percMovCorrect_v5_',num2str(adaptTrialGrp),'trials.mat']),'file')
    % If this has not already been done, iterate over session and extract
    % percent correct and fraction of trials with movement
    for i = 1:numel(sessionsv5)
        dataDir = dir(fullfile(behavFold,sessionsv5(i).name,sessionsv5(i).date,'Behav*',[sessionsv5(i).name,'_',sessionsv5(i).date,'.mat']));
        datastruct = load(fullfile(dataDir.folder,dataDir.name)); datastruct = datastruct.datastruct;
        trials = datastruct.paqdata.behav.trials;
        inits = datastruct.inits;
        sets = Rot2p_changetrials(datastruct.sets);
        gaindiv = 37.5/360*4.5*pi/sets.rotgain;
        if sets.vrfs == 60
            smspan = inits.paqfs*.04; % 40 ms smoothing span
        elseif sets. vrfs == 30
            smspan = inits.paqfs*.06; % 60 ms smoothing span
        end
        if strcmpi(sessionsv5(i).name,'dk105') & strcmp(sessionsv5(i).date,'2018-03-07')
            sets.changebacktrial = sets.changetrial+10*(sets.changebackmult-1);
        end
        trials = trials(1:find(cat(1,trials.realtrial) > 0,1, 'last'));

        behavstability = Rot2p_behavstability(trials,inits,sets,0);

        sessionsv5(i).nTrials = numel(trials);
        sessionsv5(i).nTrialsReal = sum(cat(1,trials.realtrial));
        sessionsv5(i).nTrialsBase = sum(cat(1,trials(1:sets.changetrial).realtrial));
        sessionsv5(i).nTrialsEarly = sum(cat(1,trials(sets.changetrial+1:sets.changetrial+adaptTrialGrp).realtrial));
        sessionsv5(i).nTrialsLate = sum(cat(1,trials(sets.changetrial+121-adaptTrialGrp:sets.changetrial+120).realtrial));

        if ~sets.changeback
            sessionsv5(i).nTrialsWash = 0;
        else
            nWashTrials = numel(trials)-sets.changebacktrial;
            if nWashTrials > adaptTrialGrp; nWashTrials = adaptTrialGrp; end
            sessionsv5(i).nTrialsWash = sum(cat(1,trials(sets.changebacktrial+1:sets.changebacktrial+nWashTrials).realtrial));
        end

        sessionsv5(i).iTrialsReal = false(sessionsv5(i).nTrials,1);
        sessionsv5(i).movBool = false(sessionsv5(i).nTrials,1);
        sessionsv5(i).corrBool = false(sessionsv5(i).nTrials,1);

        % Determine if movement happened and if it resulted in correct trial
        for j = 1:numel(trials)
            x0 = trials(j).data(inits.paqfs*inits.mintime+1:end-inits.paqfs*inits.mintime,10);
            x1 = smooth(x0,smspan/numel(x0),'lowess');
            v1 = diff(x1)*inits.paqfs*-gaindiv; v1(end+1) = v1(end);
            if numel(v1) > 15*inits.paqfs; v1 = v1(1:15*inits.paqfs); end
            if numel(v1)/inits.paqfs <= sets.responsetime
                sessionsv5(i).iTrialsReal(j) = true;
            end
            if any(abs(v1) > 1.5) % movement threshold is 1.5 cm/s of wheel translation
                %         if any(posTrace > sets.rotgain*360/45/pi()) % movement threshold is 1 mm of wheel translation
                sessionsv5(i).movBool(j) = true;
            end
            if sessionsv5(i).iTrialsReal(j) & trials(j).correct
                sessionsv5(i).corrBool(j) = true;
            end
        end
        sessionsv5(i).movPerc = sum(sessionsv5(i).movBool)/sessionsv5(i).nTrials;
        if sessionsv5(i).nTrialsReal > 0
            sessionsv5(i).corrPerc = sum(sessionsv5(i).corrBool)/sessionsv5(i).nTrialsReal;
            sessionsv5(i).corrPercBase = sum(sessionsv5(i).corrBool(1:sets.changetrial))/sessionsv5(i).nTrialsBase;
            sessionsv5(i).corrPercEarly = sum(sessionsv5(i).corrBool(sets.changetrial+1:sets.changetrial+adaptTrialGrp))/sessionsv5(i).nTrialsEarly;
            sessionsv5(i).corrPercLate = sum(sessionsv5(i).corrBool(sets.changetrial+121-adaptTrialGrp:sets.changetrial+120))/sessionsv5(i).nTrialsLate;
            if ~sets.changeback
                sessionsv5(i).corrPercWash = nan;
            else
                sessionsv5(i).corrPercWash = sum(sessionsv5(i).corrBool(sets.changebacktrial+1:sets.changebacktrial+nWashTrials))/sessionsv5(i).nTrialsWash;
            end
            xEndAll = nan(size(sessionsv5(i).iTrialsReal)); xEndAll(behavstability.params.trialinds) = behavstability.params.xendall;

            sessionsv5(i).xEndMu_Base = nanmean(xEndAll(1:sets.changetrial));
            sessionsv5(i).xEndMed_Base = nanmedian(xEndAll(1:sets.changetrial));
            sessionsv5(i).xEndMu_Early = nanmean(xEndAll(sets.changetrial+1:sets.changetrial+adaptTrialGrp));
            sessionsv5(i).xEndMed_Early = nanmedian(xEndAll(sets.changetrial+1:sets.changetrial+adaptTrialGrp));
            sessionsv5(i).xEndMu_Late = nanmean(xEndAll(sets.changetrial+121-adaptTrialGrp:sets.changetrial+120));
            sessionsv5(i).xEndMed_Late = nanmedian(xEndAll(sets.changetrial+121-adaptTrialGrp:sets.changetrial+120));
            if ~sets.changeback
                sessionsv5(i).xEndMu_Wash = nan;
                sessionsv5(i).xEndMed_Wash = nan;
            else
                sessionsv5(i).xEndMu_Wash = nanmean(xEndAll(sets.changebacktrial+1:sets.changebacktrial+nWashTrials));
                sessionsv5(i).xEndMed_Wash = nanmedian(xEndAll(sets.changebacktrial+1:sets.changebacktrial+nWashTrials));
            end
        else
            sessionsv5(i).corrPerc = 0;
        end
        sessionsv5(i).sets = sets;
        sessionsv5(i).params = rmfield(behavstability.params,{'trialnum','blocksize','rhoU','rhoC','rhoUmu','rhoUsd','rhoCmu','rhoCsd','rhoOmu','rhoOsd','numC'});
        sessionsv5(i).offXall = behavstability.offxall;
        sessionsv5(i).movVall = behavstability.movvall;

        switch lower(sessionsv5(i).name)
            case 'dk103'
                sessionsv5(i).matchColor = matchColors(1,:);
            case 'dk105'
                sessionsv5(i).matchColor = matchColors(2,:);
            case 'dk169'
                sessionsv5(i).matchColor = matchColors(3,:);
            case 'dk171'
                sessionsv5(i).matchColor = matchColors(4,:);
            case 'dk194'
                sessionsv5(i).matchColor = matchColors(5,:);
            case 'dk199'
                sessionsv5(i).matchColor = matchColors(6,:);
            case 'dk052'
                sessionsv5(i).matchColor = matchColors(7,:);
            case 'dk063'
                sessionsv5(i).matchColor = matchColors(8,:);
            case 'dk070'
                sessionsv5(i).matchColor = matchColors(9,:);
            case 'dk138'
                sessionsv5(i).matchColor = matchColors(10,:);
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

    for i = 1:numel(sessionsv5)
        sets = sessionsv5(i).sets;

        sessionsv5(i).nTrialsBase = sum(sessionsv5(i).iTrialsReal(1:sets.changetrial));
        sessionsv5(i).nTrialsEarly = sum(sessionsv5(i).iTrialsReal(sets.changetrial+1:sets.changetrial+adaptTrialGrp));
        sessionsv5(i).nTrialsLate = sum(sessionsv5(i).iTrialsReal(sets.changetrial+121-adaptTrialGrp:sets.changetrial+120));
        %         sessionsv5(i).nTrialsBase = sum(cat(1,trials(1:sets.changetrial).realtrial));
        %         sessionsv5(i).nTrialsEarly = sum(cat(1,trials(sets.changetrial+1:sets.changetrial+adaptTrialGrp).realtrial));
        %         sessionsv5(i).nTrialsLate = sum(cat(1,trials(sets.changetrial+121-adaptTrialGrp:sets.changetrial+120).realtrial));

        if ~sets.changeback
            sessionsv5(i).nTrialsWash = 0;
        else
            nWashTrials = sessionsv5(i).nTrials-sets.changebacktrial;
            if nWashTrials > adaptTrialGrp; nWashTrials = adaptTrialGrp; end
            sessionsv5(i).nTrialsWash = sum(sessionsv5(i).iTrialsReal(sets.changebacktrial+1:sets.changebacktrial+nWashTrials));
            %             sessionsv5(i).nTrialsWash = sum(cat(1,trials(sets.changebacktrial+1:sets.changebacktrial+nWashTrials).realtrial));
        end

        if sessionsv5(i).nTrialsReal > 0
            sessionsv5(i).corrPerc = sum(sessionsv5(i).corrBool)/sessionsv5(i).nTrialsReal;
            sessionsv5(i).corrPercBase = sum(sessionsv5(i).corrBool(1:sets.changetrial))/sessionsv5(i).nTrialsBase;
            sessionsv5(i).corrPercEarly = sum(sessionsv5(i).corrBool(sets.changetrial+1:sets.changetrial+adaptTrialGrp))/sessionsv5(i).nTrialsEarly;
            sessionsv5(i).corrPercLate = sum(sessionsv5(i).corrBool(sets.changetrial+121-adaptTrialGrp:sets.changetrial+120))/sessionsv5(i).nTrialsLate;
            if ~sets.changeback
                sessionsv5(i).corrPercWash = nan;
            else
                sessionsv5(i).corrPercWash = sum(sessionsv5(i).corrBool(sets.changebacktrial+1:sets.changebacktrial+nWashTrials))/sessionsv5(i).nTrialsWash;
            end
            xEndAll = nan(size(sessionsv5(i).iTrialsReal)); xEndAll(sessionsv5(i).params.trialinds) = sessionsv5(i).params.xendall;

            sessionsv5(i).xEndMu_Base = nanmean(xEndAll(1:sets.changetrial));
            sessionsv5(i).xEndMed_Base = nanmedian(xEndAll(1:sets.changetrial));
            sessionsv5(i).xEndMu_Early = nanmean(xEndAll(sets.changetrial+1:sets.changetrial+adaptTrialGrp));
            sessionsv5(i).xEndMed_Early = nanmedian(xEndAll(sets.changetrial+1:sets.changetrial+adaptTrialGrp));
            sessionsv5(i).xEndMu_Late = nanmean(xEndAll(sets.changetrial+121-adaptTrialGrp:sets.changetrial+120));
            sessionsv5(i).xEndMed_Late = nanmedian(xEndAll(sets.changetrial+121-adaptTrialGrp:sets.changetrial+120));

            if ~sets.changeback
                sessionsv5(i).xEndMu_Wash = nan;
                sessionsv5(i).xEndMed_Wash = nan;
            else
                sessionsv5(i).xEndMu_Wash = nanmean(xEndAll(sets.changebacktrial+1:sets.changebacktrial+nWashTrials));
                sessionsv5(i).xEndMed_Wash = nanmedian(xEndAll(sets.changebacktrial+1:sets.changebacktrial+nWashTrials));
            end

        else
            sessionsv5(i).corrPerc = 0;
        end
    end

end

clear datastruct gaindiv i inits j nWashTrials sets smspan trials v1 x0 x1 xEndAll

%% Make plots

sessionsv4_1 = sessionsv4(cat(1,sessionsv4.day1));
sessionsv4_n = sessionsv4(cat(1,sessionsv4.dayN));

% First figure: fraction of trials with movement
figout(1) = figure;
subplot(2,2,1); hold on;
smov(1) = scatter(rand(size(sessionsv1))/5+0.9,cat(1,sessionsv1.movPerc),72,cat(1,sessionsv1.matchColor)/2,'filled');
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

% First figure: fraction of trials with movement
figout(2) = figure;
subplot(2,2,1); hold on;
bar(1,mean(cat(1,sessionsv1.movPerc)),0.6,'EdgeColor','none','FaceColor',[.6 0 0],'FaceAlpha',0.5);
errorbar(1, mean(cat(1,sessionsv1.movPerc)), std(cat(1,sessionsv1.movPerc))/sqrt(numel(sessionsv1)),'k','CapSize',0,'LineWidth',2);
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
    print(figout(1),fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');

    figname = sprintf('percMovCorrect_%i_trials_bar',adaptTrialGrp);
    savefig(figout(2), fullfile(figFold,[figname,'.fig']));
    saveas(figout(2),fullfile(figFold,[figname,'.png']));
    print(figout(2),fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');

end

if closebool
    close all
end
end