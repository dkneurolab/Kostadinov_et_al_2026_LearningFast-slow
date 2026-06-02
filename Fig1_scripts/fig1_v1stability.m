function sessionsv1 = fig1_v1stability(sessionsv1,dataFold,paperFold,savebool)

% Set directories for behavioural data and output figures
behavFold = fullfile(dataFold,'04_Behav-processed');
figFold = fullfile(paperFold,'Figures','Fig1','Sub1_v1Learning');

%%
for iSesh = 1:numel(sessionsv1)
    dataDir = dir(fullfile(behavFold,sessionsv1(iSesh).name,sessionsv1(iSesh).date,'Behav*',[sessionsv1(iSesh).name,'_',sessionsv1(iSesh).date,'.mat']));
    datastruct = load(fullfile(dataDir.folder,dataDir.name)); datastruct = datastruct.datastruct;
    trials = datastruct.paqdata.behav.trials; trials = trials(cat(1,trials.trialside) > 0);
    inits = datastruct.inits;
    sets = datastruct.sets;
    
    gaindiv = sets.start1/360*4.5*pi/sets.rotgain;
    if sets.vrfs == 60
        smspan = inits.paqfs*.04; % 40 ms smoothing span
    elseif sets. vrfs == 30
        smspan = inits.paqfs*.06; % 60 ms smoothing span
    end
    
    movvAll = nan*zeros(4*inits.paqfs/50,numel(trials));
    xStopAll = nan*zeros(numel(trials),1);
    tRew = nan*zeros(numel(trials),1);
    for jTrial = 1:numel(trials)        
        v0 = trials(jTrial).data(inits.paqfs*inits.mintime+1:end-inits.paqfs*inits.mintime,end)*-gaindiv;
        if numel(v0) > 15*inits.paqfs; v0 = v0(1:15*inits.paqfs); end
        iMov = threshdetect_dk(v0,1.5,1);
        if ~isempty(iMov)
            iMov = iMov(1)+inits.paqfs*inits.mintime;
            x0 = trials(jTrial).data(:,10);
            x0start = x0(inits.paqfs*inits.mintime+1:end);
            if any(x0start < -3) && any(x0start > 3)
                v1start = diff(x0start);
                wrapInds = find(abs(v1start) > 3);
                for iWraps = 1:numel(wrapInds)
                    v1start(wrapInds(iWraps)) = mean([v1start(wrapInds(iWraps)-1) v1start(wrapInds(iWraps)+1)]);
                end
                x0start = [x0start(1); cumsum(v1start)+x0start(1)];
                x0 = [x0(1:inits.paqfs*inits.mintime); x0start];
            end
            x1 = smooth(x0,smspan/numel(x0),'lowess');
            v1 = diff(x1)*inits.paqfs; v1 = smooth([v1; v1(end)],100)*-gaindiv;
            %             v1 = trials(jTrial).data(:,end)*-gaindiv;
            movv = v1(iMov-2*inits.paqfs+1:iMov+2*inits.paqfs);
            movvAll(:,jTrial) = decimate(double(movv),50);
            
            % Find movement offset
            movSlow = abs(movv) < .5;
            movSlowMu = movmean(movSlow,inits.paqfs/10);
            iStop = threshdetect_dk(movSlowMu(end/2+1:end),1,1);
            if ~isempty(iStop)
                xStopAll(jTrial) = x1(iMov+iStop(1));                 
            end
        end
        
        % was there a reward and was it a correct trial?
        trialDur = diff(trials(jTrial).inds)/inits.paqfs;
        if trialDur < sets.responsetime & trials(jTrial).reward
            tRew(jTrial) = trials(jTrial).inds(2)/inits.paqfs;        
        end
        
    end
    
    sessionsv1(iSesh).leftTrials = cat(1,datastruct.paqdata.behav.trials.trialside) > 0;
    sessionsv1(iSesh).movTrials = ~isnan(movvAll(1,:)');
    sessionsv1(iSesh).percMov = sum(sessionsv1(iSesh).movTrials)/numel(sessionsv1(iSesh).movTrials);
    sessionsv1(iSesh).movAll = movvAll(151:300,sessionsv1(iSesh).movTrials);
    sessionsv1(iSesh).xStopTrials = ~isnan(xStopAll);
    sessionsv1(iSesh).xStopAll = xStopAll(~isnan(xStopAll));
    sessionsv1(iSesh).xStart = sets.start1;
    sessionsv1(iSesh).xGain = sets.rotgain;
    sessionsv1(iSesh).rewTrials = ~isnan(tRew);
    sessionsv1(iSesh).rewIntervals = diff(tRew(~isnan(tRew)));
    
    
end

if savebool
    if ~exist(figFold,'dir'); mkdir(figFold); end
    save(fullfile(figFold,'v1LearningCurves.mat'),'sessionsv1');
end