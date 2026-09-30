function [behavSesh,endProbCoarse,endProbFine] = fig1_adaptBehav_sortTrials(behavSesh,trialInds,trialsGood,behavStability,histBins_coarse,histBins_fine)

%% Concatenate endpos data
% Baseline
behavSesh.endPos.iTrialsBase        = trialInds.iTrialsBase;
behavSesh.endPos.endPosBase         = cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsBase)).endpos);
behavSesh.endPos.endPosRealBase     = cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsBase)).endposreal);
% Early
behavSesh.endPos.iTrialsEarly       = trialInds.iTrialsEarly;
behavSesh.endPos.endPosEarly        = cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsEarly)).endpos);
behavSesh.endPos.endPosRealEarly    = cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsEarly)).endposreal);
% Mid
behavSesh.endPos.iTrialsMid         = trialInds.iTrialsMid;
behavSesh.endPos.endPosMid          = cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsMid)).endpos);
behavSesh.endPos.endPosRealMid      = cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsMid)).endposreal);
% Late
behavSesh.endPos.iTrialsLate        = trialInds.iTrialsLate;
behavSesh.endPos.endPosLate         = cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsLate)).endpos);
behavSesh.endPos.endPosRealLate     = cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsLate)).endposreal);
% Washout - first few trials
behavSesh.endPos.iTrialsWash_early  = trialInds.iTrialsWash_early;
behavSesh.endPos.endPosWash_early   = cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsWash_early)).endpos);
behavSesh.endPos.endPosRealWash_early   = cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsWash_early)).endposreal);
% Washout - later trials
behavSesh.endPos.iTrialsWash_late   = trialInds.iTrialsWash_late;
behavSesh.endPos.endPosWash_late    = cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsWash_late)).endpos);
behavSesh.endPos.endPosRealWash_late    = cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsWash_late)).endposreal);


%% Collate outcomes
behavSesh.endPosMu = [mean(behavSesh.endPos.endPosBase); mean(behavSesh.endPos.endPosEarly); mean(behavSesh.endPos.endPosMid); ...
    mean(behavSesh.endPos.endPosLate); mean(behavSesh.endPos.endPosWash_early); mean(behavSesh.endPos.endPosWash_late)];
behavSesh.endPosRealMu = [mean(behavSesh.endPos.endPosRealBase); mean(behavSesh.endPos.endPosRealEarly); mean(behavSesh.endPos.endPosRealMid); ...
    mean(behavSesh.endPos.endPosRealLate); mean(behavSesh.endPos.endPosRealWash_early); mean(behavSesh.endPos.endPosRealWash_late)];
behavSesh.percC = 100*[mean(cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsBase)).correct));...
    mean(cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsEarly)).correct));...
    mean(cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsMid)).correct));...
    mean(cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsLate)).correct));...
    mean(cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsWash_early)).correct));...
    mean(cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsWash_late)).correct))];

%% Bin endpos for histograms

endProbCoarse = zeros(numel(histBins_coarse)-1,6);
endProbFine = zeros(numel(histBins_fine)-1,6);

% Coarse binning
endProbCoarse(:,1) = 100*histcounts(-behavSesh.endPos.endPosBase,histBins_coarse,'Normalization','probability')';
endProbCoarse(:,2) = 100*histcounts(-behavSesh.endPos.endPosEarly,histBins_coarse,'Normalization','probability')';
endProbCoarse(:,3) = 100*histcounts(-behavSesh.endPos.endPosMid,histBins_coarse,'Normalization','probability')';
endProbCoarse(:,4) = 100*histcounts(-behavSesh.endPos.endPosLate,histBins_coarse,'Normalization','probability')';
endProbCoarse(:,5) = 100*histcounts(-behavSesh.endPos.endPosWash_early,histBins_coarse,'Normalization','probability')'; 
endProbCoarse(:,6) = 100*histcounts(-behavSesh.endPos.endPosWash_late,histBins_coarse,'Normalization','probability')'; 

% Fine binning
endProbFine(:,1) = 100*histcounts(-behavSesh.endPos.endPosBase,histBins_fine,'Normalization','probability')';
endProbFine(:,2) = 100*histcounts(-behavSesh.endPos.endPosEarly,histBins_fine,'Normalization','probability')';
endProbFine(:,3) = 100*histcounts(-behavSesh.endPos.endPosMid,histBins_fine,'Normalization','probability')';
endProbFine(:,4) = 100*histcounts(-behavSesh.endPos.endPosLate,histBins_fine,'Normalization','probability')';
endProbFine(:,5) = 100*histcounts(-behavSesh.endPos.endPosWash_early,histBins_fine,'Normalization','probability')';
endProbFine(:,6) = 100*histcounts(-behavSesh.endPos.endPosWash_late,histBins_fine,'Normalization','probability')'; 

%% Concatenate positional traces and movement velocities
% Baseline - no need to splil real and adapt
behavSesh.endTrace.endTraceBase = behavStability.offxall(:,ismember(trialInds.iTrialsGood,trialInds.iTrialsBase));
% Early - split real and adapt
endTraces0 = behavStability.offxall(:,ismember(trialInds.iTrialsGood,trialInds.iTrialsEarly));
behavSesh.endTrace.endTraceEarly_real = endTraces0;
for iTrace = 1:size(endTraces0,2)
    endTrace = endTraces0(:,iTrace);
    endVel = diff(endTrace)*behavStability.adaptGain;
    behavSesh.endTrace.endTraceEarly(:,iTrace) = [endTrace(1); cumsum(endVel)+endTrace(1)];
end
% Mid - split real and adapt
endTraces0 = behavStability.offxall(:,ismember(trialInds.iTrialsGood,trialInds.iTrialsMid));
behavSesh.endTrace.endTraceMid_real = endTraces0;
for iTrace = 1:size(endTraces0,2)
    endTrace = endTraces0(:,iTrace);
    endVel = diff(endTrace)*behavStability.adaptGain;
    behavSesh.endTrace.endTraceMid(:,iTrace) = [endTrace(1); cumsum(endVel)+endTrace(1)];
end
% Late - split real and adapt
endTraces0 = behavStability.offxall(:,ismember(trialInds.iTrialsGood,trialInds.iTrialsLate));
behavSesh.endTrace.endTraceLate_real = endTraces0;
for iTrace = 1:size(endTraces0,2)
    endTrace = endTraces0(:,iTrace);
    endVel = diff(endTrace)*behavStability.adaptGain;
    behavSesh.endTrace.endTraceLate(:,iTrace) = [endTrace(1); cumsum(endVel)+endTrace(1)];
end
% Wash early - no need to splil real and adapt
behavSesh.endTrace.endTraceWash_early = behavStability.offxall(:,ismember(trialInds.iTrialsGood,trialInds.iTrialsWash_early));
if isempty(behavSesh.endTrace.endTraceWash_early); behavSesh.endTrace.endTraceWash_early = nan(size(behavSesh.endTrace.endTraceWash_early,1),1); end
% Wash late - no need to split
behavSesh.endTrace.endTraceWash_late = behavStability.offxall(:,ismember(trialInds.iTrialsGood,trialInds.iTrialsWash_late));
if isempty(behavSesh.endTrace.endTraceWash_late); behavSesh.endTrace.endTraceWash_late = nan(size(behavSesh.endTrace.endTraceWash_late,1),1); end

behavSesh.endTrace.endTraceMu = [mean(behavSesh.endTrace.endTraceBase(1:end/2,:),2,'omitnan'); mean(behavSesh.endTrace.endTraceEarly(1:end/2,:),2,'omitnan'); ...
    mean(behavSesh.endTrace.endTraceMid(1:end/2,:),2,'omitnan'); mean(behavSesh.endTrace.endTraceLate(1:end/2,:),2,'omitnan'); ...
    mean(behavSesh.endTrace.endTraceWash_early(1:end/2,:),2,'omitnan'); mean(behavSesh.endTrace.endTraceWash_late(1:end/2,:),2,'omitnan')];
behavSesh.endTrace.endTraceSD = [std(behavSesh.endTrace.endTraceBase(1:end/2,:),[],2,'omitnan'); std(behavSesh.endTrace.endTraceEarly(1:end/2,:),[],2,'omitnan'); ...
    std(behavSesh.endTrace.endTraceMid(1:end/2,:),[],2,'omitnan'); std(behavSesh.endTrace.endTraceLate(1:end/2,:),[],2,'omitnan'); ...
    std(behavSesh.endTrace.endTraceWash_early(1:end/2,:),[],2,'omitnan'); std(behavSesh.endTrace.endTraceWash_late(1:end/2,:),[],2,'omitnan')];
behavSesh.endTrace.endTraceSEM = [std(behavSesh.endTrace.endTraceBase(1:end/2,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsBase)); ...
    std(behavSesh.endTrace.endTraceEarly(1:end/2,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsEarly)); ...
    std(behavSesh.endTrace.endTraceMid(1:end/2,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsMid)); ...
    std(behavSesh.endTrace.endTraceLate(1:end/2,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsLate)); ...
    std(behavSesh.endTrace.endTraceWash_early(1:end/2,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsWash_early)); ...
    std(behavSesh.endTrace.endTraceWash_late(1:end/2,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsWash_late))];
behavSesh.endTrace.endTraceRealMu = [mean(behavSesh.endTrace.endTraceBase(1:end/2,:),2,'omitnan'); mean(behavSesh.endTrace.endTraceEarly_real(1:end/2,:),2,'omitnan'); ...
    mean(behavSesh.endTrace.endTraceMid_real(1:end/2,:),2,'omitnan'); mean(behavSesh.endTrace.endTraceLate_real(1:end/2,:),2,'omitnan'); ...
    mean(behavSesh.endTrace.endTraceWash_early(1:end/2,:),2,'omitnan'); mean(behavSesh.endTrace.endTraceWash_late(1:end/2,:),2,'omitnan')];
behavSesh.endTrace.endTraceRealSD = [std(behavSesh.endTrace.endTraceBase(1:end/2,:),[],2,'omitnan'); std(behavSesh.endTrace.endTraceEarly_real(1:end/2,:),[],2,'omitnan'); ...
    std(behavSesh.endTrace.endTraceMid_real(1:end/2,:),[],2,'omitnan'); std(behavSesh.endTrace.endTraceLate_real(1:end/2,:),[],2,'omitnan'); ...
    std(behavSesh.endTrace.endTraceWash_early(1:end/2,:),[],2,'omitnan'); std(behavSesh.endTrace.endTraceWash_late(1:end/2,:),[],2,'omitnan')];
behavSesh.endTrace.endTraceRealSEM = [std(behavSesh.endTrace.endTraceBase(1:end/2,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsBase)); ...
    std(behavSesh.endTrace.endTraceEarly_real(1:end/2,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsEarly)); ...
    std(behavSesh.endTrace.endTraceMid_real(1:end/2,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsMid)); ...
    std(behavSesh.endTrace.endTraceLate_real(1:end/2,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsLate)); ...
    std(behavSesh.endTrace.endTraceWash_early(1:end/2,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsWash_early)); ...
    std(behavSesh.endTrace.endTraceWash_late(1:end/2,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsWash_late))];



behavSesh.endTrace.movVelBase = behavStability.movvall(:,ismember(trialInds.iTrialsGood,trialInds.iTrialsBase))*-behavStability.gainDiv;
behavSesh.endTrace.movVelEarly = behavStability.movvall(:,ismember(trialInds.iTrialsGood,trialInds.iTrialsEarly))*-behavStability.gainDiv;
behavSesh.endTrace.movVelMid = behavStability.movvall(:,ismember(trialInds.iTrialsGood,trialInds.iTrialsMid))*-behavStability.gainDiv;
behavSesh.endTrace.movVelLate = behavStability.movvall(:,ismember(trialInds.iTrialsGood,trialInds.iTrialsLate))*-behavStability.gainDiv;
behavSesh.endTrace.movVelWash_early = behavStability.movvall(:,ismember(trialInds.iTrialsGood,trialInds.iTrialsWash_early))*-behavStability.gainDiv;
behavSesh.endTrace.movVelWash_late = behavStability.movvall(:,ismember(trialInds.iTrialsGood,trialInds.iTrialsWash_late))*-behavStability.gainDiv;
if isempty(behavSesh.endTrace.movVelWash_early); behavSesh.endTrace.movVelWash_early = nan(size(behavSesh.endTrace.movVelWash_early,1),1); end
if isempty(behavSesh.endTrace.movVelWash_late); behavSesh.endTrace.movVelWash_late = nan(size(behavSesh.endTrace.movVelWash_late,1),1); end

behavSesh.endTrace.movVelMu = [mean(behavSesh.endTrace.movVelBase(end/4+1:end/4*3,:),2,'omitnan'); mean(behavSesh.endTrace.movVelEarly(end/4+1:end/4*3,:),2,'omitnan'); ...
    mean(behavSesh.endTrace.movVelMid(end/4+1:end/4*3,:),2,'omitnan'); mean(behavSesh.endTrace.movVelLate(end/4+1:end/4*3,:),2,'omitnan'); ...
    mean(behavSesh.endTrace.movVelWash_early(end/4+1:end/4*3,:),2,'omitnan'); mean(behavSesh.endTrace.movVelWash_late(end/4+1:end/4*3,:),2,'omitnan')];
behavSesh.endTrace.movVelSD = [std(behavSesh.endTrace.movVelBase(end/4+1:end/4*3,:),[],2,'omitnan'); std(behavSesh.endTrace.movVelEarly(end/4+1:end/4*3,:),[],2,'omitnan'); ...
    std(behavSesh.endTrace.movVelMid(end/4+1:end/4*3,:),[],2,'omitnan'); std(behavSesh.endTrace.movVelLate(end/4+1:end/4*3,:),[],2,'omitnan'); ...
    std(behavSesh.endTrace.movVelWash_early(end/4+1:end/4*3,:),[],2,'omitnan'); std(behavSesh.endTrace.movVelWash_late(end/4+1:end/4*3,:),[],2,'omitnan')];
behavSesh.endTrace.movVelSEM = [std(behavSesh.endTrace.movVelBase(end/4+1:end/4*3,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsBase)); ...
    std(behavSesh.endTrace.movVelEarly(end/4+1:end/4*3,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsEarly)); ...
    std(behavSesh.endTrace.movVelMid(end/4+1:end/4*3,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsMid)); ...
    std(behavSesh.endTrace.movVelLate(end/4+1:end/4*3,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsLate)); ...
    std(behavSesh.endTrace.movVelWash_early(end/4+1:end/4*3,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsWash_early)); ...
    std(behavSesh.endTrace.movVelWash_late(end/4+1:end/4*3,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsWash_late))];


end