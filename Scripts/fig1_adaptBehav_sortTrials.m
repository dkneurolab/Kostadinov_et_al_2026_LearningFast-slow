function [behavSesh,endProbCoarse,endProbFine] = fig1_adaptBehav_sortTrials(behavSesh,trialInds,trialsGood,behavStability,histBins_coarse,histBins_fine)

%% Concatenate endpos data
% Baseline
behavSesh.endPos.iTrialsBase = trialInds.iTrialsBase;
behavSesh.endPos.endPosBase = cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsBase)).endpos);
behavSesh.endPos.endPosRealBase = cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsBase)).endposreal);
% Early
behavSesh.endPos.iTrialsEarly = trialInds.iTrialsEarly;
behavSesh.endPos.endPosEarly = cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsEarly)).endpos);
behavSesh.endPos.endPosRealEarly = cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsEarly)).endposreal);
% Mid
behavSesh.endPos.iTrialsMid = trialInds.iTrialsMid;
behavSesh.endPos.endPosMid = cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsMid)).endpos);
behavSesh.endPos.endPosRealMid = cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsMid)).endposreal);
% Late
behavSesh.endPos.iTrialsLate = trialInds.iTrialsLate;
behavSesh.endPos.endPosLate = cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsLate)).endpos);
behavSesh.endPos.endPosRealLate = cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsLate)).endposreal);
% Washout
behavSesh.endPos.iTrialsWash = trialInds.iTrialsWash;
behavSesh.endPos.endPosWash = cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsWash)).endpos);
behavSesh.endPos.endPosRealWash = cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsWash)).endposreal);

%% Collate outcomes
behavSesh.endPosMu = [mean(behavSesh.endPos.endPosBase); mean(behavSesh.endPos.endPosEarly); mean(behavSesh.endPos.endPosMid); ...
    mean(behavSesh.endPos.endPosLate); mean(behavSesh.endPos.endPosWash)];
behavSesh.endPosRealMu = [mean(behavSesh.endPos.endPosRealBase); mean(behavSesh.endPos.endPosRealEarly); mean(behavSesh.endPos.endPosRealMid); ...
    mean(behavSesh.endPos.endPosRealLate); mean(behavSesh.endPos.endPosRealWash)];
behavSesh.percC = 100*[mean(cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsBase)).correct));...
    mean(cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsEarly)).correct));...
    mean(cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsMid)).correct));...
    mean(cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsLate)).correct));...
    mean(cat(1,trialsGood(ismember(trialInds.iTrialsGood,trialInds.iTrialsWash)).correct))];

%% Bin endpos for histograms

endProbCoarse = zeros(numel(histBins_coarse)-1,5);
endProbFine = zeros(numel(histBins_fine)-1,5);

% Coarse binning
endProbCoarse(:,1) = 100*histcounts(-behavSesh.endPos.endPosBase,histBins_coarse,'Normalization','probability')';
endProbCoarse(:,2) = 100*histcounts(-behavSesh.endPos.endPosEarly,histBins_coarse,'Normalization','probability')';
endProbCoarse(:,3) = 100*histcounts(-behavSesh.endPos.endPosMid,histBins_coarse,'Normalization','probability')';
endProbCoarse(:,4) = 100*histcounts(-behavSesh.endPos.endPosLate,histBins_coarse,'Normalization','probability')';
endProbCoarse(:,5) = 100*histcounts(-behavSesh.endPos.endPosWash,histBins_coarse,'Normalization','probability')'; 

% Fine binning
endProbFine(:,1) = 100*histcounts(-behavSesh.endPos.endPosBase,histBins_fine,'Normalization','probability')';
endProbFine(:,2) = 100*histcounts(-behavSesh.endPos.endPosEarly,histBins_fine,'Normalization','probability')';
endProbFine(:,3) = 100*histcounts(-behavSesh.endPos.endPosMid,histBins_fine,'Normalization','probability')';
endProbFine(:,4) = 100*histcounts(-behavSesh.endPos.endPosLate,histBins_fine,'Normalization','probability')';
endProbFine(:,5) = 100*histcounts(-behavSesh.endPos.endPosWash,histBins_fine,'Normalization','probability')';

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
% Wash - no need to splil real and adapt
behavSesh.endTrace.endTraceWash = behavStability.offxall(:,ismember(trialInds.iTrialsGood,trialInds.iTrialsWash));
if isempty(behavSesh.endTrace.endTraceWash); behavSesh.endTrace.endTraceWash = nan(size(behavSesh.endTrace.endTraceWash,1),1); end

behavSesh.endTrace.endTraceMu = [mean(behavSesh.endTrace.endTraceBase(1:end/2,:),2,'omitnan'); mean(behavSesh.endTrace.endTraceEarly(1:end/2,:),2,'omitnan'); ...
    mean(behavSesh.endTrace.endTraceMid(1:end/2,:),2,'omitnan'); mean(behavSesh.endTrace.endTraceLate(1:end/2,:),2,'omitnan'); mean(behavSesh.endTrace.endTraceWash(1:end/2,:),2,'omitnan')];
behavSesh.endTrace.endTraceSD = [std(behavSesh.endTrace.endTraceBase(1:end/2,:),[],2,'omitnan'); std(behavSesh.endTrace.endTraceEarly(1:end/2,:),[],2,'omitnan'); ...
    std(behavSesh.endTrace.endTraceMid(1:end/2,:),[],2,'omitnan'); std(behavSesh.endTrace.endTraceLate(1:end/2,:),[],2,'omitnan'); std(behavSesh.endTrace.endTraceWash(1:end/2,:),[],2,'omitnan')];
behavSesh.endTrace.endTraceSEM = [std(behavSesh.endTrace.endTraceBase(1:end/2,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsBase)); ...
    std(behavSesh.endTrace.endTraceEarly(1:end/2,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsEarly)); ...
    std(behavSesh.endTrace.endTraceMid(1:end/2,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsMid)); ...
    std(behavSesh.endTrace.endTraceLate(1:end/2,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsLate)); ...
    std(behavSesh.endTrace.endTraceWash(1:end/2,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsWash))];
behavSesh.endTrace.endTraceRealMu = [mean(behavSesh.endTrace.endTraceBase(1:end/2,:),2,'omitnan'); mean(behavSesh.endTrace.endTraceEarly_real(1:end/2,:),2,'omitnan'); ...
    mean(behavSesh.endTrace.endTraceMid_real(1:end/2,:),2,'omitnan'); mean(behavSesh.endTrace.endTraceLate_real(1:end/2,:),2,'omitnan'); mean(behavSesh.endTrace.endTraceWash(1:end/2,:),2,'omitnan')];
behavSesh.endTrace.endTraceRealSD = [std(behavSesh.endTrace.endTraceBase(1:end/2,:),[],2,'omitnan'); std(behavSesh.endTrace.endTraceEarly_real(1:end/2,:),[],2,'omitnan'); ...
    std(behavSesh.endTrace.endTraceMid_real(1:end/2,:),[],2,'omitnan'); std(behavSesh.endTrace.endTraceLate_real(1:end/2,:),[],2,'omitnan'); std(behavSesh.endTrace.endTraceWash(1:end/2,:),[],2,'omitnan')];
behavSesh.endTrace.endTraceRealSEM = [std(behavSesh.endTrace.endTraceBase(1:end/2,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsBase)); ...
    std(behavSesh.endTrace.endTraceEarly_real(1:end/2,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsEarly)); ...
    std(behavSesh.endTrace.endTraceMid_real(1:end/2,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsMid)); ...
    std(behavSesh.endTrace.endTraceLate_real(1:end/2,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsLate)); ...
    std(behavSesh.endTrace.endTraceWash(1:end/2,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsWash))];



behavSesh.endTrace.movVelBase = behavStability.movvall(:,ismember(trialInds.iTrialsGood,trialInds.iTrialsBase))*-behavStability.gainDiv;
behavSesh.endTrace.movVelEarly = behavStability.movvall(:,ismember(trialInds.iTrialsGood,trialInds.iTrialsEarly))*-behavStability.gainDiv;
behavSesh.endTrace.movVelMid = behavStability.movvall(:,ismember(trialInds.iTrialsGood,trialInds.iTrialsMid))*-behavStability.gainDiv;
behavSesh.endTrace.movVelLate = behavStability.movvall(:,ismember(trialInds.iTrialsGood,trialInds.iTrialsLate))*-behavStability.gainDiv;
behavSesh.endTrace.movVelWash = behavStability.movvall(:,ismember(trialInds.iTrialsGood,trialInds.iTrialsWash))*-behavStability.gainDiv;
if isempty(behavSesh.endTrace.movVelWash); behavSesh.endTrace.movVelWash = nan(size(behavSesh.endTrace.movVelWash,1),1); end

behavSesh.endTrace.movVelMu = [mean(behavSesh.endTrace.movVelBase(end/4+1:end/4*3,:),2,'omitnan'); mean(behavSesh.endTrace.movVelEarly(end/4+1:end/4*3,:),2,'omitnan'); ...
    mean(behavSesh.endTrace.movVelMid(end/4+1:end/4*3,:),2,'omitnan'); mean(behavSesh.endTrace.movVelLate(end/4+1:end/4*3,:),2,'omitnan'); mean(behavSesh.endTrace.movVelWash(end/4+1:end/4*3,:),2,'omitnan')];
behavSesh.endTrace.movVelSD = [std(behavSesh.endTrace.movVelBase(end/4+1:end/4*3,:),[],2,'omitnan'); std(behavSesh.endTrace.movVelEarly(end/4+1:end/4*3,:),[],2,'omitnan'); ...
    std(behavSesh.endTrace.movVelMid(end/4+1:end/4*3,:),[],2,'omitnan'); std(behavSesh.endTrace.movVelLate(end/4+1:end/4*3,:),[],2,'omitnan'); std(behavSesh.endTrace.movVelWash(end/4+1:end/4*3,:),[],2,'omitnan')];
behavSesh.endTrace.movVelSEM = [std(behavSesh.endTrace.movVelBase(end/4+1:end/4*3,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsBase)); ...
    std(behavSesh.endTrace.movVelEarly(end/4+1:end/4*3,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsEarly)); ...
    std(behavSesh.endTrace.movVelMid(end/4+1:end/4*3,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsMid)); ...
    std(behavSesh.endTrace.movVelLate(end/4+1:end/4*3,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsLate)); ...
    std(behavSesh.endTrace.movVelWash(end/4+1:end/4*3,:),[],2,'omitnan')/sqrt(numel(trialInds.iTrialsWash))];


% behavSesh.endTrace.movVelRealMu = [mean(behavSesh.endTrace.movVelBase(end/4+1:end/4*3,:),2,'omitnan'); mean(behavSesh.endTrace.movVelEarly(end/4+1:end/4*3,:),2,'omitnan'); ...
%     mean(behavSesh.endTrace.movVelMid(end/4+1:end/4*3,:),2,'omitnan'); mean(behavSesh.endTrace.movVelLate(end/4+1:end/4*3,:),2,'omitnan'); mean(behavSesh.endTrace.movVelWash(end/4+1:end/4*3,:),2,'omitnan')];
% behavSesh.endTrace.movVelRealSD = [std(behavSesh.endTrace.movVelBase(end/4+1:end/4*3,:),[],2,'omitnan'); std(behavSesh.endTrace.movVelEarly(end/4+1:end/4*3,:),[],2,'omitnan'); ...
%     std(behavSesh.endTrace.movVelMid(end/4+1:end/4*3,:),[],2,'omitnan'); std(behavSesh.endTrace.movVelLate(end/4+1:end/4*3,:),[],2,'omitnan'); std(behavSesh.endTrace.movVelWash(end/4+1:end/4*3,:),[],2,'omitnan')];













end