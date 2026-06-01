function matchData = fig1_roiMatch_load(matchData,alignFold,imFold,behavFold,sessionpair)
%% Load in ROI matching info

try
    pairDir = dir(fullfile(imFold,matchData.name,sprintf('REG_%s_%s',erase(matchData.date1(end-4:end),'-'),erase(matchData.date2(end-4:end),'-'))));
    alignData = dir(fullfile(pairDir(1).folder,'*_analysis.mat'));
    alignData = load(fullfile(alignData.folder,alignData.name));
    alignData = alignData.regi;
    indOrder = [1,2];    
catch
    pairDir = dir(fullfile(imFold,matchData.name,sprintf('REG_%s_%s',erase(matchData.date2(end-4:end),'-'),erase(matchData.date1(end-4:end),'-'))));
    alignData = dir(fullfile(pairDir(1).folder,'*_analysis.mat'));
    alignData = load(fullfile(alignData.folder,alignData.name));
    alignData = alignData.regi;
    indOrder = [2,1];    
    bb = matchData.nRoi1;
    matchData.nRoi1 = matchData.nRoi2;
    matchData.nRoi2 = bb;
end

tForm = dir(fullfile(pairDir(1).folder,'*_transform.mat'));
tForm = load(fullfile(tForm.folder,tForm.name));

matchData.iROIs1 = alignData.rois.idcs(:,indOrder(1));
matchData.iROIs2 = alignData.rois.idcs(:,indOrder(2));

%% Load datastructs and RFs
dataFold1 = fullfile(alignFold,sprintf('Version_%i',sessionpair(1).version),...
    sessionpair(1).name,sprintf('%s_%s',sessionpair(1).date,sessionpair(1).fov));
dataStruct1 = load(fullfile(dataFold1,'eventsdata.mat'));
spkmat1 = load(fullfile(dataFold1,'spkmat.mat'));
imStruct1 = dir(fullfile(imFold,sessionpair(1).name,sessionpair(1).date,sessionpair(1).fov,'*final.mat'));
imStruct1 = load(fullfile(imStruct1(numel(imStruct1)).folder,imStruct1(numel(imStruct1)).name));
imStruct1.ops = imStruct1.datdk.ops; imStruct1.final = imStruct1.datdk.final; 
imStruct1.zF = zscore(imStruct1.datdk.finalparams.dff_30(1:dataStruct1.sets.imfs*100,matchData.iROIs1));
imStruct1.spks = spkmat1.spkmat(1:dataStruct1.sets.imfs*100,matchData.iROIs1); imStruct1.final(1).xycoords = [];
imStruct1 = rmfield(imStruct1,'datdk'); imStruct1.final = rmfield(imStruct1.final,{'F','rois','xycoords'});
spkmat0 = spkmat1.spkmat; spkmat0(spkmat0(:) < 0.5) = nan;
spksize = mean(spkmat0,'omitnan'); spksize(isnan(spksize)) = 1;
spkmat1 = spkmat1.spkmat./spksize;
for iIm = 1:numel(imStruct1.final)
    imStruct1.final(iIm).xpix = int32(imStruct1.final(iIm).xpix);
    imStruct1.final(iIm).ypix = int32(imStruct1.final(iIm).ypix);
end
rfs1 = load(fullfile(dataFold1,'RFs','RFs.mat')); rfs1 = rfs1.RFs;
behavFold1 = fullfile(behavFold,sessionpair(1).name,sessionpair(1).date,sprintf('Behaviour_%is',dataStruct1.sets.mintime));
behavStruct1 = load(fullfile(behavFold1,sprintf('%s_%s.mat',sessionpair(1).name,sessionpair(1).date))); behavStruct1 = behavStruct1.datastruct;

dataFold2 = fullfile(alignFold,sprintf('Version_%i',sessionpair(2).version),...
    sessionpair(2).name,sprintf('%s_%s',sessionpair(2).date,sessionpair(2).fov));
dataStruct2 = load(fullfile(dataFold2,'eventsdata.mat'));
spkmat2 = load(fullfile(dataFold2,'spkmat.mat'));
imStruct2 = dir(fullfile(imFold,sessionpair(2).name,sessionpair(2).date,sessionpair(2).fov,'*final.mat'));
imStruct2 = load(fullfile(imStruct2(numel(imStruct2)).folder,imStruct2(numel(imStruct2)).name));
imStruct2.ops = imStruct2.datdk.ops; imStruct2.final = imStruct2.datdk.final;
imStruct2.zF = zscore(imStruct2.datdk.finalparams.dff_30(1:dataStruct2.sets.imfs*100,matchData.iROIs2));
imStruct2.spks = spkmat2.spkmat(1:dataStruct2.sets.imfs*100,matchData.iROIs2); imStruct2.final(1).xycoords = [];
imStruct2 = rmfield(imStruct2,'datdk'); imStruct2.final = rmfield(imStruct2.final,{'F','rois','xycoords'});
spkmat0 = spkmat2.spkmat; spkmat0(spkmat0(:) < 0.5) = nan;
spksize = mean(spkmat0,'omitnan'); spksize(isnan(spksize)) = 1;
spkmat2 = spkmat2.spkmat./spksize;
for iIm = 1:numel(imStruct2.final)
    imStruct2.final(iIm).xpix = int32(imStruct2.final(iIm).xpix);
    imStruct2.final(iIm).ypix = int32(imStruct2.final(iIm).ypix);
end
rfs2 = load(fullfile(dataFold2,'RFs','RFs.mat')); rfs2 = rfs2.RFs;
behavFold2 = fullfile(behavFold,sessionpair(2).name,sessionpair(2).date,sprintf('Behaviour_%is',dataStruct2.sets.mintime));
behavStruct2 = load(fullfile(behavFold2,sprintf('%s_%s.mat',sessionpair(2).name,sessionpair(2).date))); behavStruct2 = behavStruct2.datastruct;

matchData.imData = cat(2,imStruct1,imStruct2);
matchData.tForm = tForm.t; matchData.tForm.tform_direction = matchData.tForm.tform_direction(indOrder)+1;

%% Take relevant data from structures - onsets
% Use same smoothing for everyone:
smspan = dataStruct1.sets.fs*.06; % 60 ms smoothing span

imfs1 = dataStruct1.sets.imfs;
paqfs1 = dataStruct1.sets.fs;
trials1 = dataStruct1.savedata.data_trials;
gaindiv1 = behavStruct1.sets.start1/360*4.5*pi/dataStruct1.sets.vrgain;
matchData.on1.spk = single(trials1.spk.spkONL(end/2-2*imfs1+1:end/2+2*imfs1,:,:));
% if dataStruct1.sets.vrfs == 60
%     smspan = dataStruct1.sets.fs*.04; % 40 ms smoothing span
% elseif dataStruct1.sets.vrfs == 30
%     smspan = dataStruct1.sets.fs*.06; % 60 ms smoothing span
% end
for i = 1:size(trials1.VR.VR_FONL,3)
    matchData.on1.lick(:,i) = decimate(trials1.VR.VR_FONL(end/2-2*paqfs1+1:end/2+2*paqfs1,2,i),50);
    x0 = trials1.VR.VR_FONL(:,3,i);
    x1 = smooth(x0,smspan/numel(x0),'lowess');
    v1 = diff(x1)*dataStruct1.sets.fs*-gaindiv1; v1(end+1) = v1(end);
    matchData.on1.mov(:,i) = decimate(v1(end/2-2*paqfs1+1:end/2+2*paqfs1),50);
end

imfs2 = dataStruct2.sets.imfs;
paqfs2 = dataStruct2.sets.fs;
trials2 = dataStruct2.savedata.data_trials;
gaindiv2 = behavStruct2.sets.start1/360*4.5*pi/dataStruct2.sets.vrgain;
matchData.on2.spk = single(trials2.spk.spkONL(end/2-2*imfs2+1:end/2+2*imfs2,:,:));
for i = 1:size(trials2.VR.VR_FONL,3)
    matchData.on2.lick(:,i) = single(decimate(trials2.VR.VR_FONL(end/2-2*paqfs2+1:end/2+2*paqfs2,2,i),50));
    x0 = trials2.VR.VR_FONL(:,3,i);
    x1 = smooth(x0,smspan/numel(x0),'lowess');
    v1 = diff(x1)*dataStruct2.sets.fs*-gaindiv2; v1(end+1) = v1(end);
    matchData.on2.mov(:,i) = single(decimate(v1(end/2-2*paqfs2+1:end/2+2*paqfs2),50));
%     matchData.on2.mov(:,i) = decimate(trials2.VR.VR_FONL(end/2-2*paqfs2+1:end/2+2*paqfs2,4,i),50);
end

%% Take revelant data from structure - movements and trial reward
if isfield(trials1.spk,'spkLCmov')
    % If there is movement analysis data, use it - if there's one there's
    % all of them, so just stay in here for both movement and trial reward
    
    % Movement onset
    spkLocal = cat(3,trials1.spk.spkLCmov,trials1.spk.spkLUmov,trials1.spk.spkLOmov);
    behavLocal = cat(3,trials1.VR.VRLCmov,trials1.VR.VRLUmov,trials1.VR.VRLOmov);
    iTrialBad = isnan(squeeze(spkLocal(1,1,:)));
    spkLocal(:,:,iTrialBad) = []; behavLocal(:,:,iTrialBad) = [];
    matchData.mov1.spk = single(spkLocal(end/2-2*imfs1+1:end/2+2*imfs1,:,:));
    for iMov = 1:size(behavLocal,3)
        matchData.mov1.lick(:,iMov) = single(decimate(behavLocal(end/2-2*paqfs1+1:end/2+2*paqfs1,2,iMov),50));
        x0 = behavLocal(:,3,iMov);
        x1 = smooth(x0,smspan/numel(x0),'lowess');
        v1 = diff(x1)*dataStruct1.sets.fs*-gaindiv1; v1(end+1) = v1(end);
        matchData.mov1.mov(:,iMov) = single(decimate(v1(end/2-2*paqfs1+1:end/2+2*paqfs1),50));
    end
    
    % Trial reward
    spkLocal = trials1.spk.spkLC;
    behavLocal = trials1.VR.VR_FLC;
    iTrialBad = isnan(squeeze(spkLocal(1,1,:)));
    spkLocal(:,:,iTrialBad) = []; behavLocal(:,:,iTrialBad) = [];
    matchData.rewT1.spk = single(spkLocal(end/2-2*imfs1+1:end/2+2*imfs1,:,:));
    for iRew = 1:size(behavLocal,3)
        matchData.rewT1.lick(:,iRew) = single(decimate(behavLocal(end/2-2*paqfs1+1:end/2+2*paqfs1,2,iRew),50));
        x0 = behavLocal(:,3,iRew);
        x1 = smooth(x0,smspan/numel(x0),'lowess');
        v1 = diff(x1)*dataStruct1.sets.fs*-gaindiv1; v1(end+1) = v1(end);
        matchData.rewT1.mov(:,iRew) = single(decimate(v1(end/2-2*paqfs1+1:end/2+2*paqfs1),50));
    end
    
else
    % If there is no movement data, it's probably a v1 session and you need
    % to get the movement another way
    iFrames = [0;cumsum(sessionpair(1).fpblock)];
    
    % Use moveri to find trials with movement
    tr2p = behavStruct1.paqdata.imaging.trials2p;
    iMoveri = strcmpi(fieldnames(tr2p),'moveri');
    movTrials = struct2cell(tr2p); movTrials = squeeze(movTrials(iMoveri,:,:));
    tr2p = tr2p(~cellfun(@isempty,movTrials));
%     spkLocal = zeros(4*imfs1,matchData.nRoi1,numel(tr2p));
    behavLocal = zeros(4*paqfs1,2,numel(tr2p));
    for iMov = 1:numel(tr2p)
        iGlobal = tr2p(iMov).moveri(1)+tr2p(iMov).trialsession(1);
        iSpkLocal = tr2p(iMov).frames2p(find(tr2p(iMov).frames2p(:,2)>iGlobal,1),1);        
        iSpkGlobal = iFrames(tr2p(iMov).block2p)+iSpkLocal;
        matchData.mov1.spk(:,:,iMov) = single(spkmat1(iSpkGlobal-2*imfs1+1:iSpkGlobal+2*imfs1,:));
%         matchData.mov1.spk(:,:,iMov) = single(spkLocal(:,:,iMov));
        behavLocal(:,:,iMov) = tr2p(iMov).data(tr2p(iMov).moveri(1)-2*paqfs1+1:tr2p(iMov).moveri(1)+2*paqfs1,[6,10]);
        matchData.mov1.lick(:,iMov) = single(decimate(behavLocal(:,1,iMov),50));
        x0 = behavLocal(:,2,iMov);
        x1 = smooth(x0,smspan/numel(x0),'lowess');
        v1 = diff(x1)*dataStruct1.sets.fs*-gaindiv1; v1(end+1) = v1(end);
        matchData.mov1.mov(:,iMov) = single(decimate(v1(end/2-2*paqfs1+1:end/2+2*paqfs1),50));
    end
    
    % Use off2pLC to find trials with reward
    off2p = behavStruct1.paqdata.imaging.off2pLC;
    for iRew = 1:numel(off2p)
        iSpkGlobal = double(iFrames(off2p(iRew).block2p))+off2p(iRew).frames2p(end/2-2*imfs1+1:end/2+2*imfs1,1);
        matchData.rewT1.spk(:,:,iRew) = spkmat1(iSpkGlobal,:);
        matchData.rewT1.lick(:,iRew) = single(decimate(double(off2p(iRew).data(end/2-2*paqfs1+1:end/2+2*paqfs1,6)),50));
        x0 = double(off2p(iRew).data(end/2-2*paqfs1+1:end/2+2*paqfs1,10));
        x1 = smooth(x0,smspan/numel(x0),'lowess');
        v1 = diff(x1)*dataStruct1.sets.fs*-gaindiv1; v1(end+1) = v1(end);
        matchData.rewT1.mov(:,iRew) = single(decimate(v1(end/2-2*paqfs1+1:end/2+2*paqfs1),50));
    end
end

if isfield(trials2.spk,'spkLCmov')
    % If there is movement analysis data, use it - if there's one there's
    % all of them, so just stay in here for both movement and trial reward
    
    % Movement onset
    spkLocal = cat(3,trials2.spk.spkLCmov,trials2.spk.spkLUmov,trials2.spk.spkLOmov);
    behavLocal = cat(3,trials2.VR.VRLCmov,trials2.VR.VRLUmov,trials2.VR.VRLOmov);
    iTrialBad = isnan(squeeze(spkLocal(1,1,:)));
    spkLocal(:,:,iTrialBad) = []; behavLocal(:,:,iTrialBad) = [];
    matchData.mov2.spk = single(spkLocal(end/2-2*imfs2+1:end/2+2*imfs2,:,:));
    for iMov = 1:size(behavLocal,3)
        matchData.mov2.lick(:,iMov) = single(decimate(behavLocal(end/2-2*paqfs2+1:end/2+2*paqfs2,2,iMov),50));
        x0 = behavLocal(:,3,iMov);
        x1 = smooth(x0,smspan/numel(x0),'lowess');
        v1 = diff(x1)*dataStruct2.sets.fs*-gaindiv2; v1(end+1) = v1(end);
        matchData.mov2.mov(:,iMov) = single(decimate(v1(end/2-2*paqfs2+1:end/2+2*paqfs2),50));
    end
    
    % Trial reward
    spkLocal = trials2.spk.spkLC;
    behavLocal = trials2.VR.VR_FLC;
    iTrialBad = isnan(squeeze(spkLocal(1,1,:)));
    spkLocal(:,:,iTrialBad) = []; behavLocal(:,:,iTrialBad) = [];
    matchData.rewT2.spk = single(spkLocal(end/2-2*imfs2+1:end/2+2*imfs2,:,:));
    for iRew = 1:size(behavLocal,3)
        matchData.rewT2.lick(:,iRew) = single(decimate(behavLocal(end/2-2*paqfs2+1:end/2+2*paqfs2,2,iRew),50));
        x0 = behavLocal(:,3,iRew);
        x1 = smooth(x0,smspan/numel(x0),'lowess');
        v1 = diff(x1)*dataStruct2.sets.fs*-gaindiv2; v1(end+1) = v1(end);
        matchData.rewT2.mov(:,iRew) = single(decimate(v1(end/2-2*paqfs2+1:end/2+2*paqfs2),50));
    end
    
else
    % If there is no movement data, it's probably a v1 session and you need
    % to get the movement another way
    iFrames = [0;cumsum(sessionpair(2).fpblock)];
    
    % Use moveri to find trials with movement
    tr2p = behavStruct2.paqdata.imaging.trials2p;
    iMoveri = strcmpi(fieldnames(tr2p),'moveri');
    movTrials = struct2cell(tr2p); movTrials = squeeze(movTrials(iMoveri,:,:));
    tr2p = tr2p(~cellfun(@isempty,movTrials));
%     spkLocal = zeros(4*imfs2,matchData.nRoi2,numel(tr2p));
    behavLocal = zeros(4*paqfs2,2,numel(tr2p));
    for iMov = 1:numel(tr2p)
        iGlobal = tr2p(iMov).moveri(1)+tr2p(iMov).trialsession(1);
        iSpkLocal = tr2p(iMov).frames2p(find(tr2p(iMov).frames2p(:,2)>iGlobal,1),1);        
        iSpkGlobal = iFrames(tr2p(iMov).block2p)+iSpkLocal;
        matchData.mov2.spk(:,:,iMov) = single(spkmat2(iSpkGlobal-2*imfs2+1:iSpkGlobal+2*imfs2,:));
%         matchData.mov2.spk = single(spkLocal(:,:,iMov));
        behavLocal(:,:,iMov) = tr2p(iMov).data(tr2p(iMov).moveri(1)-2*paqfs2+1:tr2p(iMov).moveri(1)+2*paqfs2,[6,10]);
        matchData.mov2.lick(:,iMov) = single(decimate(behavLocal(:,1,iMov),50));
        x0 = behavLocal(:,2,iMov);
        x1 = smooth(x0,smspan/numel(x0),'lowess');
        v1 = diff(x1)*dataStruct2.sets.fs*-gaindiv2; v1(end+1) = v1(end);
        matchData.mov2.mov(:,iMov) = single(decimate(v1(end/2-2*paqfs2+1:end/2+2*paqfs2),50));
    end
    
    % Use off2pLC to find trials with reward
    off2p = behavStruct2.paqdata.imaging.off2pLC;
    for iRew = 1:numel(off2p)
        iSpkGlobal = double(iFrames(off2p(iRew).block2p))+off2p(iRew).frames2p(end/2-2*imfs2+1:end/2+2*imfs2,1);
        matchData.rewT2.spk(:,:,iRew) = spkmat2(iSpkGlobal,:);
        matchData.rewT2.lick(:,iRew) = single(decimate(double(off2p(iRew).data(end/2-2*paqfs2+1:end/2+2*paqfs2,6)),50));
        x0 = double(off2p(iRew).data(end/2-2*paqfs2+1:end/2+2*paqfs2,10));
        x1 = smooth(x0,smspan/numel(x0),'lowess');
        v1 = diff(x1)*dataStruct2.sets.fs*-gaindiv2; v1(end+1) = v1(end);
        matchData.rewT2.mov(:,iRew) = single(decimate(v1(end/2-2*paqfs2+1:end/2+2*paqfs2),50));
    end
end

% matchData.mov1.spk = single(rfs1.mover_trial.SPK);
% matchData.mov1.lick = single(rfs1.mover_trial.lick);
% matchData.mov1.mov = single(rfs1.mover_trial.mov*-gaindiv1);
% 
% matchData.mov2.spk = single(rfs2.mover_trial.SPK);
% matchData.mov2.lick = single(rfs2.mover_trial.lick);
% matchData.mov2.mov = single(rfs2.mover_trial.mov*-gaindiv2);
% 
% %% Take revelant data from structure - trial reward
% matchData.rewT1.spk = single(rfs1.trialrew.SPK);
% matchData.rewT1.lick = single(rfs1.trialrew.lick);
% matchData.rewT1.mov = single(rfs1.trialrew.mov*-gaindiv1);
% 
% matchData.rewT2.spk = single(rfs2.trialrew.SPK);
% matchData.rewT2.lick = single(rfs2.trialrew.lick);
% matchData.rewT2.mov = single(rfs2.trialrew.mov*-gaindiv2);

%% Take revelant data from structure - cued reward
matchData.rewC1.spk = single(rfs1.cuedrew.SPK);
matchData.rewC1.lick = single(rfs1.cuedrew.lick);
matchData.rewC1.mov = single(rfs1.cuedrew.mov*-gaindiv1);

matchData.rewC2.spk = single(rfs2.cuedrew.SPK);
matchData.rewC2.lick = single(rfs2.cuedrew.lick);
matchData.rewC2.mov = single(rfs2.cuedrew.mov*-gaindiv2);

%% Take revelant data from structure - rand reward
matchData.rewR1.spk = single(rfs1.randrew.SPK);
matchData.rewR1.lick = single(rfs1.randrew.lick);
matchData.rewR1.mov = single(rfs1.randrew.mov*-gaindiv2);

matchData.rewR2.spk = single(rfs2.randrew.SPK);
matchData.rewR2.lick = single(rfs2.randrew.lick);
matchData.rewR2.mov = single(rfs2.randrew.mov*-gaindiv2);

end
