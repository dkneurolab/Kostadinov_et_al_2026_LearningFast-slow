function fig2_GLMmain(inits)
%% Unpack settings from inits file
t0 = tic;
set(0,'DefaultFigureWindowStyle','docked');

% Common settings that may be used multiple times
savebool            = inits.savebool;
plotbool            = inits.plotbool;
spkhistbool         = inits.spkhistbool;
xvalfold            = inits.xvalfold;

% Basis function settings
BSsets.basisdur     = inits.basisdur;
BSsets.basiswid     = []; % set later inits.basiswid;
BSsets.basisdurSPK  = inits.basisdurSPK;
BSsets.nBasesSPK    = inits.nBasesSPK;
BSsets.nloffsetSPK  = inits.nloffsetSPK;


%% Get path locations and load stuff

if strcmp(inits.computer,'pc')
    % Get settings file in order to figure out version number
    dummypath           = dir(fullfile('E:\Analysis_current\04_Behav-processed',inits.name,inits.date,'Behav*'));
    behavpath           = fullfile(dummypath.folder,dummypath.name); clear dummypath
    behavroot           = fullfile('\\live.rd.ucl.ac.uk\ritd-ag-project-rd00g3-vkost96\Data\Data - DKrig',inits.date(1:4),'Behaviour',inits.date(1:7),inits.date(1:10));    
    cd(behavroot);
    setstr = ['*',num2str(inits.vrnum),'_settings.m'];
    setfile = dir(setstr); clear setstr;
    copyfile(setfile.name,fullfile(behavpath,'settings.m')); clear setfile;
    cd(behavpath);
    setsBehav = struct; setsBehav = settings(setsBehav); delete(fullfile(behavpath,'settings.m'));            
    
    % Get initial paths for everything
    impath              = fullfile('E:\Analysis_current\02_Imaging',inits.name,inits.date,inits.fov);
    imFile0             = dir(fullfile(impath,'*final.mat'));
    imFile              = fullfile(imFile0(end).folder,imFile0(end).name); clear imFile0;    
    alignpath           = fullfile('E:\Analysis_current\01_Behav+imaging',['Version_',num2str(setsBehav.version)],inits.name,[inits.date,'_',inits.fov]);
    imroot              = fullfile('\\live.rd.ucl.ac.uk\ritd-ag-project-rd00g3-vkost96\Data\Data - DKrig',inits.date(1:4),inits.name,inits.date,inits.fov);
    
    % Set target path
    targetpath          = fullfile(alignpath,'GLM'); %['GLM_',inits.cellgrping]);
    if ~exist(targetpath,'dir'); mkdir(targetpath); end
    
    if isfile(fullfile(targetpath,'fpblock.mat'))
        fpblock = load(fullfile(targetpath,'fpblock.mat'));
        fpblock = fpblock.fpblock;
        fprintf('Frames per block file found, moving on.\n\n');                
    else
        fprintf('Getting frames per block...');
        if ~isempty(inits.fov)
            cd(imroot)
            imrootFold = dir(imroot);
            % ignore root and up-level folders as well as other extraneous files
            imrootFold(ismember({imrootFold(:).name}, {'.', '..','.DS_Store','Thumbs.db'})) = [];
            i = cellfun(@(x) contains(x,'.tif'),{imrootFold(:).name});
            imrootFold2p = imrootFold(cat(1,imrootFold.isdir));
            imrootFold = imrootFold(i);
            fpblock = zeros(numel(imrootFold),1);
            for j = 1:numel(imrootFold)
                tempinfo = imfinfo(imrootFold(j).name);
                fpblock(j) = numel(tempinfo);
            end
            cd(imrootFold2p.name)
            imrootFold2pplane = dir('plane*');
            cd(imrootFold2pplane(1).name);
            Fall = load('Fall.mat');
            framediv = Fall.ops.nplanes*Fall.ops.nchannels;
            fpblock = fpblock/double(framediv);
            clear tempinfo imrootFold2p imrootFold2pplane framediv Fall
        else
            impath = fullfile('E:\Analysis_current\02_Imaging',inits.name,inits.date);
            imrootFold = dir(impath);
            % ignore root and up-level folders as well as other extraneous files
            imrootFold(ismember({imrootFold(:).name}, {'.', '..','.DS_Store','Thumbs.db'})) = [];
            impath0b = fullfile(imrootFold.folder,imrootFold.name);
            cd(impath0b)
            imrootFold = dir([inits.name,'_',inits.date,'*']);
            [~,i] = max(cat(1,imrootFold.datenum));
            imstruct = load(imrootFold(i).name);
            fpblock = zeros(size(imstruct.datdk.traceind,1),1);
            for j = 1:numel(fpblock)
                fpblock(j) = numel(imstruct.datdk.traceind{j,2});
            end
        end
        clear imrootFold i j impath0a impath0b imstruct
        save(fullfile(targetpath,'fpblock.mat'),'fpblock');
        fprintf('done.\n');
    end
    inits.fpblock = fpblock;
elseif strcmp(inits.computer,'mac')
    % Get settings file in order to figure out version number
    behavroot           = fullfile('/Volumes/ritd-ag-project-rd00g3-vkost96/Data/Data - DKrig',inits.date(1:4),'Behaviour',inits.date(1:7),inits.date(1:10));
%     cd(behavroot);
    setstr = ['*',num2str(inits.vrnum),'_settings.m'];
    setfile = dir(fullfile(behavroot,setstr)); clear setstr;
    copyfile(fullfile(setfile.folder,setfile.name),fullfile(setfile.folder,'settings.m')); clear setfile;
    setsBehav = struct; cd(behavroot);
    setsBehav = settings(setsBehav); 
    delete('settings.m');
    
    % Get initial paths for everything
    dummypath    = dir(fullfile('/Users/dimitar/Desktop/CF_learning_paper/04_Behav-processed',inits.name,inits.date,'Behav*'));
    behavpath    = fullfile(dummypath.folder,dummypath.name); clear dummypath
    alignpath    = fullfile('/Users/dimitar/Desktop/CF_learning_paper/01_Behav+imaging',['Version_',num2str(setsBehav.version)],inits.name,[inits.date,'_',inits.fov]);
    imroot      = fullfile('/Volumes/ritd-ag-project-rd00g3-vkost96/Data/Data - DKrig',inits.date(1:4),inits.name,inits.date,inits.fov);
    figFold     = fullfile(inits.dataFold,'paperFigs','Fig2','Sub1_GLM',sprintf('%s_%s_v%i_%s',inits.name,inits.fov,setsBehav.version,inits.date(end-4:end)));
    figFold2     = fullfile(inits.paperFold,'Figures','Fig2','Sub1_GLM',sprintf('%s_%s_v%i_%s',inits.name,inits.fov,setsBehav.version,inits.date(end-4:end)));
    
    % Set target path
    targetpath          = fullfile(alignpath,'GLM');
    if ~exist(targetpath,'dir'); mkdir(targetpath); end
    if ~exist(figFold,'dir'); mkdir(figFold); end
    if ~exist(figFold2,'dir'); mkdir(figFold2); end
    
    if isfile(fullfile(targetpath,'fpblock.mat'))
        fpblock = load(fullfile(targetpath,'fpblock.mat'))  ;
        fpblock = double(fpblock.fpblock);
        fprintf('Frames per block file found, moving on.\n\n');
    else
        fprintf('Getting frames per block...');
        impath = fullfile('/Users/dimitar/Desktop/CF_learning_paper/02_Imaging',inits.name,inits.date,inits.fov);
        Fall = load(fullfile(impath,'Fall.mat'));
        if ~isempty(inits.fov)
            cd(imroot)
            imrootFold = dir(imroot);
            % ignore root and up-level folders as well as other extraneous files
            imrootFold(ismember({imrootFold(:).name}, {'.', '..','.DS_Store','Thumbs.db'})) = [];
            i = cellfun(@(x) contains(x,'.tif'),{imrootFold(:).name});
            imrootFold2p = imrootFold(cat(1,imrootFold.isdir));
            imrootFold = imrootFold(i);
            fpblock = zeros(numel(imrootFold),1);
            for j = 1:numel(imrootFold)
                tempinfo = imfinfo(imrootFold(j).name);
                fpblock(j) = numel(tempinfo);
            end
            cd(imrootFold2p.name)
            imrootFold2pplane = dir('plane*');
            cd(imrootFold2pplane(1).name);
            Fall = load('Fall.mat');
            framediv = Fall.ops.nplanes*Fall.ops.nchannels;
            fpblock = fpblock/double(framediv);
            clear tempinfo imrootFold2p imrootFold2pplane framediv Fall
        else
            impath = fullfile('/Users/dimitar/Desktop/CF_learning_paper/02_Imaging',inits.name,inits.date);
            imrootFold = dir(impath);
            % ignore root and up-level folders as well as other extraneous files
            imrootFold(ismember({imrootFold(:).name}, {'.', '..','.DS_Store','Thumbs.db'})) = [];
            impath0b = fullfile(imrootFold.folder,imrootFold.name);
            cd(impath0b)
            imrootFold = dir([inits.name,'_',inits.date,'*']);
            [~,i] = max(cat(1,imrootFold.datenum));
            imstruct = load(imrootFold(i).name);
            fpblock = zeros(size(imstruct.datdk.traceind,1),1);
            for j = 1:numel(fpblock)
                fpblock(j) = numel(imstruct.datdk.traceind{j,2});
            end
            clear imrootFold i j impath0a impath0b imstruct
        end
        save(fullfile(targetpath,'fpblock.mat'),'fpblock');
        fprintf('done.\n');
        
    end
    inits.fpblock = fpblock;
end

%% Load data
% 1. Display recordign session
fprintf(['Mouse name: ',inits.name,'\n']);
fprintf(['Recording date: ',inits.date,'\n']);
fprintf(['VR version: ',num2str(setsBehav.version),', World: ',num2str(setsBehav.currentWorld),'\n\n']);

% 2. Extract relevant bits from behavstruct
fprintf('Loading behavioural data structure...');
behavstruct = load(fullfile(behavpath,[inits.name,'_',inits.date,'.mat']));
trials2p = behavstruct.datastruct.paqdata.imaging.trials2p; % behavioural data
rfs2p = behavstruct.datastruct.paqdata.imaging.RF; % behavioural data
paqdata = behavstruct.datastruct.dataraw;
clear behavstruct behavpath behavroot
fprintf('done.\n');

% 3. Extract spkmat - either grouped or ungrouped
fprintf('Loading fluorescence and spike matrix data structure...');
datdk = load(imFile);
dFF = datdk.datdk.finalparams.dff_30;
zF = zscore(dFF);
spkstruct = load(fullfile(alignpath,'spkmat.mat'));
if inits.grpspks
    spkmat = spkstruct.spkmat;
else
    spkmat = zeros(size(spkstruct.spkmat));
    for i = 1:numel(spkstruct.spki)
        locali = round(spkstruct.spki{i});
        for j = 1:numel(locali)
            spkmat(locali(j),i) = spkmat(locali(j),i)+1;
        end
    end
end
nCells = size(spkmat,2);
spksize = zeros(nCells,1);
for i = 1:nCells
    spksize(i) = mean(spkmat(spkmat(:,i) > 0,i));
    spkmat(:,i) = spkmat(:,i)/spksize(i);
    if inits.smoothspk
        x = -10:10;
        y = 1/sqrt(2*pi()*inits.smoothspksd^2)*exp(-x.^2/(2*inits.smoothspksd^2));
        y(1:10) = 0; y = y./sum(y);
        spkmat(:,i) = conv(spkmat(:,i),y,'same');
    end
end


clear spkstruct i j locali spksize
fprintf('done.\n');

% 4. Extract trialsonly and setstructure in order to get UCO designation, etc.
fprintf('Loading trial-sorted imaging data structure...');
trialstruct = load(fullfile(alignpath,'trialsonly.mat'));
setsImg = trialstruct.sets;
trialstruct = trialstruct.trialstructs;
setsBehav = Rot2p_changetrials(setsBehav);
if setsBehav.version == 1
    trialout = [];
else
    trialout.trialsLU0 = cat(1,trialstruct.LUstruct.trialind);
    trialout.trialsLC0 = cat(1,trialstruct.LCstruct.trialind);
    trialout.trialsLO0 = cat(1,trialstruct.LOstruct.trialind);
end
clear impath trialstruct
fprintf('done.\n\n');

%% Housekeeping stuff
% Define pre and post-trial durations:
setsLocal.trialbuffsec = (trials2p(1).inds(1) - trials2p(1).trialsession(1))/setsImg.fs;
setsLocal.imfs = setsImg.imfs;
setsLocal.paqfs = setsImg.fs;
setsLocal.prebuffsec = inits.prebuffsec;
setsLocal.postbuffsec = inits.postbuffsec;
setsLocal.prebuff = setsImg.fs*(setsLocal.trialbuffsec - inits.prebuffsec);
setsLocal.postbuff = setsImg.fs*(setsLocal.trialbuffsec - inits.postbuffsec);
setsLocal.version = setsBehav.version;
BSsets.basiswid = 1000/setsImg.imfs;
inits.imfs = setsImg.imfs;
inits.paqfs = setsImg.fs;

clear paqfs imfs prebuffsec postbuffsec

%% Convert the raw data into the experiment structure
% 1. Discrete trial vars
[discretevars] = aligndiscrete_dk(trials2p, setsLocal,setsBehav);

% 2. Continuous + calculate the movement onset timepoint & add into the discrete variable
[continvars, discretevars] = aligncont_dk(trials2p, discretevars, setsLocal);

% 3. Discrete reward vars
[cuedRdata0, randRdata0] = alignrewRFs_dk(rfs2p,paqdata,setsLocal, fpblock);

% 4. Combine discrete & continuous +  imaging frames in 1 block, note that
% cued and random reward data already has this structure
[trialdata0] = makeexptstruct_dk(trials2p, discretevars, continvars, fpblock);

% 5a. Add in spikes
% get average spikes for each microzone without the current cell
if strcmp(inits.cellgrping,'zones')
    [~, zoneinds, mzones] = microzones_dk(spkmat, setsImg);
elseif strcmp(inits.cellgrping,'cells')
    zoneinds = cell(nCells,1);
    mzones = zoneinds;
    for i = 1:nCells
        zoneinds{i} = i;
        mzones{i} = ['m',num2str(i)];
    end
end

% 5b. Add spking data (zones or single cells) to data strutures
[trialdata] = addinzonespikes_dk(trialdata0,mzones,zoneinds,spkmat,zF); % add in spiking data
[cuedRdata] = addinzonespikes_dk(cuedRdata0,mzones,zoneinds,spkmat,zF); % add in spiking data
[randRdata] = addinzonespikes_dk(randRdata0,mzones,zoneinds,spkmat,zF); % add in spiking data

[zoneFR0] = calczonefrs_dk(spkmat,zoneinds);
zoneFR = zoneFR0(:,1:numel(mzones));

% 6. Visualise the data so far (in imaging frames)
if plotbool
    exnums = 1:floor(numel(trials2p)/10):numel(trials2p); if numel(exnums) > 10; exnums = exnums(1:10); end
    visdataIMframes_dk(exnums,trialdata,mzones);
end

clear discretevars continvars spk zoneFR0 spkmat trialdata0 trials2p paqdata cuedRdata0 randRdata0 zoneInds i rfs2p 

%% Make basis fxns

% 1. Take 99th percentile of velocity trace as normalization value
velstd = prctile(abs(cat(1,trialdata.wheelvel)),99);

% 2. Make all the basis functions - this will build basis fxns for both cosine and boxcar basis functions
trialbasis = zonebasisfxns_dk(trialdata,velstd, BSsets,mzones,zoneFR,setsLocal,spkhistbool);
rewbasis = zonebasisfxns_dk_rewards(cuedRdata,randRdata,velstd, BSsets,mzones,zoneFR, setsLocal,spkhistbool);

% 3. Visualise basis functions for each covariate (maybe)
if plotbool
    visBSfunctions_dk(trialbasis.expt(1).frames, trialbasis.expt(1).trialon, 'ON', trialbasis.bs(1).TrialON, 'BS ON','discrete');
    visBSfunctions_dk(trialbasis.expt(1).frames, trialbasis.expt(1).wheelvelNL, 'velL', trialbasis.bs(1).BSwheelvelL, 'BS velL','cont');
    if spkhistbool
        visBSfunctions_dk(trialbasis.expt(1).frames, trialbasis.expt(1).m1, 'M1', trialbasis.bs(1).m1_spk, 'M1spike','cont');
    end
    visBSfunctions_dk(trialbasis.expt(1).frames, trialbasis.expt(1).trialon, 'ON', trialbasis.bsbox(1).TrialON, 'BS ON','discrete');
    visBSfunctions_dk(trialbasis.expt(1).frames, trialbasis.expt(1).reward, 'OFF', trialbasis.bsbox(1).RewTcorrect, 'BS reward','discrete');
    if spkhistbool
        visBSfunctions_dk(trialbasis.expt(1).frames, trialbasis.expt(1).m1, 'M1', trialbasis.bsbox(1).m1_spk, 'M1spike','cont');
    end
end
clear vtest velstd zoneFR clear trialdata velstd

%% Design matrix test: compile all the basis functions into 1

% 1. First run for cosines
[DM,DMdata] = makeDesignMatrix_dk(trialbasis, 1:numel(trialbasis.expt),'bs');
[DMrew,DMdatarew] = makeDesignMatrix_dk(rewbasis, 1:numel(rewbasis.expt),'bs');
% 2. Next run for boxcars
[DMbox,DMboxdata] = makeDesignMatrix_dk(trialbasis, 1:numel(trialbasis.expt),'bsbox');
[DMboxrew,DMboxdatarew] = makeDesignMatrix_dk(rewbasis, 1:numel(rewbasis.expt),'bsbox');

% Ytrialall = zeros(sum(cat(1,trialbasis.expt.duration)),numel(mzones));
% Yrewall = zeros(sum(cat(1,rewbasis.expt.duration)),numel(mzones));
% 
% for i = 1:numel(mzones)
%     Ytrialall(:,i) = getSpikes_dk(trialbasis,1:numel(trialbasis.expt),mzones{i});
%     Yrewall(:,i) = getSpikes_dk(rewbasis,1:numel(rewbasis.expt),mzones{i});
% end

Ytrialall = cat(1, trialbasis.expt.spkData);
Yrewall = cat(1, rewbasis.expt.spkData);

Ytrialall_zF = cat(1, trialbasis.expt.zFData);
Yrewall_zF = cat(1, rewbasis.expt.zFData);


% Visualise the design matrix (for the first 20 trials)
if plotbool
    %     endTrialIndices = cumsum([trialbasis.expt(:).duration]);
    X = DM(1:DMdata.endTrialIndices(20),:);
    mv = max(abs(X), [], 1); mv(isnan(mv)) = 1; % get max value of each column (ie. basis function) & change any nan values to 1
    X = bsxfun(@times, X, 1 ./ mv); % multiply the DM with 1/the max value of each basis function - normalise to scale?
    figure; imagesc(X'); % the y-axis being the variables & the x-axis time/trials
    % colormap('redblue')
    colormap gray;
    title('Cosine bases');
    Xbox = DMbox(1:DMboxdata.endTrialIndices(11),:);
    mvbox = max(abs(Xbox), [], 1); mvbox(isnan(mvbox)) = 1; % get max value of each column (ie. basis function) & change any nan values to 1
    Xbox = bsxfun(@times, Xbox, 1 ./ mvbox); % multiply the DM with 1/the max value of each basis function - normalise to scale?
    Xbox(isnan(Xbox)) = 0;
    figure;
    imagesc(Xbox',[-1 1]); % the y-axis being the variables & the x-axis time/trials
    colormap('redblue');
%     colormap gray;
    hold on;
    title('Boxcar bases');
    nBaseEnds = cumsum(DMboxdata.nBases);
    yticks((nBaseEnds+([0; nBaseEnds(1:end-1)]))/2);
    yticklabels(DMboxdata.fields);
    xticks((DMboxdata.endTrialIndices(2:11)+DMboxdata.endTrialIndices(1:10))/2);
    xticklabels(1:10);
    for i = 1:numel(nBaseEnds); yline(nBaseEnds(i)+0.5,'--','Color',[0 .4 .4],'LineWidth',1); end
    for i = 2:11; xline(DMboxdata.endTrialIndices(i)+0.5,'--k','LineWidth',1); end
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','in');
    xlabel('Trial number');
    ylabel('Task predictors');
end

if savebool
    save(fullfile(targetpath,sprintf('setfiles_%s.mat',inits.cellgrping)),'inits','BSsets','setsLocal','setsImg','mzones','trialout','setsBehav');
    save(fullfile(targetpath,sprintf('DM+spikes_%s.mat',inits.cellgrping)),'DM','DMdata','DMrew','DMdatarew','DMbox','DMboxdata','DMboxrew','DMboxdatarew','Ytrialall','Yrewall','Ytrialall_zF','Yrewall_zF');
    save(fullfile(targetpath,sprintf('basisstructs_%s.mat',inits.cellgrping)),'trialbasis','rewbasis');
end

%% Make GLM data structure
if setsBehav.version == 1
    [GLMdata, GLMholdout,DMdataBoxSmall] = make_GLMdata_v1_dk(DMdata, DMboxdata, inits, trialbasis, rewbasis, mzones,xvalfold,trialout);
else
    [GLMdata, GLMholdout,DMdataBoxSmall] = make_GLMdata_dk(DMdata, DMboxdata, inits, trialbasis, rewbasis, mzones,xvalfold,trialout);
end

if savebool
    save(fullfile(targetpath,sprintf('GLMdata_only_%s.mat',inits.cellgrping)),'GLMdata', 'GLMholdout','DMdataBoxSmall'); 
end

%% Get best fits using Ridge regression
%{
figfold = fullfile(targetpath,sprintf('SummaryFigs_%s',inits.cellgrping)); if ~exist(figfold,'dir'); mkdir(figfold); end

% mzones = fieldnames(GLMdata);
EV_cos = zeros(xvalfold,numel(mzones));
EV_box = zeros(xvalfold,numel(mzones));
figfold = fullfile(targetpath,'Summary figs'); if ~exist(figfold,'dir'); mkdir(figfold); end

if inits.cosbool
    bestBeta_cos = zeros(size(GLMdata(1).DMtest,2)+1,numel(mzones),xvalfold);
    bestPerf_cos = zeros(numel(mzones),xvalfold);
    bestLam_cos = zeros(numel(mzones),xvalfold);
    varExp_cos = zeros(numel(mzones),numel(lambda),xvalfold);
    for i = 1:xvalfold
        [bestBeta, bestPerf, bestLam, varExp] = runRidgeRegression_dk2(GLMdata(i),lambda,'cos');
        bestBeta_cos(:,:,i) = bestBeta;
        bestPerf_cos(:,i) = bestPerf;
        bestLam_cos(:,i) = bestLam;
        varExp_cos(:,:,i) = varExp;
    end
    EVmu_cos = mean(varExp_cos,3);
    [bestCosPerf, bestCosLi] = max(EVmu_cos, [], 2);
    bestBoxLam = lambda(bestBoxLi);
end


if inits.boxbool
    
    bestBeta_box_xval = zeros(size(GLMdata(1).DMtestbox,2)+1,numel(mzones),xvalfold);
    EVbox = zeros(numel(mzones),numel(lambda),xvalfold);
    for i = 1:xvalfold
        [~, ~, ~, varExp] = runRidgeRegression_dk2(GLMdata(i),lambda,'box');
        EVbox(:,:,i) = varExp;
    end
    EVbox_Lmu = mean(EVbox,3);
    [~, bestBoxLi] = max(EVbox_Lmu, [], 2);
    
    for i = 1:xvalfold
        bestBeta_box_xval(:,:,i) = runRidgeRegression_dk2best(GLMdata(i),lambda,'box',bestBoxLi);
    end
    
    EVbox_best = zeros(numel(mzones),xvalfold);
    for j = 1:numel(mzones)
        EVbox_best(j,:) = EVbox(j,bestBoxLi(j),:);
    end
    EVbox_bestmu = mean(EVbox_best,2);
    
    fw_box = figure;
    nBases = DMdataBoxSmall.nBases;
    cumnBases = [1; 1+cumsum(nBases)];
    basenames = cat(1,{'Offset'},DMdataBoxSmall.fields);
    
    for i = 1:numel(mzones)
        wboxtemp = squeeze(bestBeta_box_xval(:,i,:));
        xwbox = (0.5:1:size(wboxtemp,1))';
        wboxmean = nanmean(wboxtemp,2);
        wboxsd = nanstd(wboxtemp,[],2);
        wboxylim = max(sum([abs(wboxmean) abs(wboxsd)],2))*1.1;
        if isnan(wboxylim); wboxylim = 0.5; end;
        if strcmpi(inits.cellgrping,'zones')
            subplot(numel(mzones),1,i);
        elseif strcmpi(inits.cellgrping,'cells')
            subplot(ceil(numel(mzones)/5),5,i);
        end 
        bar(xwbox,wboxmean,'EdgeColor','none'); hold on;
        errorbar(xwbox,wboxmean,wboxsd,'k','CapSize',0,'Marker','none','LineStyle','none');
        ylim([-wboxylim wboxylim]);
        xlim([0 numel(wboxmean)]);
        box off
        plot([cumnBases'; cumnBases'],[-wboxylim*ones(1,numel(cumnBases)); wboxylim*ones(1,numel(cumnBases))],'k-')
        xticks(cumnBases);
        xticklabels({});
        if mod(i,6) == 1; title(['Microzone number ',mzones{i}(2:end)]); end
        if i == numel(mzones)
%             suptitle('Summary of microzone beta values, boxcar filters')
            namelocs = mean(cat(2,[0; cumnBases(1:end-1)], cumnBases),2);
            text(namelocs,repmat(-1.5*wboxylim,size(namelocs)),basenames,'HorizontalAlignment','center');
        end
    end
    
    if savebool
        savefig(fw_box,fullfile(figfold,'Betas_box.fig'),'compact');
        print(fw_box,fullfile(figfold,'Betas_box.eps'), '-depsc', '-painters');
        saveas(fw_box,fullfile(figfold,'Betas_box.png'));
        save(fullfile(targetpath,sprintf('GLMoptims_%s.mat',inits.cellgrping)), 'bestBeta_box_xval', 'bestBoxLi', 'lambda', 'EVbox', 'EVbox_best', 'EVbox_Lmu','EVbox_bestmu');
    end
    
end

clear varExp

%}
%% Calculate trial-averaged activity and fits
%{
figfold2 = fullfile(targetpath,sprintf('MeanFigs_%s',inits.cellgrping)); if ~exist(figfold2,'dir'); mkdir(figfold2); end
% figfold3 = fullfile(targetpath,'Mean figs2'); if ~exist(figfold3,'dir'); mkdir(figfold3); end
inits.restptime = setsBehav.responsetime;

meanstruct = struct;
meanstruct.fullData = plotHoldout(trialbasis,rewbasis,mzones,lambda(bestBoxLi),...
    GLMholdout,DMdataBoxSmall,inits,figfold2,savebool);

baseBlocks = false(numel(DMdataBoxSmall.iPred),numel(DMdataBoxSmall.nBases));
nBaseEnds = [0; cumsum(DMdataBoxSmall.nBases)];

for i = 1:numel(DMdataBoxSmall.nBases)
    baseBlocks(nBaseEnds(i)+1:nBaseEnds(i+1),i) = true;
    meanstruct.(['no',DMdataBoxSmall.fields{i}]) = plotHoldout_LeaveOut(trialbasis,rewbasis,...
        mzones,lambda(bestBoxLi),GLMholdout,DMdataBoxSmall,inits,baseBlocks(:,i),i,figfold2,savebool);
end

% meanstruct2.fullData = meanstruct.fullData;
% for i = 1:numel(DMdataBoxSmall.nBases)
%     baseBlocks(nBaseEnds(i)+1:nBaseEnds(i+1),i) = true;
%     meanstruct2.(['no',DMdataBoxSmall.fields{i}]) = plotHoldout_LeaveOut_v2(trialbasis,rewbasis,...
%         mzones,lambda(bestBoxLi),GLMholdout,DMdataBoxSmall,inits,baseBlocks(:,i),i,figfold3,savebool);
% end

if savebool
    save(fullfile(targetpath,sprintf('meanstruct_%s.mat',inits.cellgrping)),'meanstruct');
    %     save(fullfile(targetpath,'meanstruct2.mat'),'meanstruct2');
end

cd(targetpath);

%}
toc(t0);

end