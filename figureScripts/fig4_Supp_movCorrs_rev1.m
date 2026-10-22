function fig4_Supp_movCorrs_rev1(matchSpont, glmParams, simThresh, nMin, xMin, paperFold, savebool)
% Loop through each v1-v4n match and identify similar movement trajectories
% Make heatmaps of spk responses in those FOVs to compare directly

%% Loading stuff
% Load in v1-4n matched data
v1v4nMatch = load(fullfile(paperFold,'Figures','Fig4','Sub1_v1v4n','match_v1_v4n.mat'));
v1v4nMatch = v1v4nMatch.matchV14;
% Load in v41-4n matched data
v41v4nMatch = load(fullfile(paperFold,'Figures','Fig4','Sub1_v41v4n','match_v4_v4.mat'));
v41v4nMatch = v41v4nMatch.matchV14;

% Load in GLM info Load in v4 sigSummary
glmFold0 = fullfile(paperFold,'Figures','Fig3','Sub3_SummaryBars');
glmFold = fullfile(glmFold0,sprintf('%s_%s_alpha%i',glmParams.dataType,glmParams.distr,glmParams.alpha*100));
glmSummary = load(fullfile(glmFold,'sigSummary_v4n.mat')); glmSummary = glmSummary.glmSummary;



% Time settings
imfs = 30;
paqfs = 100;
t2p = (1:120)'/imfs - 2 - 1/imfs/2;
tpaq = (1:400)'/paqfs - 2 - 1/paqfs/2;

%% Work through v1-v4 pairs
% Use similarity threshold and make sure that there are a minimum of nMin movements
v1v4nMov = struct;
v41v4nMov = struct;
for ifov = 1:numel(v1v4nMatch)
    % Find the matchSpont row for this FOV, and get the indices for the matched ROIs
    delDays = 0;
    for jSpont = 1:numel(matchSpont)
        if strcmpi(v1v4nMatch(ifov).name, matchSpont(jSpont).name) & matchSpont(jSpont).v1 == 1 && matchSpont(jSpont).v2 == 4 && ...
                ((v1v4nMatch(ifov).ilobv == true & strcmpi(matchSpont(jSpont).fov, 'lobv')) || (v1v4nMatch(ifov).isim2 == true & strcmpi(matchSpont(jSpont).fov, 'sim2')))
            delDaysLocal = datetime(matchSpont(jSpont).date2)-datetime(matchSpont(jSpont).date1);            
            if delDaysLocal > delDays % Take longest interval between v1 and v4 as the v1 vs v4n session
                iroi_v1v4 = matchSpont(jSpont).iROIs2;
                delDays = delDaysLocal;
            end
        elseif strcmpi(v41v4nMatch(ifov).name, matchSpont(jSpont).name) & matchSpont(jSpont).v1 == 4 && matchSpont(jSpont).v2 == 4 && ...
                ((v41v4nMatch(ifov).ilobv == true & strcmpi(matchSpont(jSpont).fov, 'lobv')) || (v41v4nMatch(ifov).isim2 == true & strcmpi(matchSpont(jSpont).fov, 'sim2')))
            if datetime(matchSpont(jSpont).date2) > datetime(matchSpont(jSpont).date1) % Make sure that you get the correct indices (later session == v4n, earlier == v41
                iroi_v41v4n = matchSpont(jSpont).iROIs2;
            else
                iroi_v41v4n = matchSpont(jSpont).iROIs1;
            end            
        end
    end
    % Get indices of sig mov rois:
    for jGlm = 1:numel(glmSummary)
        if strcmpi(v1v4nMatch(ifov).name, glmSummary(jGlm).name) && v1v4nMatch(ifov).ilobv == glmSummary(jGlm).iLobV && v1v4nMatch(ifov).isim2 == glmSummary(jGlm).iSim2
            sigBoolMov = false(size(glmSummary(jGlm).sigBool));
            sigBoolMov(1:numel(glmSummary(jGlm).sigBoolMov)) = glmSummary(jGlm).sigBoolMov;
        end
    end

    % V1-V4 stuff first
    seshlocal = v1v4nMatch(ifov).mov;
    xMov = cumsum(seshlocal(1).mov)/100;
    ixMov = find(max(xMov) > xMin); % only movements > xMin (e.g. 0.5 cm) displacement to be taken
    iv1 = [];
    iv4 = [];
    if sum(ixMov) >= nMin
        movCorrs = corr(seshlocal(1).mov(:,ixMov),seshlocal(2).mov);
        for jmov = 1:size(movCorrs,1)
            if any(movCorrs(jmov,:) >= simThresh)
                [~,imax] = max(movCorrs(jmov,:));
                iv1 = [iv1; ixMov(jmov)]; %#ok<AGROW>
                iv4 = [iv4; imax]; %#ok<AGROW>
            end
        end
        seshlocal(1).iMatched = iv1;
        seshlocal(2).iMatched = iv4;
    end
    v1v4nMov(ifov).name = v1v4nMatch(ifov).name;
    v1v4nMov(ifov).ilobv = v1v4nMatch(ifov).ilobv;
    v1v4nMov(ifov).isim2 = v1v4nMatch(ifov).isim2;
    v1v4nMov(ifov).iMatched_v1 = iv1;
    v1v4nMov(ifov).spk_v1 = seshlocal(1).spk(:,:,iv1);
    v1v4nMov(ifov).lick_v1 = seshlocal(1).lick(:,iv1);
    v1v4nMov(ifov).mov_v1 = seshlocal(1).mov(:,iv1);
    v1v4nMov(ifov).iMatched_v4 = iv4;
    v1v4nMov(ifov).spk_v4 = seshlocal(2).spk(:,:,iv4);
    v1v4nMov(ifov).lick_v4 = seshlocal(2).lick(:,iv4);
    v1v4nMov(ifov).mov_v4 = seshlocal(2).mov(:,iv4);

    v1v4nMov(ifov).sigmov_v4 = sigBoolMov(iroi_v1v4);

    if ~isempty(iv1)
        % Calculate correlations between rois, as well as nearest
        % neighbour, and random correlation pair
        rng(2349);
        mov_v1 = mean(v1v4nMov(ifov).spk_v1(end/2-imfs+1:end/2+imfs,:,:),3,'omitnan');
        mov_v4 = mean(v1v4nMov(ifov).spk_v4(end/2-imfs+1:end/2+imfs,:,:),3,'omitnan');
        movcorr_1_4 = corr(mov_v1, mov_v4);
        nrois = size(mov_v1,2);
        iDiag = (1:nrois+1:nrois^2)';
        iShift = iDiag+1; iShift(end) = 1;
        iRand = randperm(nrois^2,nrois)';
        iRandMatch = find(iDiag == iRand);
        iRand(iRandMatch) = iRand(iRandMatch)+1; iRand(iRand > nrois^2) = 1;
        v1v4nMov(ifov).movcorrs = [movcorr_1_4(iDiag) movcorr_1_4(iShift) movcorr_1_4(iRand)];

        % Significantly modulated v1 cells
        spkbase_1 = squeeze(mean(v1v4nMov(ifov).spk_v1(1:imfs,:,:),1,'omitnan'));
        spkmod_1 = squeeze(mean(v1v4nMov(ifov).spk_v1(end/2-0.3*imfs+1:end/2,:,:),1,'omitnan'));
        if size(spkbase_1,1) == 1 & size(spkbase_1,2) > 1
            spkbase_1 = spkbase_1';
            spkmod_1 = spkmod_1';
        end
        sig_v1 = false(size(spkbase_1,1),1);
        for k = 1:size(spkbase_1,1)
            [~,sig_v1(k)] = signrank(spkbase_1(k,:),spkmod_1(k,:));
        end
        v1v4nMov(ifov).sig_v1 = sig_v1;

        % Significantly modulated v4n cells
        spkbase_4 = squeeze(mean(v1v4nMov(ifov).spk_v4(1:imfs,:,:),1,'omitnan'));
        spkmod_4 = squeeze(mean(v1v4nMov(ifov).spk_v4(end/2-0.3*imfs+1:end/2,:,:),1,'omitnan'));
        if size(spkbase_4,1) == 1 & size(spkbase_4,2) > 1
            spkbase_4 = spkbase_4';
            spkmod_4 = spkmod_4';
        end
        sig_v4 = false(size(spkbase_4,1),1);
        for k = 1:size(spkbase_4,1)
            [~,sig_v4(k)] = signrank(spkbase_4(k,:),spkmod_4(k,:));
        end
        v1v4nMov(ifov).sig_v4 = sig_v4;
        v1v4nMov(ifov).nMatch = numel(iv1);

    else
        v1v4nMov(ifov).movcorrs = zeros(0,3);
        v1v4nMov(ifov).sig_v1 = [];
        v1v4nMov(ifov).sig_v4 = [];
        v1v4nMov(ifov).nMatch = 0;
    end


    % V41-V4 stuff second
    seshlocal = v41v4nMatch(ifov).mov;
    xMov = cumsum(seshlocal(1).mov)/100;
    ixMov = find(max(xMov) > xMin); % only movements > xMin (e.g. 0.5 cm) displacement to be taken
    iv41 = [];
    iv4 = [];
    if sum(ixMov) >= nMin
        movCorrs = corr(seshlocal(1).mov(:,ixMov),seshlocal(2).mov);
        for jmov = 1:size(movCorrs,1)
            if any(movCorrs(jmov,:) >= simThresh)
                [~,imax] = max(movCorrs(jmov,:));
                iv41 = [iv41; ixMov(jmov)]; %#ok<AGROW>
                iv4 = [iv4; imax]; %#ok<AGROW>
            end
        end
        seshlocal(1).iMatched = iv41;
        seshlocal(2).iMatched = iv4;
    end
    v41v4nMov(ifov).name = v41v4nMatch(ifov).name;
    v41v4nMov(ifov).ilobv = v41v4nMatch(ifov).ilobv;
    v41v4nMov(ifov).isim2 = v41v4nMatch(ifov).isim2;
    v41v4nMov(ifov).iMatched_v41 = iv41;
    v41v4nMov(ifov).spk_v41 = seshlocal(1).spk(:,:,iv41);
    v41v4nMov(ifov).lick_v41 = seshlocal(1).lick(:,iv41);
    v41v4nMov(ifov).mov_v41 = seshlocal(1).mov(:,iv41);
    v41v4nMov(ifov).iMatched_v4n = iv4;
    v41v4nMov(ifov).spk_v4n = seshlocal(2).spk(:,:,iv4);
    v41v4nMov(ifov).lick_v4n = seshlocal(2).lick(:,iv4);
    v41v4nMov(ifov).mov_v4n = seshlocal(2).mov(:,iv4);
    v41v4nMov(ifov).sigmov_v4 = sigBoolMov(iroi_v41v4n);

    if ~isempty(iv41)
        % Calculate correlations between rois, as well as nearest
        % neighbour, and random correlation pair
        rng(2349);
        mov_v41 = mean(v41v4nMov(ifov).spk_v41(end/2-imfs+1:end/2+imfs,:,:),3,'omitnan');
        mov_v4n = mean(v41v4nMov(ifov).spk_v4n(end/2-imfs+1:end/2+imfs,:,:),3,'omitnan');
        movcorr_41_4n = corr(mov_v41, mov_v4n);
        nrois = size(mov_v41,2);
        iDiag = (1:nrois+1:nrois^2)';
        iShift = iDiag+1; iShift(end) = 1;
        iRand = randperm(nrois^2,nrois)';
        iRandMatch = find(iDiag == iRand);
        iRand(iRandMatch) = iRand(iRandMatch)+1; iRand(iRand > nrois^2) = 1;
        v41v4nMov(ifov).movcorrs = [movcorr_41_4n(iDiag) movcorr_41_4n(iShift) movcorr_41_4n(iRand)];

        % Significantly modulated v41 cells
        spkbase_41 = squeeze(mean(v41v4nMov(ifov).spk_v41(1:imfs,:,:),1,'omitnan'));
        spkmod_41 = squeeze(mean(v41v4nMov(ifov).spk_v41(end/2-0.3*imfs+1:end/2,:,:),1,'omitnan'));
        if size(spkbase_41,1) == 1 & size(spkbase_41,2) > 1
            spkbase_41 = spkbase_41';
            spkmod_41 = spkmod_41';
        end
        sig_v41 = false(size(spkbase_41,1),1);
        for k = 1:size(spkbase_41,1)
            [~,sig_v41(k)] = signrank(spkbase_41(k,:),spkmod_41(k,:));
        end
        v41v4nMov(ifov).sig_v41 = sig_v41;

        % Significantly modulated v4n cells
        spkbase_4n = squeeze(mean(v41v4nMov(ifov).spk_v4n(1:imfs,:,:),1,'omitnan'));
        spkmod_4n = squeeze(mean(v41v4nMov(ifov).spk_v4n(end/2-0.3*imfs+1:end/2,:,:),1,'omitnan'));
        if size(spkbase_4n,1) == 1 & size(spkbase_4n,2) > 1
            spkbase_4n = spkbase_4n';
            spkmod_4n = spkmod_4n';
        end
        sig_v4n = false(size(spkbase_4n,1),1);
        for k = 1:size(spkbase_4n,1)
            [~,sig_v4n(k)] = signrank(spkbase_4n(k,:),spkmod_4n(k,:));
        end
        v41v4nMov(ifov).sig_v4n = sig_v4n;
        v41v4nMov(ifov).nMatch = numel(iv41);
    else
        v41v4nMov(ifov).movcorrs = zeros(0,3);
        v41v4nMov(ifov).sig_v41 = [];
        v41v4nMov(ifov).sig_v4n = [];
        v41v4nMov(ifov).nMatch = 0;

    end

    f1 = figure;
    if v1v4nMov(ifov).ilobv; fov = 'lobv'; else; fov = 'sim2'; end
    sgtitle(sprintf('%s: %s matched movement aligned activity',v1v4nMov(ifov).name, fov),'FontSize',24);
    % v1-v4: v1 stuff
    spkdata = v1v4nMov(ifov).spk_v1;
    movdata = v1v4nMov(ifov).mov_v1;
    lickdata = v1v4nMov(ifov).lick_v1;
    sigBoolMov = v1v4nMov(ifov).sigmov_v4;
    firstLabel = false(4,1);
    if ~isempty(spkdata)
        subplot(4,4,1);
        % nroi = size(spkdata,2);
        nroi = sum(sigBoolMov);
        imagesc(t2p,[1,nroi],mean(spkdata(:,sigBoolMov,:),3)',[0 10]); xlim([-1 1]);
        xline(0,'--w', 'LineWidth',1)
        colormap viridis; box off;
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
        if ~firstLabel(1)
            firstLabel(1) = true;
            ylabel('PC dendrite number')
        end
        title(sprintf('v1-v4n: v1 (%i)',nroi))
        subplot(4,4,5);
        % nroi = size(spkdata,2);
        nroi2 = sum(~sigBoolMov);
        imagesc(t2p,[1,nroi2],mean(spkdata(:,~sigBoolMov,:),3)',[0 10]); xlim([-1 1]);
        xline(0,'--w', 'LineWidth',1)
        colormap viridis; box off;
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
        title(sprintf('(%i)',nroi2))

        % Mean recording
        subplot(6,4,13);
        h = boundedline(t2p, mean(mean(spkdata(:,sigBoolMov,:),3),2), std(mean(spkdata(:,sigBoolMov,:),3),[],2),'-k','alpha'); set(h, 'linewidth', 1.5)
        h2 = boundedline(t2p, mean(mean(spkdata(:,~sigBoolMov,:),3),2), std(mean(spkdata(:,~sigBoolMov,:),3),[],2),'-b','alpha'); set(h2, 'linewidth', 1.5)
        xline(0,'--k', 'LineWidth',1)
        axis([-1 1 -2 8]); yticks([-2:2:8])
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
        if ~firstLabel(2)
            firstLabel(2) = true;
            ylabel ('Event rate (spk/s)')
        end
        % Licking
        subplot(6,4,17);
        h = boundedline(tpaq,mean(lickdata,2),std(lickdata,[],2),'-k','alpha'); set(h, 'linewidth', 1.5)
        xline(0,'--k', 'LineWidth',1)
        axis([-1 1 -.1 1.1]); yticks([0:.25:1])
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
        if ~firstLabel(3)
            firstLabel(3) = true;
            ylabel ('Licking prob.')
        end
        % Movement
        subplot(6,4,21);
        h = boundedline(tpaq,mean(movdata,2),std(movdata,[],2),'-k','alpha'); set(h, 'linewidth', 1.5)
        xline(0,'--k', 'LineWidth',1)
        axis([-1 1 -2.5 15]); yticks([-2.5:2.5:15])
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
        if ~firstLabel(4)
            firstLabel(4) = true;
            ylabel ({'Wheel movement','(cm/s)'})
        end
    end

    % v1-v4: v4 stuff
    spkdata = v1v4nMov(ifov).spk_v4;
    movdata = v1v4nMov(ifov).mov_v4;
    lickdata = v1v4nMov(ifov).lick_v4;
    if ~isempty(spkdata)
        subplot(4,4,2);
        % nroi = size(spkdata,2);
        nroi = sum(sigBoolMov);
        imagesc(t2p,[1,nroi],mean(spkdata(:,sigBoolMov,:),3)',[0 10]); xlim([-1 1]);
        xline(0,'--w', 'LineWidth',1)
        colormap viridis; box off;
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
        if ~firstLabel(1)
            firstLabel(1) = true;
            ylabel('PC dendrite number')
        end
        title(sprintf('v1-v4n: v4 (%i)',nroi))
        subplot(4,4,6);
        % nroi = size(spkdata,2);
        nroi2 = sum(~sigBoolMov);
        imagesc(t2p,[1,nroi2],mean(spkdata(:,~sigBoolMov,:),3)',[0 10]); xlim([-1 1]);
        xline(0,'--w', 'LineWidth',1)
        colormap viridis; box off;
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
        title(sprintf('(%i)',nroi2))

        % Mean recording
        subplot(6,4,14);
        h = boundedline(t2p, mean(mean(spkdata(:,sigBoolMov,:),3),2), std(mean(spkdata(:,sigBoolMov,:),3),[],2),'-k','alpha'); set(h, 'linewidth', 1.5)
        h2 = boundedline(t2p, mean(mean(spkdata(:,~sigBoolMov,:),3),2), std(mean(spkdata(:,~sigBoolMov,:),3),[],2),'-b','alpha'); set(h2, 'linewidth', 1.5)
        xline(0,'--k', 'LineWidth',1)
        axis([-1 1 -2 8]); yticks([-2:2:8])
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
        if ~firstLabel(2)
            firstLabel(2) = true;
            ylabel ('Event rate (spk/s)')
        end
        % Licking
        subplot(6,4,18);
        h = boundedline(tpaq,mean(lickdata,2),std(lickdata,[],2),'-k','alpha'); set(h, 'linewidth', 1.5)
        xline(0,'--k', 'LineWidth',1)
        axis([-1 1 -.1 1.1]); yticks([0:.25:1])
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
        if ~firstLabel(3)
            firstLabel(3) = true;
            ylabel ('Licking prob.')
        end
        % Movement
        subplot(6,4,22);
        h = boundedline(tpaq,mean(movdata,2),std(movdata,[],2),'-k','alpha'); set(h, 'linewidth', 1.5)
        xline(0,'--k', 'LineWidth',1)
        axis([-1 1 -2.5 15]); yticks([-2.5:2.5:15])
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
        if ~firstLabel(4)
            firstLabel(4) = true;
            ylabel ({'Wheel movement','(cm/s)'})
        end
    end

    % v41-v4: v41 stuff
    spkdata = v41v4nMov(ifov).spk_v41;
    movdata = v41v4nMov(ifov).mov_v41;
    lickdata = v41v4nMov(ifov).lick_v41;
    sigBoolMov = v41v4nMov(ifov).sigmov_v4;
    if ~isempty(spkdata)
        subplot(4,4,3);
        % nroi = size(spkdata,2);
        nroi = sum(sigBoolMov);
        imagesc(t2p,[1,nroi],mean(spkdata(:,sigBoolMov,:),3)',[0 10]); xlim([-1 1]);
        xline(0,'--w', 'LineWidth',1)
        colormap viridis; box off;
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
        if ~firstLabel(1)
            firstLabel(1) = true;
            ylabel('PC dendrite number')
        end
        title(sprintf('v41-v4n: v41 (%i)',nroi))
        subplot(4,4,7);
        % nroi = size(spkdata,2);
        nroi2 = sum(~sigBoolMov);
        imagesc(t2p,[1,nroi2],mean(spkdata(:,~sigBoolMov,:),3)',[0 10]); xlim([-1 1]);
        xline(0,'--w', 'LineWidth',1)
        colormap viridis; box off;
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
        title(sprintf('(%i)',nroi2))

        % Mean recording
        subplot(6,4,15);
        h = boundedline(t2p, mean(mean(spkdata(:,sigBoolMov,:),3),2), std(mean(spkdata(:,sigBoolMov,:),3),[],2),'-k','alpha'); set(h, 'linewidth', 1.5)
        h2 = boundedline(t2p, mean(mean(spkdata(:,~sigBoolMov,:),3),2), std(mean(spkdata(:,~sigBoolMov,:),3),[],2),'-b','alpha'); set(h2, 'linewidth', 1.5)
        xline(0,'--k', 'LineWidth',1)
        axis([-1 1 -2 8]); yticks([-2:2:8])
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
        if ~firstLabel(2)
            firstLabel(2) = true;
            ylabel ('Event rate (spk/s)')
        end
        % Licking
        subplot(6,4,19);
        h = boundedline(tpaq,mean(lickdata,2),std(lickdata,[],2),'-k','alpha'); set(h, 'linewidth', 1.5)
        xline(0,'--k', 'LineWidth',1)
        axis([-1 1 -.1 1.1]); yticks([0:.25:1])
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
        if ~firstLabel(3)
            firstLabel(3) = true;
            ylabel ('Licking prob.')
        end
        % Movement
        subplot(6,4,23);
        h = boundedline(tpaq,mean(movdata,2),std(movdata,[],2),'-k','alpha'); set(h, 'linewidth', 1.5)
        xline(0,'--k', 'LineWidth',1)
        axis([-1 1 -2.5 15]); yticks([-2.5:2.5:15])
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
        if ~firstLabel(4)
            firstLabel(4) = true;
            ylabel ({'Wheel movement','(cm/s)'})
        end
    end

    % v41-v4n: v4n stuff
    spkdata = v41v4nMov(ifov).spk_v4n;
    movdata = v41v4nMov(ifov).mov_v4n;
    lickdata = v41v4nMov(ifov).lick_v4n;    
    if ~isempty(spkdata)
        subplot(4,4,4);
        % nroi = size(spkdata,2);
        nroi = sum(sigBoolMov);
        imagesc(t2p,[1,nroi],mean(spkdata(:,sigBoolMov,:),3)',[0 10]); xlim([-1 1]);
        xline(0,'--w', 'LineWidth',1)
        colormap viridis; box off;
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
        if ~firstLabel(1)
            firstLabel(1) = true;
            ylabel('PC dendrite number')
        end
        title(sprintf('v41-v4n: v4n (%i)',nroi))
        subplot(4,4,8);
        % nroi = size(spkdata,2);
        nroi2 = sum(~sigBoolMov);
        imagesc(t2p,[1,nroi2],mean(spkdata(:,~sigBoolMov,:),3)',[0 10]); xlim([-1 1]);
        xline(0,'--w', 'LineWidth',1)
        colormap viridis; box off;
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
        title(sprintf('(%i)',nroi2))

        % Mean recording
        subplot(6,4,16);
        h = boundedline(t2p, mean(mean(spkdata(:,sigBoolMov,:),3),2), std(mean(spkdata(:,sigBoolMov,:),3),[],2),'-k','alpha'); set(h, 'linewidth', 1.5)
        h2 = boundedline(t2p, mean(mean(spkdata(:,~sigBoolMov,:),3),2), std(mean(spkdata(:,~sigBoolMov,:),3),[],2),'-b','alpha'); set(h2, 'linewidth', 1.5)
        xline(0,'--k', 'LineWidth',1)
        axis([-1 1 -2 8]); yticks([-2:2:8])
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
        if ~firstLabel(2)
            firstLabel(2) = true;
            ylabel ('Event rate (spk/s)')
        end
        % Licking
        subplot(6,4,20);
        h = boundedline(tpaq,mean(lickdata,2),std(lickdata,[],2),'-k','alpha'); set(h, 'linewidth', 1.5)
        xline(0,'--k', 'LineWidth',1)
        axis([-1 1 -.1 1.1]); yticks([0:.25:1])
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
        if ~firstLabel(3)
            firstLabel(3) = true;
            ylabel ('Licking prob.')
        end
        % Movement
        subplot(6,4,24);
        h = boundedline(tpaq,mean(movdata,2),std(movdata,[],2),'-k','alpha'); set(h, 'linewidth', 1.5)
        xline(0,'--k', 'LineWidth',1)
        axis([-1 1 -2.5 15]); yticks([-2.5:2.5:15])
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
        if ~firstLabel(4)
            firstLabel(4) = true;
            ylabel ({'Wheel movement','(cm/s)'})
        end
    end

    if savebool
        figFold = fullfile(paperFold,'Figures','Fig4','Sub1_v1v41v4n_movCorrs_rev1');
        if ~exist(figFold,'dir'); mkdir(figFold); end
        figname = sprintf('%s_%s_movCorrMatch',v1v4nMov(ifov).name, fov);
        savefig(f1, fullfile(figFold,[figname,'.fig']));
        saveas(f1,fullfile(figFold,[figname,'.png']));
        print(f1,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');

    end

    clear firstLabel fov iDiag iRand iRandMatch iShift iv1 iv4 iv41 ixMov k spkdata movdata lickdata ...
        nroi h figname
end


%% Summary across mice

spkv1v4_1 = [];
movv1v4_1 = [];
lickv1v4_1 = [];
spkv1v4_4 = [];
movv1v4_4 = [];
lickv1v4_4 = [];
movcorrs_v1v4 = [];
sig_v1v4 = [];

for ifov = 1:numel(v1v4nMov)
    if v1v4nMov(ifov).nMatch >= nMin
        spkv1v4_1 = cat(2, spkv1v4_1, mean(v1v4nMov(ifov).spk_v1,3,'omitnan')); % [spkv1v4_1 mean(v1v4nMov(ifov).spk_v1,3,'omitnan')];
        movv1v4_1 = cat(2, movv1v4_1, mean(v1v4nMov(ifov).mov_v1,2,'omitnan')); % [movv1v4_1 mean(v1v4nMov(ifov).mov_v1,2,'omitnan')];
        lickv1v4_1 = cat(2, lickv1v4_1, mean(v1v4nMov(ifov).lick_v1,2,'omitnan')); % [lickv1v4_1 mean(v1v4nMov(ifov).lick_v1,2,'omitnan')];

        spkv1v4_4 = cat(2, spkv1v4_4, mean(v1v4nMov(ifov).spk_v4,3,'omitnan')); % [spkv1v4_4 mean(v1v4nMov(ifov).spk_v4,3,'omitnan')];
        movv1v4_4 = cat(2, movv1v4_4, mean(v1v4nMov(ifov).mov_v4,2,'omitnan')); % [movv1v4_4 mean(v1v4nMov(ifov).mov_v4,2,'omitnan')];
        lickv1v4_4 = cat(2, lickv1v4_4, mean(v1v4nMov(ifov).lick_v4,2,'omitnan')); % [lickv1v4_4 mean(v1v4nMov(ifov).lick_v4,2,'omitnan')];

        movcorrs_v1v4 = cat(1, movcorrs_v1v4, v1v4nMov(ifov).movcorrs); % mean(v1v4nMov(ifov).spk_v4,3,'omitnan')];
        % sig_v1v4 = cat(1, sig_v1v4, [v1v4nMov(ifov).sig_v1 v1v4nMov(ifov).sig_v4]);
        sig_v1v4 = cat(1, sig_v1v4, [sum(v1v4nMov(ifov).sig_v1) sum(v1v4nMov(ifov).sig_v4)] / numel(v1v4nMov(ifov).sig_v1));
    end
end

summary_v1v4 = v2struct(spkv1v4_1, movv1v4_1, lickv1v4_1, spkv1v4_4, movv1v4_4, lickv1v4_4, movcorrs_v1v4, sig_v1v4);

spkv41v4n_1 = [];
movv41v4n_1 = [];
lickv41v4n_1 = [];
spkv41v4n_n = [];
movv41v4n_n = [];
lickv41v4n_n = [];
movcorrs_v41v4n = [];
sig_v41v4n = [];

for ifov = 1:numel(v41v4nMov)
    if v41v4nMov(ifov).nMatch >= nMin
        spkv41v4n_1 = cat(2, spkv41v4n_1, mean(v41v4nMov(ifov).spk_v41,3,'omitnan')); % [spkv1v4_1 mean(v1v4nMov(ifov).spk_v1,3,'omitnan')];
        movv41v4n_1 = cat(2, movv41v4n_1, mean(v41v4nMov(ifov).mov_v41,2,'omitnan')); % [movv1v4_1 mean(v1v4nMov(ifov).mov_v1,2,'omitnan')];
        lickv41v4n_1 = cat(2, lickv41v4n_1, mean(v41v4nMov(ifov).lick_v41,2,'omitnan')); % [lickv1v4_1 mean(v1v4nMov(ifov).lick_v1,2,'omitnan')];

        spkv41v4n_n = cat(2, spkv41v4n_n, mean(v41v4nMov(ifov).spk_v4n,3,'omitnan')); % [spkv1v4_4 mean(v1v4nMov(ifov).spk_v4,3,'omitnan')];
        movv41v4n_n = cat(2, movv41v4n_n, mean(v41v4nMov(ifov).mov_v4n,2,'omitnan')); % [movv1v4_4 mean(v1v4nMov(ifov).mov_v4,2,'omitnan')];
        lickv41v4n_n = cat(2, lickv41v4n_n, mean(v41v4nMov(ifov).lick_v4n,2,'omitnan')); % [lickv1v4_4 mean(v1v4nMov(ifov).lick_v4,2,'omitnan')];

        movcorrs_v41v4n = cat(1, movcorrs_v41v4n, v41v4nMov(ifov).movcorrs); % mean(v1v4nMov(ifov).spk_v4,3,'omitnan')];
        % sig_v1v4 = cat(1, sig_v1v4, [v1v4nMov(ifov).sig_v1 v1v4nMov(ifov).sig_v4]);
        sig_v41v4n = cat(1, sig_v41v4n, [sum(v41v4nMov(ifov).sig_v41) sum(v41v4nMov(ifov).sig_v4n)] / numel(v41v4nMov(ifov).sig_v41));
    end
end

summary_v41v4n = v2struct(spkv41v4n_1, movv41v4n_1, lickv41v4n_1, spkv41v4n_n, movv41v4n_n, lickv41v4n_n, movcorrs_v41v4n, sig_v41v4n);

% Stats:
[pCorr, hCorr] = ranksum(movcorrs_v1v4(:,1), movcorrs_v41v4n(:,1));
[pSig14, hSig14] = signrank(sig_v1v4(:,1),sig_v1v4(:,2));
[pSig414n, hSig414n] = signrank(sig_v41v4n(:,1),sig_v41v4n(:,2));

statsAll = v2struct(pCorr,hCorr,pSig14,hSig14,pSig414n,hSig414n);

% Pca sorting based on v4n
coeff = pca(spkv1v4_4(end/2-imfs/2+1:end/2+imfs/2,:));
[~, iSort_v1v4] = sort(coeff(:,1),'descend');
coeff = pca(spkv41v4n_n(end/2-imfs/2+1:end/2+imfs/2,:));
[~, iSort_v41v4n] = sort(coeff(:,1),'descend');


% Plot in one big figure
f2 = figure;
sgtitle('Summary lobv + sim2: matched movement aligned activity','FontSize',24);
% v1-v4: v1 stuff
spkdata = spkv1v4_1(:,iSort_v1v4);
movdata = movv1v4_1;
lickdata = lickv1v4_1;
subplot(2,4,1);
nroi = size(spkdata,2);
imagesc(t2p,[1,nroi],mean(spkdata,3)',[0 10]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel('PC dendrite number')

title('v1-v4n: v1')
% Mean recording
subplot(6,4,13);
h = boundedline(t2p, mean(mean(spkdata,3),2), std(mean(spkdata,3),[],2),'-k','alpha'); set(h, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -3 12]); yticks([-3:3:12])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel ('Event rate (spk/s)')

% Licking
subplot(6,4,17);
h = boundedline(tpaq,mean(lickdata,2),std(lickdata,[],2),'-k','alpha'); set(h, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -.1 1.1]); yticks([0:.25:1])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel ('Licking prob.')

% Movement
subplot(6,4,21);
h = boundedline(tpaq,mean(movdata,2),std(movdata,[],2),'-k','alpha'); set(h, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2.5 7.5]); yticks([-2.5:2.5:7.5])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel ({'Wheel movement','(cm/s)'})

% v1-v4: v4 stuff
spkdata = spkv1v4_4(:,iSort_v1v4);
movdata = movv1v4_4;
lickdata = lickv1v4_4;
subplot(2,4,2);
nroi = size(spkdata,2);
imagesc(t2p,[1,nroi],mean(spkdata,3)',[0 10]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('v1-v4n: v4')
% Mean recording
subplot(6,4,14);
h = boundedline(t2p, mean(mean(spkdata,3),2), std(mean(spkdata,3),[],2),'-k','alpha'); set(h, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -3 12]); yticks([-3:3:12])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
% Licking
subplot(6,4,18);
h = boundedline(tpaq,mean(lickdata,2),std(lickdata,[],2),'-k','alpha'); set(h, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -.1 1.1]); yticks([0:.25:1])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
% Movement
subplot(6,4,22);
h = boundedline(tpaq,mean(movdata,2),std(movdata,[],2),'-k','alpha'); set(h, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2.5 7.5]); yticks([-2.5:2.5:7.5])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');

% v41-v4: v41 stuff
spkdata = spkv41v4n_1(:,iSort_v41v4n);
movdata = movv41v4n_1;
lickdata = lickv41v4n_1;
subplot(2,4,3);
nroi = size(spkdata,2);
imagesc(t2p,[1,nroi],mean(spkdata,3)',[0 10]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('v41-v4n: v41')
% Mean recording
subplot(6,4,15);
h = boundedline(t2p, mean(mean(spkdata,3),2), std(mean(spkdata,3),[],2),'-k','alpha'); set(h, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -3 12]); yticks([-3:3:12])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
% Licking
subplot(6,4,19);
h = boundedline(tpaq,mean(lickdata,2),std(lickdata,[],2),'-k','alpha'); set(h, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -.1 1.1]); yticks([0:.25:1])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
% Movement
subplot(6,4,23);
h = boundedline(tpaq,mean(movdata,2),std(movdata,[],2),'-k','alpha'); set(h, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2.5 7.5]); yticks([-2.5:2.5:7.5])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');

% v41-v4n: v4n stuff
spkdata = spkv41v4n_n(:,iSort_v41v4n);
movdata = movv41v4n_n;
lickdata = lickv41v4n_n;
subplot(2,4,4);
nroi = size(spkdata,2);
imagesc(t2p,[1,nroi],mean(spkdata,3)',[0 10]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('v41-v4n: v4n')
% Mean recording
subplot(6,4,16);
h = boundedline(t2p, mean(mean(spkdata,3),2), std(mean(spkdata,3),[],2),'-k','alpha'); set(h, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -3 12]); yticks([-3:3:12])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
% Licking
subplot(6,4,20);
h = boundedline(tpaq,mean(lickdata,2),std(lickdata,[],2),'-k','alpha'); set(h, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -.1 1.1]); yticks([0:.25:1])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
% Movement
subplot(6,4,24);
h = boundedline(tpaq,mean(movdata,2),std(movdata,[],2),'-k','alpha'); set(h, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2.5 7.5]); yticks([-2.5:2.5:7.5])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');



if savebool
    figFold = fullfile(paperFold,'Figures','Fig4','Sub1_v1v41v4n_movCorrs_rev1');
    if ~exist(figFold,'dir'); mkdir(figFold); end
    figname = 'Summary_movCorrMatch';
    savefig(f2, fullfile(figFold,[figname,'.fig']));
    saveas(f2,fullfile(figFold,[figname,'.png']));
    print(f2,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');

    save(fullfile(figFold,'v1_41_4n_movMatch.mat'),'v1v4nMov', 'v41v4nMov', 'summary_v1v4','summary_v41v4n', 'statsAll');

end


end
