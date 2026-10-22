function fig3_Supp_traditionalRFs(dataFold, paperFold, glmParams, savebool)
% Loop through each dataset and make heatmaps task-aligned activity using
% RFs scripts - which should have isolate events

%% Settings stuff
% Load in v4 sigSummary
glmFold0 = fullfile(paperFold,'Figures','Fig3','Sub3_SummaryBars');
glmFold = fullfile(glmFold0,sprintf('%s_%s_alpha%i',glmParams.dataType,glmParams.distr,glmParams.alpha*100));
glmSummary = load(fullfile(glmFold,'sigSummary_v4n.mat')); glmSummary = glmSummary.glmSummary;

% RF parent folder
rfFold0 = fullfile(dataFold,'01_Behav+imaging','Version_4');

% Initialize new struct
seshv4rfs = struct('name',{glmSummary.name}, 'date',{glmSummary.date}, 'fov',{glmSummary.fov}, ...
    'sigBool', {glmSummary.sigBool}, 'sigBoolVis', {glmSummary.sigBoolVis}, 'sigBoolMov', {glmSummary.sigBoolMov}, 'sigBoolRew', {glmSummary.sigBoolRew}, 'sigBoolLick', {glmSummary.sigBoolLick});

sigBool_glm_lobv = cat(1,seshv4rfs(cat(1,glmSummary.iLobV)).sigBool);
sigBool_glm_sim2 = cat(1,seshv4rfs(cat(1,glmSummary.iSim2)).sigBool);
fprintf(sprintf('Sig Lobv ROIs: %i of %i, %.0f percent\n',sum(sigBool_glm_lobv),numel(sigBool_glm_lobv), 100 * sum(sigBool_glm_lobv) / numel(sigBool_glm_lobv)))
fprintf(sprintf('Sig Sim2 ROIs: %i of %i, %.0f percent\n',sum(sigBool_glm_sim2),numel(sigBool_glm_sim2), 100 * sum(sigBool_glm_sim2) / numel(sigBool_glm_sim2)))

% Make t-axes - assum 2 second around each event and 30, 100 hz sampling
imfs = 30;
paqfs = 100;
t2p = (1:120)'/imfs - 2 - 1/imfs/2;
tpaq = (1:400)'/paqfs - 2 - 1/paqfs/2;
%% Loop through individual datasets

% Make empty arrays for master matrices
onSpk_lobv_sig = []; onSpk_lobv_non = []; onMov_lobv = []; onLick_lobv = [];
movSpk_lobv_sig = []; movSpk_lobv_non = []; movMov_lobv = []; movLick_lobv = [];
rewSpk_lobv_sig = []; rewSpk_lobv_non = []; rewMov_lobv = []; rewLick_lobv = [];
lickSpk_lobv_sig = []; lickSpk_lobv_non = []; lickMov_lobv = []; lickLick_lobv = [];
onSpk_sim2_sig = []; onSpk_sim2_non = []; onMov_sim2 = []; onLick_sim2 = [];
movSpk_sim2_sig = []; movSpk_sim2_non = []; movMov_sim2 = []; movLick_sim2 = [];
rewSpk_sim2_sig = []; rewSpk_sim2_non = []; rewMov_sim2 = []; rewLick_sim2 = [];
lickSpk_sim2_sig = []; lickSpk_sim2_non = []; lickMov_sim2 = []; lickLick_sim2 = [];
onSpk_lobv = []; movSpk_lobv = []; rewSpk_lobv = []; lickSpk_lobv = [];
onSpk_sim2 = []; movSpk_sim2 = []; rewSpk_sim2 = []; lickSpk_sim2 = [];
onSpk_odd_lobv = []; movSpk_odd_lobv = []; rewSpk_odd_lobv = []; lickSpk_odd_lobv = [];
onSpk_odd_sim2 = []; movSpk_odd_sim2 = []; rewSpk_odd_sim2 = []; lickSpk_odd_sim2 = [];
onSpk_even_lobv = []; movSpk_even_lobv = []; rewSpk_even_lobv = []; lickSpk_even_lobv = [];
onSpk_even_sim2 = []; movSpk_even_sim2 = []; rewSpk_even_sim2 = []; lickSpk_even_sim2 = [];

% Loop through!
for iSesh = 1:numel(seshv4rfs)
    % Fix boolean zeros at the end of the sigBool arrays - should be
    % general for all arrays if an issue for one!
    if numel(seshv4rfs(iSesh).sigBool) > numel(seshv4rfs(iSesh).sigBoolVis)
        seshv4rfs(iSesh).sigBoolVis(numel(seshv4rfs(iSesh).sigBoolVis)+1:numel(seshv4rfs(iSesh).sigBool)) = false;
        seshv4rfs(iSesh).sigBoolMov(numel(seshv4rfs(iSesh).sigBoolMov)+1:numel(seshv4rfs(iSesh).sigBool)) = false;
        seshv4rfs(iSesh).sigBoolRew(numel(seshv4rfs(iSesh).sigBoolRew)+1:numel(seshv4rfs(iSesh).sigBool)) = false;
        seshv4rfs(iSesh).sigBoolLick(numel(seshv4rfs(iSesh).sigBoolLick)+1:numel(seshv4rfs(iSesh).sigBool)) = false;
    end


    % Load in RFs
    RFs = load(fullfile(rfFold0, seshv4rfs(iSesh).name,sprintf('%s_%s',seshv4rfs(iSesh).date,seshv4rfs(iSesh).fov), 'RFs', 'RFs.mat')); RFs = RFs.RFs;
    seshv4rfs(iSesh).onSpk = RFs.on.SPK;
    seshv4rfs(iSesh).onLick = RFs.on.lick;
    seshv4rfs(iSesh).onMov = -RFs.on.mov; % in cm/s
    seshv4rfs(iSesh).movSpk = RFs.mover.SPK;
    seshv4rfs(iSesh).movLick = RFs.mover.lick;
    seshv4rfs(iSesh).movMov = -RFs.mover.mov; % in cm/s
    seshv4rfs(iSesh).rewSpk = RFs.trialrew.SPK;
    seshv4rfs(iSesh).rewLick = RFs.trialrew.lick;
    seshv4rfs(iSesh).rewMov = -RFs.trialrew.mov; % in cm/s
    seshv4rfs(iSesh).lickSpk = RFs.lick.SPK;
    seshv4rfs(iSesh).lickLick = RFs.lick.lick;
    seshv4rfs(iSesh).lickMov = -RFs.lick.mov; % in cm/s    

    % On: sign-rank sig test
    spkmat = seshv4rfs(iSesh).onSpk;
    sigbool = false(size(spkmat,2),1);
    for jCell = 1:size(spkmat,2)
        spkBase = imfs*squeeze(mean(spkmat(end/2-imfs+1:end/2,jCell,:),1));
        spkMod =  imfs*squeeze(mean(spkmat(end/2+1:end/2+0.3*imfs,jCell,:),1));
        [a,b] = signrank(spkBase,spkMod);
        sigbool(jCell) = logical(b);
    end
    seshv4rfs(iSesh).sigBoolVis_wilc = sigbool;
    % mov: sign-rank sig test
    spkmat = seshv4rfs(iSesh).movSpk;
    sigbool = false(size(spkmat,2),1);
    for jCell = 1:size(spkmat,2)
        spkBase = imfs*squeeze(mean(spkmat(1:imfs,jCell,:),1));
        spkMod =  imfs*squeeze(mean(spkmat(end/2-0.3*imfs+1:end/2,jCell,:),1));
        [a,b] = signrank(spkBase,spkMod);
        sigbool(jCell) = logical(b);
    end
    seshv4rfs(iSesh).sigBoolMov_wilc = sigbool;
    % rew: sign-rank sig test
    spkmat = seshv4rfs(iSesh).rewSpk;
    sigbool = false(size(spkmat,2),1);
    for jCell = 1:size(spkmat,2)
        spkBase = imfs*squeeze(mean(spkmat(end/2-imfs+1:end/2,jCell,:),1));
        spkMod =  imfs*squeeze(mean(spkmat(end/2+1:end/2+0.2*imfs,jCell,:),1));
        [a,b] = signrank(spkBase,spkMod);
        sigbool(jCell) = logical(b);
    end
    seshv4rfs(iSesh).sigBoolRew_wilc = sigbool;
    % lick: sign-rank sig test
    spkmat = seshv4rfs(iSesh).lickSpk;
    sigbool = false(size(spkmat,2),1);
    for jCell = 1:size(spkmat,2)
        spkBase = imfs*squeeze(mean(spkmat(1:imfs,jCell,:),1));
        % spkMod =  imfs*squeeze(mean(spkmat(end/2-0.2*imfs+1:end/2+0.2*imfs,jCell,:),1));
        spkMod =  imfs*squeeze(mean(spkmat(end/2-0.2*imfs+1:end/2,jCell,:),1));
        [a,b] = signrank(spkBase,spkMod);
        sigbool(jCell) = logical(b);
    end
    seshv4rfs(iSesh).sigBoolLick_wilc = sigbool;

    % % Process on data:
    onSpk = imfs*mean(seshv4rfs(iSesh).onSpk,3);
    onSpk_odd = imfs*mean(seshv4rfs(iSesh).onSpk(:,:,1:2:end),3);
    onSpk_even = imfs*mean(seshv4rfs(iSesh).onSpk(:,:,2:2:end),3);
    onSpk_sig = imfs*mean(seshv4rfs(iSesh).onSpk(:, seshv4rfs(iSesh).sigBoolVis,:),3);
    onSpk_sig_mu = mean(onSpk_sig, 2);
    onSpk_sig_sd = std(onSpk_sig, [], 2);
    onSpk_non = imfs*mean(seshv4rfs(iSesh).onSpk(:, ~seshv4rfs(iSesh).sigBoolVis,:),3);
    onSpk_non_mu = mean(onSpk_non, 2);
    onSpk_non_sd = std(onSpk_non, [], 2);
    onLick_mu = mean(seshv4rfs(iSesh).onLick, 2);
    onLick_sd = std(seshv4rfs(iSesh).onLick, [], 2);
    onMov_mu = mean(seshv4rfs(iSesh).onMov, 2);
    onMov_sd = std(seshv4rfs(iSesh).onMov, [], 2);

    % % Process mov data:
    movSpk = imfs*mean(seshv4rfs(iSesh).movSpk,3);
    movSpk_odd = imfs*mean(seshv4rfs(iSesh).movSpk(:,:,1:2:end),3);
    movSpk_even = imfs*mean(seshv4rfs(iSesh).movSpk(:,:,2:2:end),3);
    movSpk_sig = imfs*mean(seshv4rfs(iSesh).movSpk(:, seshv4rfs(iSesh).sigBoolMov,:),3);
    movSpk_sig_mu = mean(movSpk_sig, 2);
    movSpk_sig_sd = std(movSpk_sig, [], 2);
    movSpk_non = imfs*mean(seshv4rfs(iSesh).movSpk(:, ~seshv4rfs(iSesh).sigBoolMov,:),3);
    movSpk_non_mu = mean(movSpk_non, 2);
    movSpk_non_sd = std(movSpk_non, [], 2);
    movLick_mu = mean(seshv4rfs(iSesh).movLick, 2);
    movLick_sd = std(seshv4rfs(iSesh).movLick, [], 2);
    movMov_mu = mean(seshv4rfs(iSesh).movMov, 2);
    movMov_sd = std(seshv4rfs(iSesh).movMov, [], 2);    

    % % Process rew data:
    rewSpk = imfs*mean(seshv4rfs(iSesh).rewSpk,3);
    rewSpk_odd = imfs*mean(seshv4rfs(iSesh).rewSpk(:,:,1:2:end),3);
    rewSpk_even = imfs*mean(seshv4rfs(iSesh).rewSpk(:,:,1:2:end),3);
    rewSpk_sig = imfs*mean(seshv4rfs(iSesh).rewSpk(:, seshv4rfs(iSesh).sigBoolRew,:),3);
    rewSpk_sig_mu = mean(rewSpk_sig, 2);
    rewSpk_sig_sd = std(rewSpk_sig, [], 2);
    rewSpk_non = imfs*mean(seshv4rfs(iSesh).rewSpk(:, ~seshv4rfs(iSesh).sigBoolRew,:),3);
    rewSpk_non_mu = mean(rewSpk_non, 2);
    rewSpk_non_sd = std(rewSpk_non, [], 2);
    rewLick_mu = mean(seshv4rfs(iSesh).rewLick, 2);
    rewLick_sd = std(seshv4rfs(iSesh).rewLick, [], 2);
    rewMov_mu = mean(seshv4rfs(iSesh).rewMov, 2);
    rewMov_sd = std(seshv4rfs(iSesh).rewMov, [], 2);    

    % % Process lick data:
    lickSpk = imfs*mean(seshv4rfs(iSesh).lickSpk,3);
    lickSpk_odd = imfs*mean(seshv4rfs(iSesh).lickSpk(:,:,1:2:end),3);
    lickSpk_even = imfs*mean(seshv4rfs(iSesh).lickSpk(:,:,2:2:end),3);
    lickSpk_sig = imfs*mean(seshv4rfs(iSesh).lickSpk(:, seshv4rfs(iSesh).sigBoolLick,:),3);
    lickSpk_sig_mu = mean(lickSpk_sig, 2);
    lickSpk_sig_sd = std(lickSpk_sig, [], 2);
    lickSpk_non = imfs*mean(seshv4rfs(iSesh).lickSpk(:, ~seshv4rfs(iSesh).sigBoolLick,:),3);
    lickSpk_non_mu = mean(lickSpk_non, 2);
    lickSpk_non_sd = std(lickSpk_non, [], 2);
    lickLick_mu = mean(seshv4rfs(iSesh).lickLick, 2);
    lickLick_sd = std(seshv4rfs(iSesh).lickLick, [], 2);
    lickMov_mu = mean(seshv4rfs(iSesh).lickMov, 2);
    lickMov_sd = std(seshv4rfs(iSesh).lickMov, [], 2);

    % Populate master matrices depending of lobule
    if strcmpi(seshv4rfs(iSesh).fov, 'lobv')
        onSpk_lobv = cat(2, onSpk_lobv, onSpk);
        onSpk_odd_lobv = cat(2, onSpk_odd_lobv, onSpk_odd);
        onSpk_even_lobv = cat(2, onSpk_even_lobv, onSpk_even);
        onSpk_lobv_sig = cat(2, onSpk_lobv_sig, onSpk_sig);
        onSpk_lobv_non = cat(2, onSpk_lobv_non, onSpk_non);
        onMov_lobv = cat(2, onMov_lobv, onMov_mu);
        onLick_lobv = cat(2, onLick_lobv, onLick_mu);
        movSpk_lobv = cat(2, movSpk_lobv, movSpk);
        movSpk_odd_lobv = cat(2, movSpk_odd_lobv, movSpk_odd);
        movSpk_even_lobv = cat(2, movSpk_even_lobv, movSpk_even);
        movSpk_lobv_sig = cat(2, movSpk_lobv_sig, movSpk_sig);
        movSpk_lobv_non = cat(2, movSpk_lobv_non, movSpk_non);
        movMov_lobv = cat(2, movMov_lobv, movMov_mu);
        movLick_lobv = cat(2, movLick_lobv, movLick_mu);
        rewSpk_lobv = cat(2, rewSpk_lobv, rewSpk);
        rewSpk_odd_lobv = cat(2, rewSpk_odd_lobv, rewSpk_odd);
        rewSpk_even_lobv = cat(2, rewSpk_even_lobv, rewSpk_even);
        rewSpk_lobv_sig = cat(2, rewSpk_lobv_sig, rewSpk_sig);
        rewSpk_lobv_non = cat(2, rewSpk_lobv_non, rewSpk_non);
        rewMov_lobv = cat(2, rewMov_lobv, rewMov_mu);
        rewLick_lobv = cat(2, rewLick_lobv, rewLick_mu);
        lickSpk_lobv = cat(2, lickSpk_lobv, lickSpk);
        lickSpk_odd_lobv = cat(2, lickSpk_odd_lobv, lickSpk_odd);
        lickSpk_even_lobv = cat(2, lickSpk_even_lobv, lickSpk_even);
        lickSpk_lobv_sig = cat(2, lickSpk_lobv_sig, lickSpk_sig);
        lickSpk_lobv_non = cat(2, lickSpk_lobv_non, lickSpk_non);
        lickMov_lobv = cat(2, lickMov_lobv, lickMov_mu);
        lickLick_lobv = cat(2, lickLick_lobv, lickLick_mu);
    elseif strcmpi(seshv4rfs(iSesh).fov, 'sim2')
        onSpk_sim2 = cat(2, onSpk_sim2, onSpk);
        onSpk_odd_sim2 = cat(2, onSpk_odd_sim2, onSpk_odd);
        onSpk_even_sim2 = cat(2, onSpk_even_sim2, onSpk_even);
        onSpk_sim2_sig = cat(2, onSpk_sim2_sig, onSpk_sig);
        onSpk_sim2_non = cat(2, onSpk_sim2_non, onSpk_non);
        onMov_sim2 = cat(2, onMov_sim2, onMov_mu);
        onLick_sim2 = cat(2, onLick_sim2, onLick_mu);
        movSpk_sim2 = cat(2, movSpk_sim2, movSpk);
        movSpk_odd_sim2 = cat(2, movSpk_odd_sim2, movSpk_odd);
        movSpk_even_sim2 = cat(2, movSpk_even_sim2, movSpk_even);
        movSpk_sim2_sig = cat(2, movSpk_sim2_sig, movSpk_sig);
        movSpk_sim2_non = cat(2, movSpk_sim2_non, movSpk_non);
        movMov_sim2 = cat(2, movMov_sim2, movMov_mu);
        movLick_sim2 = cat(2, movLick_sim2, movLick_mu);
        rewSpk_sim2 = cat(2,rewSpk_sim2, rewSpk);
        rewSpk_odd_sim2 = cat(2,rewSpk_odd_sim2, rewSpk_odd);
        rewSpk_even_sim2 = cat(2,rewSpk_even_sim2, rewSpk_even);
        rewSpk_sim2_sig = cat(2, rewSpk_sim2_sig, rewSpk_sig);
        rewSpk_sim2_non = cat(2, rewSpk_sim2_non, rewSpk_non);
        rewMov_sim2 = cat(2, rewMov_sim2, rewMov_mu);
        rewLick_sim2 = cat(2, rewLick_sim2, rewLick_mu);
        lickSpk_sim2 = cat(2,lickSpk_sim2, lickSpk);
        lickSpk_odd_sim2 = cat(2,lickSpk_odd_sim2, lickSpk_odd);
        lickSpk_even_sim2 = cat(2,lickSpk_even_sim2, lickSpk_even);
        lickSpk_sim2_sig = cat(2, lickSpk_sim2_sig, lickSpk_sig);
        lickSpk_sim2_non = cat(2, lickSpk_sim2_non, lickSpk_non);
        lickMov_sim2 = cat(2, lickMov_sim2, lickMov_mu);
        lickLick_sim2 = cat(2, lickLick_sim2, lickLick_mu);
    end


    % Plot heatmaps, average responses, licks, and wheel turns
    f1a = figure;    
    
    % Trial onsets:
    iRow = 1;
    subplot(7,4,[0,4]+iRow); % onSpk_sig
    imagesc(t2p,[1:size(onSpk_sig,2)],onSpk_sig',[0 8]); xlim([-1 1]);
    xline(0,'--w', 'LineWidth',1)
    colormap viridis; box off;
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel(sprintf('Sigs (n = %i)',size(onSpk_sig,2)))
    title('Trial onset')
    
    subplot(7,4,[8,12]+iRow); % onSpk_non
    imagesc(t2p,[1:size(onSpk_non,2)],onSpk_non',[0 8]); xlim([-1 1]);
    xline(0,'--w', 'LineWidth',1)
    colormap viridis; box off;
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel(sprintf('Non-sigs (n = %i)',size(onSpk_non,2)))
    
    subplot(7,4,[16]+iRow); % onSpk_traces
    h1_1 = boundedline(t2p,onSpk_sig_mu, onSpk_sig_sd, '-k','alpha'); set(h1_1, 'linewidth', 1.5)
    h1_2 = boundedline(t2p,onSpk_non_mu, onSpk_non_sd, '-b','alpha'); set(h1_2, 'linewidth', 1.5)
    xline(0,'--k', 'LineWidth',1)
    axis([-1 1 -2 8]); yticks([-2:2:8])
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel ('Event rate (spk/s)')
    
    subplot(7,4,[20]+iRow); % onLick_traces
    h1_3 = boundedline(tpaq,onLick_mu, onLick_sd, '-k','alpha'); set(h1_3, 'linewidth', 1.5)
    xline(0,'--k', 'LineWidth',1)
    axis([-2 2 -.1 1.1]); yticks([0:.25:1])
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel('Licking')

    subplot(7,4,[24]+iRow); % onMov_traces
    h1_4 = boundedline(tpaq,onMov_mu, onMov_sd, '-k','alpha'); set(h1_4, 'linewidth', 1.5)
    xline(0,'--k', 'LineWidth',1)
    axis([-2 2 -2.5 12.5]); yticks([-2.5:2.5:12.5])
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    xlabel('Time from event (s)'); ylabel('Wheel (cm/s)')

    % Wheel movements:
    iRow = 2;

    subplot(7,4,[0,4]+iRow); % movSpk_sig
    imagesc(t2p,[1:size(movSpk_sig,2)],movSpk_sig',[0 8]); xlim([-1 1]);
    xline(0,'--w', 'LineWidth',1)
    colormap viridis; box off;
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel(sprintf('Sigs (n = %i)',size(movSpk_sig,2)))
    title('Wheel movement')
    
    subplot(7,4,[8,12]+iRow); % movSpk_non
    imagesc(t2p,[1:size(movSpk_non,2)],movSpk_non',[0 8]); xlim([-1 1]);
    xline(0,'--w', 'LineWidth',1)
    colormap viridis; box off;
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel(sprintf('Non-sigs (n = %i)',size(movSpk_non,2)))
    
    subplot(7,4,[16]+iRow); % movSpk_traces
    h2_1 = boundedline(t2p,movSpk_sig_mu, movSpk_sig_sd, '-k','alpha'); set(h2_1, 'linewidth', 1.5)
    h2_2 = boundedline(t2p,movSpk_non_mu, movSpk_non_sd, '-b','alpha'); set(h2_2, 'linewidth', 1.5)
    xline(0,'--k', 'LineWidth',1)
    axis([-1 1 -2 8]); yticks([-2:2:8])
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel ('Event rate (spk/s)')
    
    subplot(7,4,[20]+iRow); % movLick_traces
    h2_3 = boundedline(tpaq,movLick_mu, movLick_sd, '-k','alpha'); set(h2_3, 'linewidth', 1.5)
    xline(0,'--k', 'LineWidth',1)
    axis([-2 2 -.1 1.1]); yticks([0:.25:1])
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel('Licking')

    subplot(7,4,[24]+iRow); % movMov_traces
    h2_4 = boundedline(tpaq,movMov_mu, movMov_sd, '-k','alpha'); set(h2_4, 'linewidth', 1.5)
    xline(0,'--k', 'LineWidth',1)
    axis([-2 2 -2.5 12.5]); yticks([-2.5:2.5:12.5])
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    xlabel('Time from event (s)'); ylabel('Wheel (cm/s)')

    % Rewards:
    iRow = 3;

    subplot(7,4,[0,4]+iRow); % rewSpk_sig
    imagesc(t2p,[1:size(rewSpk_sig,2)],rewSpk_sig',[0 8]); xlim([-1 1]);
    xline(0,'--w', 'LineWidth',1)
    colormap viridis; box off;
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel(sprintf('Sigs (n = %i)',size(rewSpk_sig,2)))
    title('Trial reward')
    
    subplot(7,4,[8,12]+iRow); % rewSpk_non
    imagesc(t2p,[1:size(rewSpk_non,2)],rewSpk_non',[0 8]); xlim([-1 1]);
    xline(0,'--w', 'LineWidth',1)
    colormap viridis; box off;
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel(sprintf('Non-sigs (n = %i)',size(rewSpk_non,2)))
    
    subplot(7,4,[16]+iRow); % rewSpk_traces
    h3_1 = boundedline(t2p,rewSpk_sig_mu, rewSpk_sig_sd, '-k','alpha'); set(h3_1, 'linewidth', 1.5)
    h3_2 = boundedline(t2p,rewSpk_non_mu, rewSpk_non_sd, '-b','alpha'); set(h3_2, 'linewidth', 1.5)
    xline(0,'--k', 'LineWidth',1)
    axis([-1 1 -2 8]); yticks([-2:2:8])
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel ('Event rate (spk/s)')
    
    subplot(7,4,[20]+iRow); % rewLick_traces
    h3_3 = boundedline(tpaq,rewLick_mu, rewLick_sd, '-k','alpha'); set(h3_3, 'linewidth', 1.5)
    xline(0,'--k', 'LineWidth',1)
    axis([-2 2 -.1 1.1]); yticks([0:.25:1])
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel('Licking')

    subplot(7,4,[24]+iRow); % rewMov_traces
    h3_4 = boundedline(tpaq,rewMov_mu, rewMov_sd, '-k','alpha'); set(h3_4, 'linewidth', 1.5)
    xline(0,'--k', 'LineWidth',1)
    axis([-2 2 -2.5 12.5]); yticks([-2.5:2.5:12.5])
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    xlabel('Time from event (s)'); ylabel('Wheel (cm/s)')

    % Licking:
    iRow = 4;

    subplot(7,4,[0,4]+iRow); % lickSpk_sig
    imagesc(t2p,[1:size(lickSpk_sig,2)],lickSpk_sig',[0 8]); xlim([-1 1]);
    xline(0,'--w', 'LineWidth',1)
    colormap viridis; box off;
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel(sprintf('Sigs (n = %i)',size(lickSpk_sig,2)))
    title('Licks')
    
    subplot(7,4,[8,12]+iRow); % lickSpk_non
    imagesc(t2p,[1:size(lickSpk_non,2)],lickSpk_non',[0 8]); xlim([-1 1]);
    xline(0,'--w', 'LineWidth',1)
    colormap viridis; box off;
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel(sprintf('Non-sigs (n = %i)',size(lickSpk_non,2)))
    
    subplot(7,4,[16]+iRow); % lickSpk_traces
    h4_1 = boundedline(t2p,lickSpk_sig_mu, lickSpk_sig_sd, '-k','alpha'); set(h4_1, 'linewidth', 1.5)
    h4_2 = boundedline(t2p,lickSpk_non_mu, lickSpk_non_sd, '-b','alpha'); set(h4_2, 'linewidth', 1.5)
    xline(0,'--k', 'LineWidth',1)
    axis([-1 1 -2 8]); yticks([-2:2:8])
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel ('Event rate (spk/s)')
    
    subplot(7,4,[20]+iRow); % lickLick_traces
    h4_3 = boundedline(tpaq,lickLick_mu, lickLick_sd, '-k','alpha'); set(h4_3, 'linewidth', 1.5)
    xline(0,'--k', 'LineWidth',1)
    axis([-2 2 -.1 1.1]); yticks([0:.25:1])
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel('Licking')

    subplot(7,4,[24]+iRow); % lickMov_traces
    h4_4 = boundedline(tpaq,lickMov_mu, lickMov_sd, '-k','alpha'); set(h4_4, 'linewidth', 1.5)
    xline(0,'--k', 'LineWidth',1)
    axis([-2 2 -2.5 12.5]); yticks([-2.5:2.5:12.5])
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    xlabel('Time from event (s)'); ylabel('Wheel (cm/s)')
    
    suptitle(sprintf('%s: %s - Expert receptive fields (GLM)',seshv4rfs(iSesh).name,seshv4rfs(iSesh).fov))    

    % Plot heatmaps, average responses, licks, and wheel turns - wilcoxon
    f1b = figure;    
    
    % Trial onsets:
    iRow = 1;
    subplot(7,4,[0,4]+iRow); % onSpk_sig
    spklocal = onSpk(:,seshv4rfs(iSesh).sigBoolVis_wilc);
    imagesc(t2p,[1:size(spklocal,2)],spklocal',[0 8]); xlim([-1 1]);
    xline(0,'--w', 'LineWidth',1)
    colormap viridis; box off;
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel(sprintf('Sigs (n = %i)',size(spklocal,2)))
    title('Trial onset')
    
    subplot(7,4,[8,12]+iRow); % onSpk_non
    spklocal2 = onSpk(:,~seshv4rfs(iSesh).sigBoolVis_wilc);
    imagesc(t2p,[1:size(spklocal2,2)],spklocal2',[0 8]); xlim([-1 1]);
    xline(0,'--w', 'LineWidth',1)
    colormap viridis; box off;
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel(sprintf('Non-sigs (n = %i)',size(spklocal2,2)))
    
    subplot(7,4,[16]+iRow); % onSpk_traces
    h1_1 = boundedline(t2p, mean(spklocal,2), std(spklocal, [], 2), '-k','alpha'); set(h1_1, 'linewidth', 1.5)
    h1_2 = boundedline(t2p, mean(spklocal2,2), std(spklocal2, [], 2), '-b','alpha'); set(h1_2, 'linewidth', 1.5)
    xline(0,'--k', 'LineWidth',1)
    axis([-1 1 -2 8]); yticks([-2:2:8])
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel ('Event rate (spk/s)')
    
    subplot(7,4,[20]+iRow); % onLick_traces
    h1_3 = boundedline(tpaq,onLick_mu, onLick_sd, '-k','alpha'); set(h1_3, 'linewidth', 1.5)
    xline(0,'--k', 'LineWidth',1)
    axis([-2 2 -.1 1.1]); yticks([0:.25:1])
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel('Licking')

    subplot(7,4,[24]+iRow); % onMov_traces
    h1_4 = boundedline(tpaq,onMov_mu, onMov_sd, '-k','alpha'); set(h1_4, 'linewidth', 1.5)
    xline(0,'--k', 'LineWidth',1)
    axis([-2 2 -2.5 12.5]); yticks([-2.5:2.5:12.5])
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    xlabel('Time from event (s)'); ylabel('Wheel (cm/s)')

    % Wheel movements:
    iRow = 2;
    subplot(7,4,[0,4]+iRow); % movSpk_sig
    spklocal = movSpk(:,seshv4rfs(iSesh).sigBoolMov_wilc);
    imagesc(t2p,[1:size(spklocal,2)],spklocal',[0 8]); xlim([-1 1]);
    xline(0,'--w', 'LineWidth',1)
    colormap viridis; box off;
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel(sprintf('Sigs (n = %i)',size(spklocal,2)))
    title('Wheel movement')
    
    subplot(7,4,[8,12]+iRow); % movSpk_non
    spklocal2 = movSpk(:,~seshv4rfs(iSesh).sigBoolMov_wilc);
    imagesc(t2p,[1:size(spklocal2,2)],spklocal2',[0 8]); xlim([-1 1]);    
    xline(0,'--w', 'LineWidth',1)
    colormap viridis; box off;
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel(sprintf('Non-sigs (n = %i)',size(spklocal2,2)))
    
    subplot(7,4,[16]+iRow); % movSpk_traces
    h2_1 = boundedline(t2p, mean(spklocal,2), std(spklocal, [], 2), '-k','alpha'); set(h2_1, 'linewidth', 1.5)
    h2_2 = boundedline(t2p, mean(spklocal2,2), std(spklocal2, [], 2), '-b','alpha'); set(h2_2, 'linewidth', 1.5)
    xline(0,'--k', 'LineWidth',1)
    axis([-1 1 -2 8]); yticks([-2:2:8])
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel ('Event rate (spk/s)')
    
    subplot(7,4,[20]+iRow); % movLick_traces
    h2_3 = boundedline(tpaq,movLick_mu, movLick_sd, '-k','alpha'); set(h2_3, 'linewidth', 1.5)
    xline(0,'--k', 'LineWidth',1)
    axis([-2 2 -.1 1.1]); yticks([0:.25:1])
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel('Licking')

    subplot(7,4,[24]+iRow); % movMov_traces
    h2_4 = boundedline(tpaq,movMov_mu, movMov_sd, '-k','alpha'); set(h2_4, 'linewidth', 1.5)
    xline(0,'--k', 'LineWidth',1)
    axis([-2 2 -2.5 12.5]); yticks([-2.5:2.5:12.5])
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    xlabel('Time from event (s)'); ylabel('Wheel (cm/s)')

    % Rewards:
    iRow = 3;
    subplot(7,4,[0,4]+iRow); % rewSpk_sig
    spklocal = rewSpk(:,seshv4rfs(iSesh).sigBoolRew_wilc);
    imagesc(t2p,[1:size(spklocal,2)],spklocal',[0 8]); xlim([-1 1]);
    xline(0,'--w', 'LineWidth',1)
    colormap viridis; box off;
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel(sprintf('Sigs (n = %i)',size(spklocal,2)))
    title('Trial reward')
    
    subplot(7,4,[8,12]+iRow); % rewSpk_non
    spklocal2 = rewSpk(:,~seshv4rfs(iSesh).sigBoolRew_wilc);
    imagesc(t2p,[1:size(spklocal2,2)],spklocal2',[0 8]); xlim([-1 1]);        
    xline(0,'--w', 'LineWidth',1)
    colormap viridis; box off;
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel(sprintf('Non-sigs (n = %i)',size(spklocal2,2)))
    
    subplot(7,4,[16]+iRow); % rewSpk_traces
    h3_1 = boundedline(t2p, mean(spklocal,2), std(spklocal, [], 2), '-k','alpha'); set(h3_1, 'linewidth', 1.5)
    h3_2 = boundedline(t2p, mean(spklocal2,2), std(spklocal2, [], 2), '-b','alpha'); set(h3_2, 'linewidth', 1.5)
    xline(0,'--k', 'LineWidth',1)
    axis([-1 1 -2 8]); yticks([-2:2:8])
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel ('Event rate (spk/s)')
    
    subplot(7,4,[20]+iRow); % rewLick_traces
    h3_3 = boundedline(tpaq,rewLick_mu, rewLick_sd, '-k','alpha'); set(h3_3, 'linewidth', 1.5)
    xline(0,'--k', 'LineWidth',1)
    axis([-2 2 -.1 1.1]); yticks([0:.25:1])
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel('Licking')

    subplot(7,4,[24]+iRow); % rewMov_traces
    h3_4 = boundedline(tpaq,rewMov_mu, rewMov_sd, '-k','alpha'); set(h3_4, 'linewidth', 1.5)
    xline(0,'--k', 'LineWidth',1)
    axis([-2 2 -2.5 12.5]); yticks([-2.5:2.5:12.5])
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    xlabel('Time from event (s)'); ylabel('Wheel (cm/s)')

    % Licking:
    iRow = 4;
    subplot(7,4,[0,4]+iRow); % lickSpk_sig
    spklocal = lickSpk(:,seshv4rfs(iSesh).sigBoolLick_wilc);
    imagesc(t2p,[1:size(spklocal,2)],spklocal',[0 8]); xlim([-1 1]);    
    xline(0,'--w', 'LineWidth',1)
    colormap viridis; box off;
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel(sprintf('Sigs (n = %i)',size(spklocal,2)))
    title('Licks')
    
    subplot(7,4,[8,12]+iRow); % lickSpk_non
    spklocal2 = lickSpk(:,~seshv4rfs(iSesh).sigBoolLick_wilc);
    imagesc(t2p,[1:size(spklocal2,2)],spklocal2',[0 8]); xlim([-1 1]);
    xline(0,'--w', 'LineWidth',1)
    colormap viridis; box off;
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel(sprintf('Non-sigs (n = %i)',size(spklocal2,2)))
    
    subplot(7,4,[16]+iRow); % lickSpk_traces
    h4_1 = boundedline(t2p,mean(spklocal,2), std(spklocal, [], 2), '-k','alpha'); set(h4_1, 'linewidth', 1.5)
    h4_2 = boundedline(t2p,mean(spklocal2,2), std(spklocal2, [], 2), '-b','alpha'); set(h4_2, 'linewidth', 1.5)
    xline(0,'--k', 'LineWidth',1)
    axis([-1 1 -2 8]); yticks([-2:2:8])
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel ('Event rate (spk/s)')
    
    subplot(7,4,[20]+iRow); % lickLick_traces
    h4_3 = boundedline(tpaq,lickLick_mu, lickLick_sd, '-k','alpha'); set(h4_3, 'linewidth', 1.5)
    xline(0,'--k', 'LineWidth',1)
    axis([-2 2 -.1 1.1]); yticks([0:.25:1])
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    ylabel('Licking')

    subplot(7,4,[24]+iRow); % lickMov_traces
    h4_4 = boundedline(tpaq,lickMov_mu, lickMov_sd, '-k','alpha'); set(h4_4, 'linewidth', 1.5)
    xline(0,'--k', 'LineWidth',1)
    axis([-2 2 -2.5 12.5]); yticks([-2.5:2.5:12.5])
    set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
    xlabel('Time from event (s)'); ylabel('Wheel (cm/s)')
    
    suptitle(sprintf('%s: %s - Expert receptive fields (Wilcoxon)',seshv4rfs(iSesh).name,seshv4rfs(iSesh).fov))    

    if savebool
        figFold = fullfile(paperFold,'Figures','Fig3','SuppFig_receptiveFields');
        if ~exist(figFold,'dir'); mkdir(figFold); end
        figname = sprintf('%s_%s_expertRFs_GLM',seshv4rfs(iSesh).name,seshv4rfs(iSesh).fov);
        savefig(f1a, fullfile(figFold,[figname,'.fig']));
        saveas(f1a,fullfile(figFold,[figname,'.png']));
        print(f1a,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');
        figname = sprintf('%s_%s_expertRFs_wilc',seshv4rfs(iSesh).name,seshv4rfs(iSesh).fov);
        savefig(f1b, fullfile(figFold,[figname,'.fig']));
        saveas(f1b,fullfile(figFold,[figname,'.png']));
        print(f1b,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');
    end

    clear h1_1 h1_2 h1_3 h1_4 h2_1 h2_2 h2_3 h2_4 h3_1 h3_2 h3_3 h3_4 h4_1 h4_2 h4_3 h4_4 ...
        onLick_mu onLick_sd onMov_mu onMov_sd onSpk_sig onSpk_sig_mu onSpk_sig_sd onSpk_non onSpk_non_mu onSpk_non_sd ...
        movLick_mu movLick_sd movMov_mu movMov_sd movSpk_sig movSpk_sig_mu movSpk_sig_sd movSpk_non movSpk_non_mu movSpk_non_sd ...
        rewLick_mu rewLick_sd rewMov_mu rewMov_sd rewSpk_sig rewSpk_sig_mu rewSpk_sig_sd rewSpk_non rewSpk_non_mu rewSpk_non_sd ...
        lickLick_mu lickLick_sd lickMov_mu lickMov_sd lickSpk_sig lickSpk_sig_mu lickSpk_sig_sd lickMov_sd lickSpk_non lickSpk_non_mu lickSpk_non_sd
end

lobvSummary = v2struct(onSpk_lobv, movSpk_lobv, rewSpk_lobv, lickSpk_lobv,...
    onSpk_lobv_sig, onSpk_lobv_non, onMov_lobv, onLick_lobv, ...
    movSpk_lobv_sig, movSpk_lobv_non, movMov_lobv, movLick_lobv, ...
    rewSpk_lobv_sig, rewSpk_lobv_non, rewMov_lobv, rewLick_lobv, ...
    lickSpk_lobv_sig, lickSpk_lobv_non, lickMov_lobv, lickLick_lobv);

sim2Summary = v2struct(onSpk_sim2, movSpk_sim2, rewSpk_sim2, lickSpk_sim2,...
    onSpk_sim2_sig, onSpk_sim2_non, onMov_sim2, onLick_sim2, ...
    movSpk_sim2_sig, movSpk_sim2_non, movMov_sim2, movLick_sim2, ...
    rewSpk_sim2_sig, rewSpk_sim2_non, rewMov_sim2, rewLick_sim2, ...
    lickSpk_sim2_sig, lickSpk_sim2_non, lickMov_sim2, lickLick_sim2);

clear h1_1 h1_2 h1_3 h1_4 h2_1 h2_2 h2_3 h2_4 h3_1 h3_2 h3_3 h3_4 h4_1 h4_2 h4_3 h4_4 a b iRow iSesh jCell %...
    % lickLick_mu lickLick_sd lickMov_mu lickMov_sd lickMov_sim2 lickSpk lickSpk_even lickSpk_odd ...
    % lickSpk_non lickSpk_non_mu lickSpk_non_sd lickSpk_sig lickSpk_sig_mu lickSpk_sig_sd

%% Lobv sorting - GLM

f2a = figure;

sigBoolVis_glm_lobv = cat(1,seshv4rfs(cat(1,glmSummary.iLobV)).sigBoolVis);
sigBoolMov_glm_lobv = cat(1,seshv4rfs(cat(1,glmSummary.iLobV)).sigBoolMov);
sigBoolRew_glm_lobv = cat(1,seshv4rfs(cat(1,glmSummary.iLobV)).sigBoolRew);
sigBoolLick_glm_lobv = cat(1,seshv4rfs(cat(1,glmSummary.iLobV)).sigBoolLick);

% Trial onsets:
iRow = 1;
coeff = pca(onSpk_odd_lobv(end/2-imfs/2+1:end/2+imfs/2,:));
[~, iSort] = sort(coeff(:,1),'descend');
iSort_sig = sigBoolVis_glm_lobv(iSort);
iSort_non = ~sigBoolVis_glm_lobv(iSort);
onSpk_even_lobv2 = onSpk_even_lobv(:,iSort);

subplot(7,4,[0,4]+iRow); % onSpk_sig
spklocal_sig = onSpk_even_lobv2(:,iSort_sig);
imagesc(t2p,[1:size(spklocal_sig,2)],spklocal_sig',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Sigs (n = %i)',sum(iSort_sig)))
title('Trial onset')

subplot(7,4,[8,12]+iRow); % onSpk_non
spklocal_non = onSpk_even_lobv2(:,iSort_non);
imagesc(t2p,[1:size(spklocal_non,2)],spklocal_non',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Non-sigs (n = %i)',sum(iSort_non)));

subplot(7,4,[16]+iRow); % onSpk_traces
h1_1 = boundedline(t2p,mean(onSpk_lobv_sig,2), std(onSpk_lobv_sig,[],2), '-k','alpha'); set(h1_1, 'linewidth', 1.5)
h1_2 = boundedline(t2p,mean(onSpk_lobv_non,2), std(onSpk_lobv_non,[],2), '-b','alpha'); set(h1_2, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 8]); yticks([-2:2:8])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel ('Event rate (spk/s)')

subplot(7,4,[20]+iRow); % onLick_traces
h1_3 = boundedline(tpaq,mean(onLick_lobv,2), std(onLick_lobv,[],2), '-k','alpha'); set(h1_3, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -.1 1.1]); yticks([0:.25:1])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel('Licking')

subplot(7,4,[24]+iRow); % onMov_traces
h1_4 = boundedline(tpaq,mean(onMov_lobv,2), std(onMov_lobv,[],2), '-k','alpha'); set(h1_4, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 4]); yticks([-2:1:4]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlabel('Time from event (s)'); ylabel('Wheel (cm/s)')

% Wheel movements:
iRow = 2;
coeff = pca(movSpk_odd_lobv(end/2-imfs/2+1:end/2+imfs/2,:));
[~, iSort] = sort(coeff(:,1),'descend');
iSort_sig = sigBoolMov_glm_lobv(iSort);
iSort_non = ~sigBoolMov_glm_lobv(iSort);
movSpk_even_lobv2 = movSpk_even_lobv(:,iSort);

subplot(7,4,[0,4]+iRow); % movSpk_sig
spklocal_sig = movSpk_even_lobv2(:,iSort_sig);
imagesc(t2p,[1:size(spklocal_sig,2)],spklocal_sig',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Sigs (n = %i)',size(spklocal_sig,2)))
title('Wheel movement')

subplot(7,4,[8,12]+iRow); % movSpk_non
spklocal_non = movSpk_even_lobv2(:,iSort_non);
imagesc(t2p,[1:size(spklocal_non,2)],spklocal_non',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Non-sigs (n = %i)',size(spklocal_non,2)))

subplot(7,4,[16]+iRow); % movSpk_traces
h2_1 = boundedline(t2p,mean(movSpk_lobv_sig,2), std(movSpk_lobv_sig,[],2), '-k','alpha'); set(h2_1, 'linewidth', 1.5)
h2_2 = boundedline(t2p,mean(movSpk_lobv_non,2), std(movSpk_lobv_non,[],2), '-b','alpha'); set(h2_2, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 8]); yticks([-2:2:8])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel ('Event rate (spk/s)')

subplot(7,4,[20]+iRow); % movLick_traces
h2_3 = boundedline(tpaq,mean(movLick_lobv,2), std(movLick_lobv,[],2), '-k','alpha'); set(h2_3, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -.1 1.1]); yticks([0:.25:1])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel('Licking')

subplot(7,4,[24]+iRow); % movMov_traces
h2_4 = boundedline(tpaq,mean(movMov_lobv,2), std(movMov_lobv,[],2), '-k','alpha'); set(h2_4, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 4]); yticks([-2:1:4]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlabel('Time from event (s)'); ylabel('Wheel (cm/s)')

% Rewards:
iRow = 3;
coeff = pca(rewSpk_odd_lobv(end/2-imfs/2+1:end/2+imfs/2,:));
[~, iSort] = sort(coeff(:,1),'descend');
iSort_sig = sigBoolRew_glm_lobv(iSort);
iSort_non = ~sigBoolRew_glm_lobv(iSort);
rewSpk_even_lobv2 = rewSpk_even_lobv(:,iSort);

subplot(7,4,[0,4]+iRow); % rewSpk_sig
spklocal_sig = rewSpk_even_lobv2(:,iSort_sig);
imagesc(t2p,[1:size(spklocal_sig,2)],spklocal_sig',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Sigs (n = %i)',size(spklocal_sig,2)))
title('Trial reward')

subplot(7,4,[8,12]+iRow); % movSpk_non
spklocal_non = rewSpk_even_lobv2(:,iSort_non);
imagesc(t2p,[1:size(spklocal_non,2)],spklocal_non',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Non-sigs (n = %i)',size(spklocal_non,2)))

subplot(7,4,[16]+iRow); % rewSpk_traces
h3_1 = boundedline(t2p,mean(rewSpk_lobv_sig,2), std(rewSpk_lobv_sig,[],2), '-k','alpha'); set(h3_1, 'linewidth', 1.5)
h3_2 = boundedline(t2p,mean(rewSpk_lobv_non,2), std(rewSpk_lobv_non,[],2), '-b','alpha'); set(h3_2, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 8]); yticks([-2:2:8])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel ('Event rate (spk/s)')

subplot(7,4,[20]+iRow); % rewLick_traces
h3_3 = boundedline(tpaq,mean(rewLick_lobv,2), std(rewLick_lobv,[],2), '-k','alpha'); set(h3_3, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -.1 1.1]); yticks([0:.25:1])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel('Licking')

subplot(7,4,[24]+iRow); % rewMov_traces
h3_4 = boundedline(tpaq,mean(rewMov_lobv,2), std(rewMov_lobv,[],2), '-k','alpha'); set(h3_4, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 4]); yticks([-2:1:4]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlabel('Time from event (s)'); ylabel('Wheel (cm/s)')

% Licking:
iRow = 4;
coeff = pca(lickSpk_odd_lobv(end/2-imfs/2+1:end/2+imfs/2,:));
[~, iSort] = sort(coeff(:,1),'descend');
iSort_sig = sigBoolLick_glm_lobv(iSort);
iSort_non = ~sigBoolLick_glm_lobv(iSort);
lickSpk_even_lobv2 = lickSpk_even_lobv(:,iSort);

subplot(7,4,[0,4]+iRow); % lickSpk_sig
spklocal_sig = lickSpk_even_lobv2(:,iSort_sig);
imagesc(t2p,[1:size(spklocal_sig,2)],spklocal_sig',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Sigs (n = %i)',size(spklocal_sig,2)))
title('Licking')

subplot(7,4,[8,12]+iRow); % movSpk_non
spklocal_non = lickSpk_even_lobv2(:,iSort_non);
imagesc(t2p,[1:size(spklocal_non,2)],spklocal_non',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Non-sigs (n = %i)',size(spklocal_non,2)))

subplot(7,4,[16]+iRow); % lickSpk_traces
h4_1 = boundedline(t2p,mean(lickSpk_lobv_sig,2), std(lickSpk_lobv_sig,[],2), '-k','alpha'); set(h4_1, 'linewidth', 1.5)
h4_2 = boundedline(t2p,mean(lickSpk_lobv_non,2), std(lickSpk_lobv_non,[],2), '-b','alpha'); set(h4_2, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 8]); yticks([-2:2:8])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel ('Event rate (spk/s)')

subplot(7,4,[20]+iRow); % lickLick_traces
h4_3 = boundedline(tpaq,mean(lickLick_lobv,2), std(lickLick_lobv,[],2), '-k','alpha'); set(h4_3, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -.1 1.1]); yticks([0:.25:1])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel('Licking')

subplot(7,4,[24]+iRow); % lickMov_traces
h4_4 = boundedline(tpaq,mean(lickMov_lobv,2), std(lickMov_lobv,[],2), '-k','alpha'); set(h4_4, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 4]); yticks([-2:1:4]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlabel('Time from event (s)'); ylabel('Wheel (cm/s)')

suptitle('Lobule V - Expert receptive fields (GLM)')

%% Lobv sorting - Wilcoxon

f2b = figure;

sigBoolVis_wilc_lobv = cat(1,seshv4rfs(cat(1,glmSummary.iLobV)).sigBoolVis_wilc);
sigBoolMov_wilc_lobv = cat(1,seshv4rfs(cat(1,glmSummary.iLobV)).sigBoolMov_wilc);
sigBoolRew_wilc_lobv = cat(1,seshv4rfs(cat(1,glmSummary.iLobV)).sigBoolRew_wilc);
sigBoolLick_wilc_lobv = cat(1,seshv4rfs(cat(1,glmSummary.iLobV)).sigBoolLick_wilc);

% Trial onsets:
iRow = 1;
coeff = pca(onSpk_odd_lobv(end/2-imfs/2+1:end/2+imfs/2,:));
[~, iSort] = sort(coeff(:,1),'descend');
iSort_sig = sigBoolVis_wilc_lobv(iSort);
iSort_non = ~sigBoolVis_wilc_lobv(iSort);
onSpk_even_lobv2 = onSpk_even_lobv(:,iSort);

subplot(7,4,[0,4]+iRow); % onSpk_sig
spklocal_sig = onSpk_even_lobv2(:,iSort_sig);
imagesc(t2p,[1:size(spklocal_sig,2)],spklocal_sig',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Sigs (n = %i)',sum(iSort_sig)))
title('Trial onset')

subplot(7,4,[8,12]+iRow); % onSpk_non
spklocal_non = onSpk_even_lobv2(:,iSort_non);
imagesc(t2p,[1:size(spklocal_non,2)],spklocal_non',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Non-sigs (n = %i)',sum(iSort_non)))

subplot(7,4,[16]+iRow); % onSpk_traces
h1_1 = boundedline(t2p,mean(onSpk_lobv(:,sigBoolVis_wilc_lobv),2), std(onSpk_lobv(:,sigBoolVis_wilc_lobv),[],2), '-k','alpha'); set(h1_1, 'linewidth', 1.5)
h1_2 = boundedline(t2p,mean(onSpk_lobv(:,~sigBoolVis_wilc_lobv),2), std(onSpk_lobv(:,~sigBoolVis_wilc_lobv),[],2), '-b','alpha'); set(h1_2, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 8]); yticks([-2:2:8])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel ('Event rate (spk/s)')

subplot(7,4,[20]+iRow); % onLick_traces
h1_3 = boundedline(tpaq,mean(onLick_lobv,2), std(onLick_lobv,[],2), '-k','alpha'); set(h1_3, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -.1 1.1]); yticks([0:.25:1])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel('Licking')

subplot(7,4,[24]+iRow); % onMov_traces
h1_4 = boundedline(tpaq,mean(onMov_lobv,2), std(onMov_lobv,[],2), '-k','alpha'); set(h1_4, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 4]); yticks([-2:1:4]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlabel('Time from event (s)'); ylabel('Wheel (cm/s)')

% Wheel movements:
iRow = 2;
coeff = pca(movSpk_odd_lobv(end/2-imfs/2+1:end/2+imfs/2,:));
[~, iSort] = sort(coeff(:,1),'descend');
iSort_sig = sigBoolMov_wilc_lobv(iSort);
iSort_non = ~sigBoolMov_wilc_lobv(iSort);
movSpk_even_lobv2 = movSpk_even_lobv(:,iSort);

subplot(7,4,[0,4]+iRow); % movSpk_sig
spklocal_sig = movSpk_even_lobv2(:,iSort_sig);
imagesc(t2p,[1:size(spklocal_sig,2)],spklocal_sig',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Sigs (n = %i)',sum(iSort_sig)))
title('Wheel movement')

subplot(7,4,[8,12]+iRow); % movSpk_non
spklocal_non = movSpk_even_lobv2(:,iSort_non);
imagesc(t2p,[1:size(spklocal_non,2)],spklocal_non',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Non-sigs (n = %i)',sum(iSort_non)))

subplot(7,4,[16]+iRow); % movSpk_traces
h2_1 = boundedline(t2p,mean(movSpk_lobv(:,sigBoolMov_wilc_lobv),2), std(movSpk_lobv(:,sigBoolMov_wilc_lobv),[],2), '-k','alpha'); set(h2_1, 'linewidth', 1.5)
h2_2 = boundedline(t2p,mean(movSpk_lobv(:,~sigBoolMov_wilc_lobv),2), std(movSpk_lobv(:,~sigBoolMov_wilc_lobv),[],2), '-b','alpha'); set(h2_2, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 8]); yticks([-2:2:8])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel ('Event rate (spk/s)')

subplot(7,4,[20]+iRow); % movLick_traces
h2_3 = boundedline(tpaq,mean(movLick_lobv,2), std(movLick_lobv,[],2), '-k','alpha'); set(h2_3, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -.1 1.1]); yticks([0:.25:1])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel('Licking')

subplot(7,4,[24]+iRow); % movMov_traces
h2_4 = boundedline(tpaq,mean(movMov_lobv,2), std(movMov_lobv,[],2), '-k','alpha'); set(h2_4, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 4]); yticks([-2:1:4]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlabel('Time from event (s)'); ylabel('Wheel (cm/s)')

% Rewards:
iRow = 3;
coeff = pca(rewSpk_odd_lobv(end/2-imfs/2+1:end/2+imfs/2,:));
[~, iSort] = sort(coeff(:,1),'descend');
iSort_sig = sigBoolRew_wilc_lobv(iSort);
iSort_non = ~sigBoolRew_wilc_lobv(iSort);
rewSpk_even_lobv2 = rewSpk_even_lobv(:,iSort);

subplot(7,4,[0,4]+iRow); % rewSpk_sig
spklocal_sig = rewSpk_even_lobv2(:,iSort_sig);
imagesc(t2p,[1:size(spklocal_sig,2)],spklocal_sig',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Sigs (n = %i)',sum(iSort_sig)))
title('Trial reward')

subplot(7,4,[8,12]+iRow); % rewSpk_non
spklocal_non = rewSpk_even_lobv2(:,iSort_non);
imagesc(t2p,[1:size(spklocal_non,2)],spklocal_non',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Non-sigs (n = %i)',sum(iSort_non)))

subplot(7,4,[16]+iRow); % rewSpk_traces
h3_1 = boundedline(t2p,mean(rewSpk_lobv(:,sigBoolRew_wilc_lobv),2), std(rewSpk_lobv(:,sigBoolRew_wilc_lobv),[],2), '-k','alpha'); set(h3_1, 'linewidth', 1.5)
h3_2 = boundedline(t2p,mean(rewSpk_lobv(:,~sigBoolRew_wilc_lobv),2), std(rewSpk_lobv(:,~sigBoolRew_wilc_lobv),[],2), '-b','alpha'); set(h3_2, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 8]); yticks([-2:2:8])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel ('Event rate (spk/s)')

subplot(7,4,[20]+iRow); % rewLick_traces
h3_3 = boundedline(tpaq,mean(rewLick_lobv,2), std(rewLick_lobv,[],2), '-k','alpha'); set(h3_3, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -.1 1.1]); yticks([0:.25:1])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel('Licking')

subplot(7,4,[24]+iRow); % rewMov_traces
h3_4 = boundedline(tpaq,mean(rewMov_lobv,2), std(rewMov_lobv,[],2), '-k','alpha'); set(h3_4, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 4]); yticks([-2:1:4]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlabel('Time from event (s)'); ylabel('Wheel (cm/s)')

% Licking:
iRow = 4;
coeff = pca(lickSpk_odd_lobv(end/2-imfs/2+1:end/2+imfs/2,:));
[~, iSort] = sort(coeff(:,1),'descend');
iSort_sig = sigBoolLick_wilc_lobv(iSort);
iSort_non = ~sigBoolLick_wilc_lobv(iSort);
lickSpk_even_lobv2 = lickSpk_even_lobv(:,iSort);

coeff_sig = pca(lickSpk_lobv(end/2-imfs/2+1:end/2+imfs/2,sigBoolLick_wilc_lobv));
coeff_non = pca(lickSpk_lobv(end/2-imfs/2+1:end/2+imfs/2,~sigBoolLick_wilc_lobv));
[~, iSort_sig] = sort(coeff_sig(:,1),'descend');
[~, iSort_non] = sort(coeff_non(:,1),'descend');

subplot(7,4,[0,4]+iRow); % lickSpk_sig
spklocal_sig = lickSpk_even_lobv2(:,iSort_sig);
imagesc(t2p,[1:size(spklocal_sig,2)],spklocal_sig(:,iSort_sig)',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Sigs (n = %i)',size(spklocal_sig,2)))
title('Licking')

subplot(7,4,[8,12]+iRow); % lickSpk_non
spklocal_non = lickSpk_even_lobv2(:,iSort_non);
imagesc(t2p,[1:size(spklocal_non,2)],spklocal_non(:,iSort_non)',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Non-sigs (n = %i)',size(spklocal_non,2)))

subplot(7,4,[16]+iRow); % lickSpk_traces
h4_1 = boundedline(t2p,mean(lickSpk_lobv(:,sigBoolLick_wilc_lobv),2), std(lickSpk_lobv(:,sigBoolLick_wilc_lobv),[],2), '-k','alpha'); set(h4_1, 'linewidth', 1.5)
h4_2 = boundedline(t2p,mean(lickSpk_lobv(:,~sigBoolLick_wilc_lobv),2), std(lickSpk_lobv(:,~sigBoolLick_wilc_lobv),[],2), '-b','alpha'); set(h4_2, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 8]); yticks([-2:2:8])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel ('Event rate (spk/s)')

subplot(7,4,[20]+iRow); % lickLick_traces
h4_3 = boundedline(tpaq,mean(lickLick_lobv,2), std(lickLick_lobv,[],2), '-k','alpha'); set(h4_3, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -.1 1.1]); yticks([0:.25:1])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel('Licking')

subplot(7,4,[24]+iRow); % lickMov_traces
h4_4 = boundedline(tpaq,mean(lickMov_lobv,2), std(lickMov_lobv,[],2), '-k','alpha'); set(h4_4, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 4]); yticks([-2:1:4]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlabel('Time from event (s)'); ylabel('Wheel (cm/s)')

suptitle('Lobule V - Expert receptive fields (Wilcoxon)');

%% Sim2 pop data

f3a = figure;

sigBoolVis_glm_sim2 = cat(1,seshv4rfs(cat(1,glmSummary.iSim2)).sigBoolVis);
sigBoolMov_glm_sim2 = cat(1,seshv4rfs(cat(1,glmSummary.iSim2)).sigBoolMov);
sigBoolRew_glm_sim2 = cat(1,seshv4rfs(cat(1,glmSummary.iSim2)).sigBoolRew);
sigBoolLick_glm_sim2 = cat(1,seshv4rfs(cat(1,glmSummary.iSim2)).sigBoolLick);

% Trial onsets:
iRow = 1;
coeff = pca(onSpk_odd_sim2(end/2-imfs/2+1:end/2+imfs/2,:));
[~, iSort] = sort(coeff(:,1),'descend');
iSort_sig = sigBoolVis_glm_sim2(iSort);
iSort_non = ~sigBoolVis_glm_sim2(iSort);
onSpk_even_sim2_2 = onSpk_even_sim2(:,iSort);

subplot(7,4,[0,4]+iRow); % onSpk_sig
spklocal_sig = onSpk_even_sim2_2(:,iSort_sig);
imagesc(t2p,[1:size(spklocal_sig,2)],spklocal_sig',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Sigs (n = %i)',sum(iSort_sig)))
title('Trial onset')

subplot(7,4,[8,12]+iRow); % onSpk_non
spklocal_non = onSpk_even_sim2_2(:,iSort_non);
imagesc(t2p,[1:size(spklocal_non,2)],spklocal_non',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Non-sigs (n = %i)',sum(iSort_non)));

subplot(7,4,[16]+iRow); % onSpk_traces
h1_1 = boundedline(t2p,mean(onSpk_sim2_sig,2), std(onSpk_sim2_sig,[],2), '-k','alpha'); set(h1_1, 'linewidth', 1.5)
h1_2 = boundedline(t2p,mean(onSpk_sim2_non,2), std(onSpk_sim2_non,[],2), '-b','alpha'); set(h1_2, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 8]); yticks([-2:2:8])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel ('Event rate (spk/s)')

subplot(7,4,[20]+iRow); % onLick_traces
h1_3 = boundedline(tpaq,mean(onLick_sim2,2), std(onLick_sim2,[],2), '-k','alpha'); set(h1_3, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -.1 1.1]); yticks([0:.25:1])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel('Licking')

subplot(7,4,[24]+iRow); % onMov_traces
h1_4 = boundedline(tpaq,mean(onMov_sim2,2), std(onMov_sim2,[],2), '-k','alpha'); set(h1_4, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 4]); yticks([-2:1:4]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlabel('Time from event (s)'); ylabel('Wheel (cm/s)')

% Wheel movements:
iRow = 2;
coeff = pca(movSpk_odd_sim2(end/2-imfs/2+1:end/2+imfs/2,:));
[~, iSort] = sort(coeff(:,1),'descend');
iSort_sig = sigBoolMov_glm_sim2(iSort);
iSort_non = ~sigBoolMov_glm_sim2(iSort);
movSpk_even_sim2_2 = movSpk_even_sim2(:,iSort);

subplot(7,4,[0,4]+iRow); % movSpk_sig
spklocal_sig = movSpk_even_sim2_2(:,iSort_sig);
imagesc(t2p,[1:size(spklocal_sig,2)],spklocal_sig',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Sigs (n = %i)',size(spklocal_sig,2)))
title('Wheel movement')

subplot(7,4,[8,12]+iRow); % movSpk_non
spklocal_non = movSpk_even_sim2_2(:,iSort_non);
imagesc(t2p,[1:size(spklocal_non,2)],spklocal_non',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Non-sigs (n = %i)',size(spklocal_non,2)))

subplot(7,4,[16]+iRow); % movSpk_traces
h2_1 = boundedline(t2p,mean(movSpk_sim2_sig,2), std(movSpk_sim2_sig,[],2), '-k','alpha'); set(h2_1, 'linewidth', 1.5)
h2_2 = boundedline(t2p,mean(movSpk_sim2_non,2), std(movSpk_sim2_non,[],2), '-b','alpha'); set(h2_2, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 8]); yticks([-2:2:8])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel ('Event rate (spk/s)')

subplot(7,4,[20]+iRow); % movLick_traces
h2_3 = boundedline(tpaq,mean(movLick_sim2,2), std(movLick_sim2,[],2), '-k','alpha'); set(h2_3, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -.1 1.1]); yticks([0:.25:1])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel('Licking')

subplot(7,4,[24]+iRow); % movMov_traces
h2_4 = boundedline(tpaq,mean(movMov_sim2,2), std(movMov_sim2,[],2), '-k','alpha'); set(h2_4, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 4]); yticks([-2:1:4]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlabel('Time from event (s)'); ylabel('Wheel (cm/s)')

% Rewards:
iRow = 3;
coeff = pca(rewSpk_odd_sim2(end/2-imfs/2+1:end/2+imfs/2,:));
[~, iSort] = sort(coeff(:,1),'descend');
iSort_sig = sigBoolRew_glm_sim2(iSort);
iSort_non = ~sigBoolRew_glm_sim2(iSort);
rewSpk_even_sim2_2 = rewSpk_even_sim2(:,iSort);

subplot(7,4,[0,4]+iRow); % rewSpk_sig
spklocal_sig = rewSpk_even_sim2_2(:,iSort_sig);
imagesc(t2p,[1:size(spklocal_sig,2)],spklocal_sig',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Sigs (n = %i)',size(spklocal_sig,2)))
title('Trial reward')

subplot(7,4,[8,12]+iRow); % rewSpk_non
spklocal_non = rewSpk_even_sim2_2(:,iSort_non);
imagesc(t2p,[1:size(spklocal_non,2)],spklocal_non',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Non-sigs (n = %i)',size(spklocal_non,2)))

subplot(7,4,[16]+iRow); % rewSpk_traces
h3_1 = boundedline(t2p,mean(rewSpk_sim2_sig,2), std(rewSpk_sim2_sig,[],2), '-k','alpha'); set(h3_1, 'linewidth', 1.5)
h3_2 = boundedline(t2p,mean(rewSpk_sim2_non,2), std(rewSpk_sim2_non,[],2), '-b','alpha'); set(h3_2, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 16]); yticks([-2:2:16])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel ('Event rate (spk/s)')

subplot(7,4,[20]+iRow); % rewLick_traces
h3_3 = boundedline(tpaq,mean(rewLick_sim2,2), std(rewLick_sim2,[],2), '-k','alpha'); set(h3_3, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -.1 1.1]); yticks([0:.25:1])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel('Licking')

subplot(7,4,[24]+iRow); % rewMov_traces
h3_4 = boundedline(tpaq,mean(rewMov_sim2,2), std(rewMov_sim2,[],2), '-k','alpha'); set(h3_4, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 4]); yticks([-2:1:4]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlabel('Time from event (s)'); ylabel('Wheel (cm/s)')

% Licking:
iRow = 4;
coeff = pca(lickSpk_odd_sim2(end/2-imfs/2+1:end/2+imfs/2,:));
[~, iSort] = sort(coeff(:,1),'descend');
iSort_sig = sigBoolLick_glm_sim2(iSort);
iSort_non = ~sigBoolLick_glm_sim2(iSort);
lickSpk_even_sim2_2 = lickSpk_even_sim2(:,iSort);

subplot(7,4,[0,4]+iRow); % lickSpk_sig
spklocal_sig = lickSpk_even_sim2_2(:,iSort_sig);
imagesc(t2p,[1:size(spklocal_sig,2)],spklocal_sig',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Sigs (n = %i)',size(spklocal_sig,2)))
title('Licking')

subplot(7,4,[8,12]+iRow); % movSpk_non
spklocal_non = lickSpk_even_sim2_2(:,iSort_non);
imagesc(t2p,[1:size(spklocal_non,2)],spklocal_non',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Non-sigs (n = %i)',size(spklocal_non,2)))

subplot(7,4,[16]+iRow); % lickSpk_traces
h4_1 = boundedline(t2p,mean(lickSpk_sim2_sig,2), std(lickSpk_sim2_sig,[],2), '-k','alpha'); set(h4_1, 'linewidth', 1.5)
h4_2 = boundedline(t2p,mean(lickSpk_sim2_non,2), std(lickSpk_sim2_non,[],2), '-b','alpha'); set(h4_2, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 8]); yticks([-2:2:8])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel ('Event rate (spk/s)')

subplot(7,4,[20]+iRow); % lickLick_traces
h4_3 = boundedline(tpaq,mean(lickLick_sim2,2), std(lickLick_sim2,[],2), '-k','alpha'); set(h4_3, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -.1 1.1]); yticks([0:.25:1])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel('Licking')

subplot(7,4,[24]+iRow); % lickMov_traces
h4_4 = boundedline(tpaq,mean(lickMov_sim2,2), std(lickMov_sim2,[],2), '-k','alpha'); set(h4_4, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 4]); yticks([-2:1:4]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlabel('Time from event (s)'); ylabel('Wheel (cm/s)')

suptitle('Lobule simplex - Expert receptive fields (GLM)')

%% Sim2 sorting - Wilcoxon

f3b = figure;

sigBoolVis_wilc_sim2 = cat(1,seshv4rfs(cat(1,glmSummary.iSim2)).sigBoolVis_wilc);
sigBoolMov_wilc_sim2 = cat(1,seshv4rfs(cat(1,glmSummary.iSim2)).sigBoolMov_wilc);
sigBoolRew_wilc_sim2 = cat(1,seshv4rfs(cat(1,glmSummary.iSim2)).sigBoolRew_wilc);
sigBoolLick_wilc_sim2 = cat(1,seshv4rfs(cat(1,glmSummary.iSim2)).sigBoolLick_wilc);


% Trial onsets:
iRow = 1;
coeff = pca(onSpk_odd_sim2(end/2-imfs/2+1:end/2+imfs/2,:));
[~, iSort] = sort(coeff(:,1),'descend');
iSort_sig = sigBoolVis_wilc_sim2(iSort);
iSort_non = ~sigBoolVis_wilc_sim2(iSort);
onSpk_even_sim2_2 = onSpk_even_sim2(:,iSort);

subplot(7,4,[0,4]+iRow); % onSpk_sig
spklocal_sig = onSpk_even_sim2_2(:,iSort_sig);
imagesc(t2p,[1:size(spklocal_sig,2)],spklocal_sig',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Sigs (n = %i)',sum(iSort_sig)))
title('Trial onset')

subplot(7,4,[8,12]+iRow); % onSpk_non
spklocal_non = onSpk_even_sim2_2(:,iSort_non);
imagesc(t2p,[1:size(spklocal_non,2)],spklocal_non',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Non-sigs (n = %i)',sum(iSort_non)))

subplot(7,4,[16]+iRow); % onSpk_traces
h1_1 = boundedline(t2p,mean(onSpk_sim2(:,sigBoolVis_wilc_sim2),2), std(onSpk_sim2(:,sigBoolVis_wilc_sim2),[],2), '-k','alpha'); set(h1_1, 'linewidth', 1.5)
h1_2 = boundedline(t2p,mean(onSpk_sim2(:,~sigBoolVis_wilc_sim2),2), std(onSpk_sim2(:,~sigBoolVis_wilc_sim2),[],2), '-b','alpha'); set(h1_2, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 8]); yticks([-2:2:8])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel ('Event rate (spk/s)')

subplot(7,4,[20]+iRow); % onLick_traces
h1_3 = boundedline(tpaq,mean(onLick_sim2,2), std(onLick_sim2,[],2), '-k','alpha'); set(h1_3, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -.1 1.1]); yticks([0:.25:1])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel('Licking')

subplot(7,4,[24]+iRow); % onMov_traces
h1_4 = boundedline(tpaq,mean(onMov_sim2,2), std(onMov_sim2,[],2), '-k','alpha'); set(h1_4, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 4]); yticks([-2:1:4]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlabel('Time from event (s)'); ylabel('Wheel (cm/s)')

% Wheel movements:
iRow = 2;
coeff = pca(movSpk_odd_sim2(end/2-imfs/2+1:end/2+imfs/2,:));
[~, iSort] = sort(coeff(:,1),'descend');
iSort_sig = sigBoolMov_wilc_sim2(iSort);
iSort_non = ~sigBoolMov_wilc_sim2(iSort);
movSpk_even_sim2_2 = movSpk_even_sim2(:,iSort);

subplot(7,4,[0,4]+iRow); % movSpk_sig
spklocal_sig = movSpk_even_sim2_2(:,iSort_sig);
imagesc(t2p,[1:size(spklocal_sig,2)],spklocal_sig',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Sigs (n = %i)',sum(iSort_sig)))
title('Wheel movement')

subplot(7,4,[8,12]+iRow); % movSpk_non
spklocal_non = movSpk_even_sim2_2(:,iSort_non);
imagesc(t2p,[1:size(spklocal_non,2)],spklocal_non',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Non-sigs (n = %i)',sum(iSort_non)))

subplot(7,4,[16]+iRow); % movSpk_traces
h2_1 = boundedline(t2p,mean(movSpk_sim2(:,sigBoolMov_wilc_sim2),2), std(movSpk_sim2(:,sigBoolMov_wilc_sim2),[],2), '-k','alpha'); set(h2_1, 'linewidth', 1.5)
h2_2 = boundedline(t2p,mean(movSpk_sim2(:,~sigBoolMov_wilc_sim2),2), std(movSpk_sim2(:,~sigBoolMov_wilc_sim2),[],2), '-b','alpha'); set(h2_2, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 8]); yticks([-2:2:8])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel ('Event rate (spk/s)')

subplot(7,4,[20]+iRow); % movLick_traces
h2_3 = boundedline(tpaq,mean(movLick_sim2,2), std(movLick_sim2,[],2), '-k','alpha'); set(h2_3, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -.1 1.1]); yticks([0:.25:1])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel('Licking')

subplot(7,4,[24]+iRow); % movMov_traces
h2_4 = boundedline(tpaq,mean(movMov_sim2,2), std(movMov_sim2,[],2), '-k','alpha'); set(h2_4, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 4]); yticks([-2:1:4]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlabel('Time from event (s)'); ylabel('Wheel (cm/s)')

% Rewards:
iRow = 3;
coeff = pca(rewSpk_odd_sim2(end/2-imfs/2+1:end/2+imfs/2,:));
[~, iSort] = sort(coeff(:,1),'descend');
iSort_sig = sigBoolRew_wilc_sim2(iSort);
iSort_non = ~sigBoolRew_wilc_sim2(iSort);
rewSpk_even_sim2_2 = rewSpk_even_sim2(:,iSort);

subplot(7,4,[0,4]+iRow); % rewSpk_sig
spklocal_sig = rewSpk_even_sim2_2(:,iSort_sig);
imagesc(t2p,[1:size(spklocal_sig,2)],spklocal_sig',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Sigs (n = %i)',sum(iSort_sig)))
title('Trial reward')

subplot(7,4,[8,12]+iRow); % rewSpk_non
spklocal_non = rewSpk_even_sim2_2(:,iSort_non);
imagesc(t2p,[1:size(spklocal_non,2)],spklocal_non',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Non-sigs (n = %i)',sum(iSort_non)))

subplot(7,4,[16]+iRow); % rewSpk_traces
h3_1 = boundedline(t2p,mean(rewSpk_sim2(:,sigBoolRew_wilc_sim2),2), std(rewSpk_sim2(:,sigBoolRew_wilc_sim2),[],2), '-k','alpha'); set(h3_1, 'linewidth', 1.5)
h3_2 = boundedline(t2p,mean(rewSpk_sim2(:,~sigBoolRew_wilc_sim2),2), std(rewSpk_sim2(:,~sigBoolRew_wilc_sim2),[],2), '-b','alpha'); set(h3_2, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 8]); yticks([-2:2:8])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel ('Event rate (spk/s)')

subplot(7,4,[20]+iRow); % rewLick_traces
h3_3 = boundedline(tpaq,mean(rewLick_sim2,2), std(rewLick_sim2,[],2), '-k','alpha'); set(h3_3, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -.1 1.1]); yticks([0:.25:1])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel('Licking')

subplot(7,4,[24]+iRow); % rewMov_traces
h3_4 = boundedline(tpaq,mean(rewMov_sim2,2), std(rewMov_sim2,[],2), '-k','alpha'); set(h3_4, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 4]); yticks([-2:1:4]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlabel('Time from event (s)'); ylabel('Wheel (cm/s)')

% Licking:
iRow = 4;
coeff = pca(lickSpk_odd_sim2(end/2-imfs/2+1:end/2+imfs/2,:));
[~, iSort] = sort(coeff(:,1),'descend');
iSort_sig = sigBoolLick_wilc_sim2(iSort);
iSort_non = ~sigBoolLick_wilc_sim2(iSort);
lickSpk_even_sim2_2 = lickSpk_even_sim2(:,iSort);

coeff_sig = pca(lickSpk_sim2(end/2-imfs/2+1:end/2+imfs/2,sigBoolLick_wilc_sim2));
coeff_non = pca(lickSpk_sim2(end/2-imfs/2+1:end/2+imfs/2,~sigBoolLick_wilc_sim2));
[~, iSort_sig] = sort(coeff_sig(:,1),'descend');
[~, iSort_non] = sort(coeff_non(:,1),'descend');

subplot(7,4,[0,4]+iRow); % lickSpk_sig
spklocal_sig = lickSpk_even_sim2_2(:,iSort_sig);
imagesc(t2p,[1:size(spklocal_sig,2)],spklocal_sig(:,iSort_sig)',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Sigs (n = %i)',size(spklocal_sig,2)))
title('Licking')

subplot(7,4,[8,12]+iRow); % lickSpk_non
spklocal_non = lickSpk_even_sim2_2(:,iSort_non);
imagesc(t2p,[1:size(spklocal_non,2)],spklocal_non(:,iSort_non)',[0 8]); xlim([-1 1]);
xline(0,'--w', 'LineWidth',1)
colormap viridis; box off;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel(sprintf('Non-sigs (n = %i)',size(spklocal_non,2)))

subplot(7,4,[16]+iRow); % lickSpk_traces
h4_1 = boundedline(t2p,mean(lickSpk_sim2(:,sigBoolLick_wilc_sim2),2), std(lickSpk_sim2(:,sigBoolLick_wilc_sim2),[],2), '-k','alpha'); set(h4_1, 'linewidth', 1.5)
h4_2 = boundedline(t2p,mean(lickSpk_sim2(:,~sigBoolLick_wilc_sim2),2), std(lickSpk_sim2(:,~sigBoolLick_wilc_sim2),[],2), '-b','alpha'); set(h4_2, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 8]); yticks([-2:2:8])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel ('Event rate (spk/s)')

subplot(7,4,[20]+iRow); % lickLick_traces
h4_3 = boundedline(tpaq,mean(lickLick_sim2,2), std(lickLick_sim2,[],2), '-k','alpha'); set(h4_3, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -.1 1.1]); yticks([0:.25:1])
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
ylabel('Licking')

subplot(7,4,[24]+iRow); % lickMov_traces
h4_4 = boundedline(tpaq,mean(lickMov_sim2,2), std(lickMov_sim2,[],2), '-k','alpha'); set(h4_4, 'linewidth', 1.5)
xline(0,'--k', 'LineWidth',1)
axis([-1 1 -2 4]); yticks([-2:1:4]);
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlabel('Time from event (s)'); ylabel('Wheel (cm/s)')

suptitle('Lobule simplex - Expert receptive fields (Wilcoxon)');


%% Make confusion matrices for different events:
onSig_lobv_glm = cat(1,seshv4rfs(cat(1,glmSummary.iLobV)).sigBoolVis);
onSig_lobv_wilc = cat(1,seshv4rfs(cat(1,glmSummary.iLobV)).sigBoolVis_wilc);
movSig_lobv_glm = cat(1,seshv4rfs(cat(1,glmSummary.iLobV)).sigBoolMov);
movSig_lobv_wilc = cat(1,seshv4rfs(cat(1,glmSummary.iLobV)).sigBoolMov_wilc);
rewSig_lobv_glm = cat(1,seshv4rfs(cat(1,glmSummary.iLobV)).sigBoolRew);
rewSig_lobv_wilc = cat(1,seshv4rfs(cat(1,glmSummary.iLobV)).sigBoolRew_wilc);
lickSig_lobv_glm = cat(1,seshv4rfs(cat(1,glmSummary.iLobV)).sigBoolLick);
lickSig_lobv_wilc = cat(1,seshv4rfs(cat(1,glmSummary.iLobV)).sigBoolLick_wilc);

onSig_sim2_glm = cat(1,seshv4rfs(cat(1,glmSummary.iSim2)).sigBoolVis);
onSig_sim2_wilc = cat(1,seshv4rfs(cat(1,glmSummary.iSim2)).sigBoolVis_wilc);
movSig_sim2_glm = cat(1,seshv4rfs(cat(1,glmSummary.iSim2)).sigBoolMov);
movSig_sim2_wilc = cat(1,seshv4rfs(cat(1,glmSummary.iSim2)).sigBoolMov_wilc);
rewSig_sim2_glm = cat(1,seshv4rfs(cat(1,glmSummary.iSim2)).sigBoolRew);
rewSig_sim2_wilc = cat(1,seshv4rfs(cat(1,glmSummary.iSim2)).sigBoolRew_wilc);
lickSig_sim2_glm = cat(1,seshv4rfs(cat(1,glmSummary.iSim2)).sigBoolLick);
lickSig_sim2_wilc = cat(1,seshv4rfs(cat(1,glmSummary.iSim2)).sigBoolLick_wilc);

confMat_lobv_on = confusionmat(onSig_lobv_glm, onSig_lobv_wilc,'Order',[1 0]) ./ [sum(onSig_lobv_glm) sum(~onSig_lobv_glm)]'*100;
confMat_lobv_mov = confusionmat(movSig_lobv_glm, movSig_lobv_wilc,'Order',[1 0]) ./ [sum(movSig_lobv_glm) sum(~movSig_lobv_glm)]'*100;
confMat_lobv_rew = confusionmat(rewSig_lobv_glm, rewSig_lobv_wilc,'Order',[1 0]) ./ [sum(rewSig_lobv_glm) sum(~rewSig_lobv_glm)]'*100;
confMat_lobv_lick = confusionmat(lickSig_lobv_glm, lickSig_lobv_wilc,'Order',[1 0]) ./ [sum(lickSig_lobv_glm) sum(~lickSig_lobv_glm)]'*100;

confMat_sim2_on = confusionmat(onSig_sim2_glm, onSig_sim2_wilc,'Order',[1 0]) ./ [sum(onSig_sim2_glm) sum(~onSig_sim2_glm)]'*100;
confMat_sim2_mov = confusionmat(movSig_sim2_glm, movSig_sim2_wilc,'Order',[1 0]) ./ [sum(movSig_sim2_glm) sum(~movSig_sim2_glm)]'*100;
confMat_sim2_rew = confusionmat(rewSig_sim2_glm, rewSig_sim2_wilc,'Order',[1 0]) ./ [sum(rewSig_sim2_glm) sum(~rewSig_sim2_glm)]'*100;
confMat_sim2_lick = confusionmat(vv, lickSig_sim2_wilc,'Order',[1 0]) ./ [sum(lickSig_sim2_glm) sum(~lickSig_sim2_glm)]'*100;

f4 = figure;
subplot(2,4,1);
imagesc(confMat_lobv_on,[0, 100]);
axis square; colormap gray; box off; colorbar
xticks([1,2]); yticks([1,2]); xticklabels({'Sig','Non-sig'}); yticklabels({'Sig','Non-sig'});
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlabel('Wilcoxon class'); ylabel('GLM class'); title('Lobv: on')
subplot(2,4,2);
imagesc(confMat_lobv_mov,[0, 100]);
axis square; colormap gray; box off; colorbar
xticks([1,2]); yticks([1,2]); xticklabels({'Sig','Non-sig'}); yticklabels({'Sig','Non-sig'});
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlabel('Wilcoxon class'); ylabel('GLM class'); title('Lobv: mov')
subplot(2,4,3);
imagesc(confMat_lobv_rew,[0, 100]);
axis square; colormap gray; box off; colorbar
xticks([1,2]); yticks([1,2]); xticklabels({'Sig','Non-sig'}); yticklabels({'Sig','Non-sig'});
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlabel('Wilcoxon class'); ylabel('GLM class'); title('Lobv: rew')
subplot(2,4,4);
imagesc(confMat_lobv_lick,[0, 100]);
axis square; colormap gray; box off; colorbar
xticks([1,2]); yticks([1,2]); xticklabels({'Sig','Non-sig'}); yticklabels({'Sig','Non-sig'});
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlabel('Wilcoxon class'); ylabel('GLM class'); title('Lobv: lick')
subplot(2,4,5);
imagesc(confMat_sim2_on,[0, 100]);
axis square; colormap gray; box off; colorbar
xticks([1,2]); yticks([1,2]); xticklabels({'Sig','Non-sig'}); yticklabels({'Sig','Non-sig'});
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlabel('Wilcoxon class'); ylabel('GLM class'); title('Sim2: on')
subplot(2,4,6);
imagesc(confMat_sim2_mov,[0, 100]);
axis square; colormap gray; box off; colorbar
xticks([1,2]); yticks([1,2]); xticklabels({'Sig','Non-sig'}); yticklabels({'Sig','Non-sig'});
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlabel('Wilcoxon class'); ylabel('GLM class'); title('Sim2: mov')
subplot(2,4,7);
imagesc(confMat_sim2_rew,[0, 100]);
axis square; colormap gray; box off; colorbar
xticks([1,2]); yticks([1,2]); xticklabels({'Sig','Non-sig'}); yticklabels({'Sig','Non-sig'});
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlabel('Wilcoxon class'); ylabel('GLM class'); title('Sim2: rew')
subplot(2,4,8);
imagesc(confMat_sim2_lick,[0, 100]);
axis square; colormap gray; box off; colorbar
xticks([1,2]); yticks([1,2]); xticklabels({'Sig','Non-sig'}); yticklabels({'Sig','Non-sig'});
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlabel('Wilcoxon class'); ylabel('GLM class'); title('Sim2: lick')

%% save maybe

if savebool

    figFold = fullfile(paperFold,'Figures','Fig3','SuppFig_receptiveFields');
    if ~exist(figFold,'dir'); mkdir(figFold); end
    figname = 'lobvAll_expertRFs_GLM';
    savefig(f2a, fullfile(figFold,[figname,'.fig']));
    saveas(f2a,fullfile(figFold,[figname,'.png']));
    print(f2a,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');
    figname = 'sim2All_expertRFs_GLM';
    savefig(f3a, fullfile(figFold,[figname,'.fig']));
    saveas(f3a,fullfile(figFold,[figname,'.png']));
    print(f3a,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');
    figname = 'lobvAll_expertRFs_wilc';
    savefig(f2b, fullfile(figFold,[figname,'.fig']));
    saveas(f2b,fullfile(figFold,[figname,'.png']));
    print(f2b,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');
    figname = 'sim2All_expertRFs_wilc';
    savefig(f3b, fullfile(figFold,[figname,'.fig']));
    saveas(f3b,fullfile(figFold,[figname,'.png']));
    print(f3b,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');
    figname = 'rf_confusionMat';
    savefig(f4, fullfile(figFold,[figname,'.fig']));
    saveas(f4,fullfile(figFold,[figname,'.png']));
    print(f4,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');

    save(fullfile(figFold,'v4_receptiveFields.mat'),'seshv4rfs','lobvSummary','sim2Summary');    

end