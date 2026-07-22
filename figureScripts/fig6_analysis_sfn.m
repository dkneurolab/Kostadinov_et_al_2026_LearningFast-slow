function fig6_analysis_sfn(matchData, paperFold, savebool)

% % Big matrix of task data in matched rois
% load('/Users/dimitar/Desktop/imData/FastSlowLearning/matchDataBehav.mat')
imfs = 30;
% Load v4 significant GLM summary:
sigSummary = load(fullfile(paperFold,'Figures/Fig3/Sub3_SummaryBars/spk_gaussian_alpha50/sigSummary.mat'));
glmSummary = sigSummary.glmSummary;

% Indices of 145 matched ROIs:
v145Matched = load(fullfile(paperFold,'Figures/Fig1/Sub4_summary/matchPCs_v145.mat'));
v145Matched = cat(2,v145Matched.lobvPCs,v145Matched.sim2PCs);

% v144 kernels:
GLMv144 = load(fullfile(paperFold,'Figures/Fig4/Sub5_fitKernels/kernelsGLMv144.mat'));
GLMv144 = GLMv144.GLMv144;

% v5 adaptation behaviour
v5Behav = load(fullfile(paperFold,'Figures/Fig5/Sub2_behavSummary/v5_adaptSummary.mat'));
for i = 1:numel(v5Behav.lobvdata); v5Behav.lobvdata(i).fov = 'lobv'; end
for i = 1:numel(v5Behav.sim2data); v5Behav.sim2data(i).fov = 'sim2'; end
v5Behav = cat(2,v5Behav.lobvdata,v5Behav.sim2data); clear i

%% Make master data struct that takes in GLMv144 and matched ROI indices and makes a single data struct
% Make separate structs for mov and rew, since the roi indices may appear
% in both and we want to avoid confusion
v145movdata = struct('name',[],'fov',[],'iv1',[],'iv4',[],'iv5',[], ...
    'v1movk',[],'v4movk',[],'v14movcorr',[],'v5movBase',[],'v5movMid',[],'v5movRatio',[]);
v145rewdata = struct('name',[],'fov',[],'iv1',[],'iv4',[],'iv5',[],...
    'v1rewk',[],'v4rewk',[],'v14rewcorr',[],'v5rewBase',[],'v5rewMid',[],'v5rewRatio',[]);

a = 0; b = 0;
for i145 = 1:numel(v145Matched)
    % Get name and lobule of v145 matched dendrites
    mouse = v145Matched(i145).name;
    fov = v145Matched(i145).fov;
    if strcmpi(fov,'lobv')
        ilobv = true;
        isim2 = false;
    else
        ilobv = false;
        isim2 = true;
    end
    % Identify the matching row in GLMv144 structure
    jMatch = 0;
    for j144 = 1:numel(GLMv144)
        if strcmpi(GLMv144(j144).name,mouse) && ...
                GLMv144(j144).ilobv == ilobv && ...
                GLMv144(j144).isim2 == isim2
            jMatch = j144;
        end
    end
    v14local = GLMv144(jMatch);
    v145rois = v145Matched(i145).v145rois;
    % Go through v14matched rois and extract movement kernels and correlate
    if isfield(v14local.GLMkernelsv1,'kernelMov')
        for iKernel = 1:numel(v14local.GLMkernelsv1)
            if ~isempty(v14local.GLMkernelsv1(iKernel).kernelMov) && ...
                    ismember(v14local.GLMkernelsv1(iKernel).iPC_v4,v145rois(:,2))
                a = a+1;
                [~,matchRow] = ismember(v14local.GLMkernelsv1(iKernel).iPC_v4,v145rois(:,2));
                v145movdata(a).name = mouse;
                v145movdata(a).fov = fov;
                v145movdata(a).iv1 = v145rois(matchRow,1);
                v145movdata(a).iv4 = v145rois(matchRow,2);
                v145movdata(a).iv5 = v145rois(matchRow,3);
                v145movdata(a).v1movk = v14local.GLMkernelsv1(iKernel).kernelMov;
                v145movdata(a).v4movk = v14local.GLMkernelsv4n(iKernel).kernelMov;
                % Compute kernel correlation
                movCorr = corrcoef(v145movdata(a).v1movk,v145movdata(a).v4movk);
                v145movdata(a).v14movcorr = movCorr(1,2);
            end
        end

    end

    % Go through v14matched rois and extract reward kernels and correlate
    if isfield(v14local.GLMkernelsv1,'kernelRew')
        for iKernel = 1:numel(v14local.GLMkernelsv1)
            if ~isempty(v14local.GLMkernelsv1(iKernel).kernelRew) && ...
                    ismember(v14local.GLMkernelsv1(iKernel).iPC_v4,v145rois(:,2))
                b = b+1;
                [~,matchRow] = ismember(v14local.GLMkernelsv1(iKernel).iPC_v4,v145rois(:,2));
                v145rewdata(b).name = mouse;
                v145rewdata(b).fov = fov;
                v145rewdata(b).iv1 = v145rois(matchRow,1);
                v145rewdata(b).iv4 = v145rois(matchRow,2);
                v145rewdata(b).iv5 = v145rois(matchRow,3);
                v145rewdata(b).v1rewk = v14local.GLMkernelsv1(iKernel).kernelRew;
                v145rewdata(b).v4rewk = v14local.GLMkernelsv4n(iKernel).kernelRew;
                % Compute kernel correlation
                rewCorr = corrcoef(v145rewdata(b).v1rewk,v145rewdata(b).v4rewk);
                v145rewdata(b).v14rewcorr = rewCorr(1,2);

            end
        end
    end
end


clear a b fov i145 iKernel ilobv isim2 j144 jMatch matchRow mouse movCorr rewCorr v145rois v14local
%% Now populate v5 data using same ROIs!
% Go through matchData and get matching roi info using GLM sigsummary of v4
% data - issue is getting the indices of v5Behav ROIs

for i145 = 1:numel(v145Matched)
    % Get name and lobule of v145 matched dendrites
    mouse = v145Matched(i145).name;
    fov = v145Matched(i145).fov;
    % Get indices of matched ROIs to cross reference to GLM significance
    for jData = 1:numel(matchData)
        if strcmpi(matchData(jData).name, mouse) && strcmpi(matchData(jData).fov, fov) && ...
                matchData(jData).v1 == 4 && matchData(jData).v2 == 5
            iv4 = matchData(jData).iROIs1;
            iv5 = matchData(jData).iROIs2;
        end
    end
    % Get logical significance for v4 GLM fits
    for kglm = 1:numel(glmSummary)
        if strcmpi(glmSummary(kglm).name,mouse) && strcmpi(glmSummary(kglm).fov,fov)
            movSig = glmSummary(kglm).sigBoolMov;
            rewSig = glmSummary(kglm).sigBoolRew;
        end
    end
    % Use ismember to get matching ROIs - should match numbers in v5Behav!
    iMov = find(movSig);
    iRew = find(rewSig);
    [~,bb] = ismember(iMov,iv4);
    bb = bb(bb > 0);
    iv5mov = iv5(bb);
    [~,bb] = ismember(iRew,iv4);
    bb = bb(bb > 0);
    iv5rew = iv5(bb);
    
    % Now get indices of subset of iv5mov and iv5rew that were matched
    v145rois = v145Matched(i145).v145rois;
    [iv5movMatched,cc] = ismember(iv5mov,v145rois(:,3));
    [iv5rewMatched,dd] = ismember(iv5rew,v145rois(:,3));
    cc = cc(cc > 0);
    dd = dd(dd > 0);
    iv5movGood = v145rois(cc,3);
    iv5rewGood = v145rois(dd,3);
    
    % Find corresponding v5 session and extract ratio of mid-adapt to
    % baseline signal
    for iv5 = 1:numel(v5Behav)
        if strcmpi(v5Behav(iv5).name,mouse) && strcmpi(v5Behav(iv5).fov,fov)
            iMatch = iv5;
        end
    end
    % spkMovBase = max(v5Behav(iMatch).spkMovBase(end/2-.3*imfs+1:end/2,iv5movMatched))';
    % spkMovMid = max(v5Behav(iMatch).spkMovMid(end/2-.3*imfs+1:end/2,iv5movMatched))';
    spkMovBase = mean(v5Behav(iMatch).spkMovBase(end/2-.3*imfs+1:end/2,iv5movMatched))';
    spkMovMid = mean(v5Behav(iMatch).spkMovMid(end/2-.3*imfs+1:end/2,iv5movMatched))';
    spkMovRatio = spkMovMid./spkMovBase;

    % spkRewBase = max(v5Behav(iMatch).spkRewBase(end/2+1:end/2+15,iv5rewMatched))';
    % spkRewMid = max(v5Behav(iMatch).spkRewMid(end/2+1:end/2+15,iv5rewMatched))';
    spkRewBase = mean(v5Behav(iMatch).spkRewBase(end/2+1:end/2+.2*imfs,iv5rewMatched))';
    spkRewMid = mean(v5Behav(iMatch).spkRewMid(end/2+1:end/2+.2*imfs,iv5rewMatched))';
    spkRewRatio = spkRewMid./spkRewBase;
    
    % Now populate the master matrix with the matching unit ratios
    for imovMatch = 1:numel(v145movdata)
        if strcmpi(v145movdata(imovMatch).name,mouse) && ...
                strcmpi(v145movdata(imovMatch).fov,fov) && ...
                ismember(v145movdata(imovMatch).iv5,iv5movGood)
            iRoi = ismember(iv5movGood,v145movdata(imovMatch).iv5);
            v145movdata(imovMatch).v5movBase = spkMovBase(iRoi);
            v145movdata(imovMatch).v5movMid = spkMovMid(iRoi);
            v145movdata(imovMatch).v5movRatio = spkMovRatio(iRoi);
        end
    end
    for irewMatch = 1:numel(v145rewdata)
        if strcmpi(v145rewdata(irewMatch).name,mouse) && ...
                strcmpi(v145rewdata(irewMatch).fov,fov) && ...
                ismember(v145rewdata(irewMatch).iv5,iv5rewGood)
            iRoi = ismember(iv5rewGood,v145rewdata(irewMatch).iv5);
            v145rewdata(irewMatch).v5rewBase = spkRewBase(iRoi);
            v145rewdata(irewMatch).v5rewMid = spkRewMid(iRoi);
            v145rewdata(irewMatch).v5rewRatio = spkRewRatio(iRoi);
        end
    end
end

% Now set empty values to nan, so it's easier to concatenate and plot
for imovMatch = 1:numel(v145movdata)
    if strcmpi(v145movdata(imovMatch).fov,'lobv')
        v145movdata(imovMatch).ilobv = true;
        v145movdata(imovMatch).isim2 = false;
    else
        v145movdata(imovMatch).ilobv = false;
        v145movdata(imovMatch).isim2 = true;
    end
    if isempty(v145movdata(imovMatch).v5movRatio)
        v145movdata(imovMatch).v5movBase = nan;
        v145movdata(imovMatch).v5movMid = nan;
        v145movdata(imovMatch).v5movRatio = nan;
    end
end

for irewMatch = 1:numel(v145rewdata)
    if strcmpi(v145movdata(irewMatch).fov,'lobv')
        v145rewdata(irewMatch).ilobv = true;
        v145rewdata(irewMatch).isim2 = false;
    else
        v145rewdata(irewMatch).ilobv = false;
        v145rewdata(irewMatch).isim2 = true;
    end
    if isempty(v145rewdata(irewMatch).v5rewRatio)
        v145rewdata(irewMatch).v5rewBase = nan;
        v145rewdata(irewMatch).v5rewMid = nan;
        v145rewdata(irewMatch).v5rewRatio = nan;
    end
end

clear bb cc dd fov i145 iMatch iMov movMatch iRew irewMatch iRoi iv4 iv5 iv5mov iv5movGood iv5rew iv5rewGood
clear iv5movMatched iv5rewMatched jData kglm mouse movSig rewSig spkMovBase spkMovMid spkMovRatio
clear spkRewBase spkRewMid spkRewRatio test v145rois


%% Now split by condition and plot scatters:

imovlobv = cat(1,v145movdata.ilobv);
imovsim2 = cat(1,v145movdata.isim2);
irewlobv = cat(1,v145rewdata.ilobv);
irewsim2 = cat(1,v145rewdata.isim2);

v14movcorr = cat(1,v145movdata.v14movcorr);
v14rewcorr = cat(1,v145rewdata.v14rewcorr);
v5movratio = cat(1,v145movdata.v5movRatio);
v5rewratio = cat(1,v145rewdata.v5rewRatio);

movGood = ~isnan(v14movcorr) & ~isnan(v5movratio) & ~isinf(abs(log10(v14movcorr))) & ~isinf(abs(log10(v5movratio)));
imovlobv = imovlobv(movGood);
imovsim2 = imovsim2(movGood);
v14movcorr = v14movcorr(movGood);
v5movratio = v5movratio(movGood);
rewGood = ~isnan(v14rewcorr) & ~isnan(v5rewratio) & ~isinf(abs(log10(v14rewcorr))) & ~isinf(abs(log10(v5rewratio)));
irewlobv = irewlobv(rewGood);
irewsim2 = irewsim2(rewGood);
v14rewcorr = v14rewcorr(rewGood);
v5rewratio = v5rewratio(rewGood);

f1 = figure;
subplot(1,4,1)
oLobvMov = regress_hd(v14movcorr(imovlobv),log10(v5movratio(imovlobv)),1);
boundedline(oLobvMov.x_fit,oLobvMov.y_fit,oLobvMov.delta,'cmap',[0 0 0],'alpha'); hold on;
scatter(v14movcorr(imovlobv),log10(v5movratio(imovlobv)),'filled');
xlabel('Naive-trained kernel correlation');
ylabel('Response scaling during adaptation');
title('Lobule V: movement-encoding PC dendrites')
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
% axis square; xlim([-1 1]); ylim([0 7]);
axis square; xlim([-1 1]); ylim([-1 1]);
fprintf('Number of good lobv mov PC dendrites: %i\n', sum(imovlobv))

subplot(1,4,2)
oSim2Mov = regress_hd(v14movcorr(imovsim2),log10(v5movratio(imovsim2)),1);
boundedline(oSim2Mov.x_fit,oSim2Mov.y_fit,oSim2Mov.delta,'cmap',[0 0 0],'alpha'); hold on;
scatter(v14movcorr(imovsim2),log10(v5movratio(imovsim2)),'filled');
xlabel('Naive-trained kernel correlation');
ylabel('Response scaling during adaptation');
title('Simplex: movement-encoding PC dendrites')
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
axis square; xlim([-1 1]); ylim([0 7]);
axis square; xlim([-1 1]); ylim([-1 1]);
fprintf('Number of good sim2 mov PC dendrites: %i\n', sum(imovsim2))

subplot(1,4,3)
oLobvRew = regress_hd(v14rewcorr(irewlobv),log10(v5rewratio(irewlobv)),1);
boundedline(oLobvRew.x_fit,oLobvRew.y_fit,oLobvRew.delta,'cmap',[0 0 0],'alpha'); hold on;
scatter(v14rewcorr(irewlobv),log10(v5rewratio(irewlobv)),'filled');
xlabel('Naive-trained kernel correlation');
ylabel('Response scaling during adaptation');
title('Lobule V: reward-encoding PC dendrites')
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
axis square; axis square; xlim([-1 1]); ylim([-1 1]);
fprintf('Number of good lobv rew PC dendrites: %i\n', sum(irewlobv))

subplot(1,4,4)
oSim2Rew = regress_hd(v14rewcorr(irewsim2),log10(v5rewratio(irewsim2)),1);
boundedline(oSim2Rew.x_fit,oSim2Rew.y_fit,oSim2Rew.delta,'cmap',[0 0 0],'alpha'); hold on;
scatter(v14rewcorr(irewsim2),log10(v5rewratio(irewsim2)),'filled');
xlabel('Naive-trained kernel correlation');
ylabel('Response scaling during adaptation');
title('Simplex: reward-encoding PC dendrites')
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
% axis square; axis square; xlim([-1 1]); ylim([0 7]);
axis square; axis square; xlim([-1 1]); ylim([-1 1]);
fprintf('Number of good sim2 rew PC dendrites: %i\n', sum(irewsim2))

statsOut = v2struct(oLobvMov, oLobvRew, oSim2Mov, oSim2Rew);

if savebool
    figname = 'fast_vs_slow_learning_log';
    savFold = fullfile(paperFold,'Figures/Fig6/Sub1_slow-fastScatter');
    if ~isfolder(savFold); mkdir(savFold); end
    savefig(f1, fullfile(savFold,[figname,'.fig']));
    saveas(f1,fullfile(savFold,[figname,'.png']));
    print(f1,fullfile(savFold,[figname,'.eps']), '-depsc', '-vector');
    save(fullfile(savFold,'fast_vs_slow_fits.mat'), 'statsOut')
end