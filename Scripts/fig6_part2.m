function fig6_part2(matchData, glmParams, dataFold, paperFold, savebool)
%% Settings:
% Decimation:
nDec = 50;
% Behav smoothing:
msSmthing = 60;% 60 ms smoothing span
% Remove redundant fields in matchData:
matchData = rmfield(matchData,{'on1','on2','rewC1','rewC2','rewR1','rewR2'});

% Load in v4 sigSummary
glmFold0 = fullfile(paperFold,'Figures','Fig3','Sub3_SummaryBars');
glmFold = fullfile(glmFold0,sprintf('%s_%s_alpha%i',glmParams.dataType,glmParams.distr,glmParams.alpha*100));
glmSummary = load(fullfile(glmFold,'sigSummary.mat')); glmSummary = glmSummary.glmSummary;

% 
pathv5Behav = fullfile(paperFold,'Figures','Fig1','Sub1_percMovCorr');
sessionsv5 = load(fullfile(pathv5Behav,'percMovCorrect_v5.mat')); sessionsv5 = sessionsv5.sessionsv5;

%% Trim
glmGood = false(numel(glmSummary),1);
v5Good = false(numel(sessionsv5),1);
v5lobv = v5Good;
glmlobv = cat(1,glmSummary.iLobV);
% Find matching good sessions
for iGlm = 1:numel(glmGood)    
    for jV5 = 1:numel(v5Good)
        if strcmpi(glmSummary(iGlm).name, sessionsv5(jV5).name) && strcmpi(glmSummary(iGlm).fov, sessionsv5(jV5).fov)
            glmGood(iGlm) = true;
            v5Good(jV5) = true;
            if strcmpi(sessionsv5(jV5).fov,'lobv')
                v5lobv(jV5) = true;
            end
        end
    end
end
% Extract them
seshv5lobv = sessionsv5(v5Good & v5lobv);
seshv5sim2 = sessionsv5(v5Good & ~v5lobv);
glmv5lobv0 = glmSummary(glmGood & glmlobv);
glmv5sim20 = glmSummary(glmGood & ~glmlobv);

% Find v4-v5 match sessions:
iMatchlobv = false(numel(matchData),1);
iMatchsim2 = false(numel(matchData),1);
for iMatch = 1:numel(matchData)
    if matchData(iMatch).v1 == 4 && matchData(iMatch).v2 == 5
        if strcmpi(matchData(iMatch).fov,'lobv')
           iMatchlobv(iMatch) = true; 
        elseif strcmpi(matchData(iMatch).fov,'sim2')
            iMatchsim2(iMatch) = true;
        end
    end
end
% Extract them
matchV5lobv = matchData(iMatchlobv);
matchV5sim2 = matchData(iMatchsim2);
v45Match = cat(2,matchV5lobv,matchV5sim2);
glm45Match = cat(2,glmv5lobv0,glmv5sim20);

clear glmSummary sessionsv5 pathv5Behav glmFold0 glmGood glmlobv iGlm iMatch iMatchlobv iMatchsim2 jV5 v5Good v5lobv

% Load v4 data and extract the same order as the matched FOVs
v4Behav0 = load(fullfile(paperFold,'Figures\Fig3\Sub5_behavCorrelates\spk_gaussian_alpha50','v4MovRewData.mat'));
v4Behav0 = v4Behav0.seshV4;
v4Param = load(fullfile(paperFold, 'Figures\Fig3\Sub5_behavCorrelates\spk_gaussian_alpha50','v4_parametrize.mat'));
v4Param = v4Param.seshSummary;
for iMatch = 1:numel(v45Match)
    for jMatch = 1:numel(v4Behav0)
        if strcmpi(v45Match(iMatch).name,v4Behav0(jMatch).name) && ...
                strcmpi(v45Match(iMatch).fov,v4Behav0(jMatch).fov) && ...
                strcmpi(v45Match(iMatch).date1,v4Behav0(jMatch).v4_Ndate)
        v4Behav(iMatch) = v4Behav0(jMatch);
        end
    end
end

% Load v5 data and extract the same order as the matched FOVs
v5Behav0 = load(fullfile(paperFold,'Figures','Fig5','Sub2_behavSummary','v5_adaptSummary.mat'));
v5Behav = cat(2,v5Behav0.lobvdata,v5Behav0.sim2data);

clear iMatch jMatch matchV5lobv matchV5sim2 seshv5lobv seshv5sim2 v4Behav0 v5Behav0
%% Calculate mov and reward ratios
v45ratios = struct;
for iMatch = 1:numel(v45Match)
    v45ratios(iMatch).name = v45Match(iMatch).name;
    v45ratios(iMatch).fov = v45Match(iMatch).fov;
    v45ratios(iMatch).ilobv = v4Behav(iMatch).ilobv;
    v45ratios(iMatch).isim2 = v4Behav(iMatch).isim2;
    
    % v4 indices (v5 already trimmed)
    iMov4 = find(ismember(v45Match(iMatch).iROIs1,find(glm45Match(iMatch).sigBoolMov)));
    iRew4 = find(ismember(v45Match(iMatch).iROIs1,find(glm45Match(iMatch).sigBoolRew)));
    
    % Reward history - v4:
    spkRew = v4Behav(iMatch).spkRew(:,iRew4,:);
    spkRewBin = logical(spkRew);
    iCorrect = repmat(v4Behav(iMatch).iUCO(:,2),1,3);
    histRew = [false false false; v4Behav(iMatch).iUCO(1:end-1,:)];
    histUCO = iCorrect & histRew;
    spkRewHist = cat(3,mean(spkRew(:,:,histUCO(:,1)),3),mean(spkRew(:,:,histUCO(:,2)),3),mean(spkRew(:,:,histUCO(:,3)),3));
    spkRewHistMax = squeeze(max(spkRewHist(end/2+1:end/2+15,:,:)));
    spkRewBinHist = cat(3,mean(spkRewBin(:,:,histUCO(:,1)),3),mean(spkRewBin(:,:,histUCO(:,2)),3),mean(spkRewBin(:,:,histUCO(:,3)),3));
    spkRewBinHistMax = squeeze(max(spkRewBinHist(end/2+1:end/2+15,:,:)));
    spkRewU = spkRew(:,:,histUCO(:,1));
    spkRewC = spkRew(:,:,histUCO(:,2));
    spkRewO = spkRew(:,:,histUCO(:,3));
    spkRewPopHist = cat(3,mean(spkRew(:,:,histUCO(:,1)),2),mean(spkRew(:,:,histUCO(:,2)),2),mean(spkRew(:,:,histUCO(:,3)),2));
    
    v45ratios(iMatch).v4spkRewOC = spkRewHistMax(:,3)./spkRewHistMax(:,2);
    v45ratios(iMatch).v4spkRewOCbin = spkRewBinHistMax(:,3)./spkRewBinHistMax(:,2);
    
    % Reward base vs mid - v5:
    spkRewHistMax = [max(v5Behav(iMatch).spkRewBase(end/2+1:end/2+15,:))' max(v5Behav(iMatch).spkRewMid(end/2+1:end/2+15,:))'];
    v45ratios(iMatch).v5spkRewMidBase = spkRewHistMax(:,2)./spkRewHistMax(:,1);
    
    iBad0 = [isnan(v45ratios(iMatch).v4spkRewOC) v45ratios(iMatch).v4spkRewOC > 10 isnan(v45ratios(iMatch).v5spkRewMidBase) v45ratios(iMatch).v5spkRewMidBase > 10];
    
    iBad = any(iBad0,2);
    v45ratios(iMatch).v4spkRewOC = v45ratios(iMatch).v4spkRewOC(~iBad);
    v45ratios(iMatch).v4spkRewOCbin = v45ratios(iMatch).v4spkRewOCbin(~iBad);
    v45ratios(iMatch).v5spkRewMidBase = v45ratios(iMatch).v5spkRewMidBase(~iBad);

end

v45lobv = v45ratios(1:5);
v45sim2 = v45ratios(6:end);
v4Rewlobv = cat(1,v45lobv.v4spkRewOC);
v4Rewsim2 = cat(1,v45sim2.v4spkRewOC);
v4RewlobvBin = cat(1,v45lobv.v4spkRewOCbin);
v4Rewsim2Bin = cat(1,v45sim2.v4spkRewOCbin);
v5Rewlobv = single(cat(1,v45lobv.v5spkRewMidBase));
v5Rewsim2 = single(cat(1,v45sim2.v5spkRewMidBase));

%% Plotting
% Expert - reward history dependence
fExpertBin = figure;
subplot(1,2,1); hold on;
axis square; axis([0 10 0 10]);refline(1,0);
for ilobv = 1:numel(v45lobv)
    scatter(v45lobv(ilobv).v4spkRewOC,v45lobv(ilobv).v4spkRewOCbin,60,'filled');
end
olobv = regress_hd(v4Rewlobv,v4RewlobvBin,1);
boundedline(olobv.x_fit,olobv.y_fit,olobv.delta,'cmap',[0 0 0],'alpha');
text(5,9,sprintf('Rsq = %.2f, p = %.2f, m = %.2f',olobv.rsq,olobv.corr_p,olobv.p(1)));
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlabel('Expert reward ratio (Overshoot/Correct)');
ylabel('Binarized expert reward ratio (Overshoot/Correct)');
title('Lobule V: single PC dendrites');

subplot(1,2,2); hold on;
axis square; axis([0 10 0 10]);refline(1,0);
for isim2 = 1:numel(v45sim2)    
    scatter(v45sim2(isim2).v4spkRewOC,v45sim2(isim2).v4spkRewOCbin,60,'filled');
end
osim2 = regress_hd(v4Rewsim2,v4Rewsim2Bin,1);
boundedline(osim2.x_fit,osim2.y_fit,osim2.delta,'cmap',[0 0 0],'alpha');
text(5,9,sprintf('Rsq = %.2f, p = %.2f, m = %.2f',osim2.rsq,osim2.corr_p,osim2.p(1)));
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlabel('Expert reward ratio (Overshoot/Correct)');
ylabel('Binarized expert reward ratio (Overshoot/Correct)');
title('Simplex: single PC dendrites');

suptitle('Reward history: magnitude vs reliability')


% Expert history dependence vs adaptation scaling
fRewCell = figure;
subplot(1,2,1); hold on;
axis square; axis([0 10 0 10]);refline(1,0);
for ilobv = 1:numel(v45lobv)
    scatter(v45lobv(ilobv).v4spkRewOC,v45lobv(ilobv).v5spkRewMidBase,60,'filled');
end
olobv = regress_hd(v4Rewlobv,v5Rewlobv,1);
boundedline(olobv.x_fit,olobv.y_fit,olobv.delta,'cmap',[0 0 0],'alpha');
text(5,9,sprintf('Rsq = %.2f, p = %.2f, m = %.2f',olobv.rsq,olobv.corr_p,olobv.p(1)));
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlabel('Expert reward ratio (Overshoot/Correct)');
ylabel('Expert reward ratio (Mid-Adaptation/Baseline)');
title('Lobule V: single PC dendrites');

subplot(1,2,2); hold on;
axis square; axis([0 10 0 10]);refline(1,0);
for isim2 = 1:numel(v45sim2)    
    scatter(v45sim2(isim2).v4spkRewOC,v45sim2(isim2).v5spkRewMidBase,60,'filled');
end
osim2 = regress_hd(v4Rewsim2,v5Rewsim2,1);
boundedline(osim2.x_fit,osim2.y_fit,osim2.delta,'cmap',[0 0 0],'alpha');
text(5,9,sprintf('Rsq = %.2f, p = %.2f, m = %.2f',osim2.rsq,osim2.corr_p,osim2.p(1)));
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
xlabel('Expert reward ratio (Overshoot/Correct)');
ylabel('Expert reward ratio (Mid-Adaptation/Baseline)');
title('Simplex: single PC dendrites');

suptitle('Reward correlations: expert and adaptation')


end
