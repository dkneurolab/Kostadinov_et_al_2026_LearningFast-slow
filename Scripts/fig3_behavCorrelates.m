function fig3_behavCorrelates(matchFOVs,glmParams,dataFold0,paperFold,plotbool,savebool)
%% Things that need to be done here:
% - Load mov and rew spkmat and behaviour - DONE
% - Identify mov- and reward-selective cells in each FOV - DONE
% - Parametrize movement by:
%     1. Amplitude
%     2. Max velocity
%     3. Max accelaration
%     4. correctness
%     5. Activity amount
% - Analyze reward history dependence!

%% Housekeeping:
% Make path to find GLM info and save stuff
glmSubDir = sprintf('GLMfull_%s_%s_alpha%i',glmParams.dataType,glmParams.distr,glmParams.alpha*100);
figFold0 = fullfile(paperFold,'Figures','Fig3','Sub5_behavCorrelates');
figFold = fullfile(figFold0,sprintf('%s_%s_alpha%i',glmParams.dataType,glmParams.distr,glmParams.alpha*100));
if ~exist(figFold,'dir'); mkdir(figFold); end

%% Load data if it exists or iterate and extract data we need!

if ~exist(fullfile(figFold,'v4MovRewData.mat'),'file')
    seshV4 = rmfield(matchFOVs,{'v1date','iv1sesh','v4_1date','iv4_1sesh','v5date','iv5sesh','iGrp'});
    nDec = 50;
    
    for iGood = 1:numel(seshV4)
        % Load trial struct
        dataFold = fullfile(dataFold0,'01_Behav+imaging','Version_4',seshV4(iGood).name,sprintf('%s_%s',seshV4(iGood).v4_Ndate,seshV4(iGood).fov));
        trialData = load(fullfile(dataFold,'trialsonly.mat'));
        trialSets = trialData.sets;
        trialData = trialData.trialstructs;
        trialStruct = cat(2,trialData.LUstruct,trialData.LCstruct,trialData.LOstruct);
        trialStruct = rmfield(trialStruct,{'FLon','zFLon','spkLon','VRLon','indmovvr','indmov2p','FLmov','zFLmov','FLoff','zFLoff',...
            'indcornegvr','indcorneg2p','Fneg','zFneg','spkneg','VRneg','indcorposvr','indcorpos2p','Fpos','zFpos','spkpos','VRpos'});
        nUCO = [numel(trialData.LUstruct); numel(trialData.LCstruct); numel(trialData.LOstruct)];
        sumUCO = cumsum(nUCO);
        gainDiv = 37.5/360*4.5*pi/trialSets.vrgain;
        smspan = trialSets.fs*.06; % 60 ms smoothing span
        vrfsDec = trialSets.fs/nDec;
        imfs = trialSets.imfs;
        emptyFlag = false(numel(trialStruct),1);
        
        for jUCO = 1:sumUCO(end)
            if isempty(trialStruct(jUCO).trialind)
                emptyFlag(jUCO) = true;
            else
                if jUCO <= sumUCO(1)
                    trialStruct(jUCO).iUCO = [1 0 0];
                elseif jUCO > sumUCO(1) && jUCO <= sumUCO(2)
                    trialStruct(jUCO).iUCO = [0 1 0];
                else
                    trialStruct(jUCO).iUCO = [0 0 1];
                end
                
                % Process movement behaviour
                x0 = trialStruct(jUCO).VRLmov(:,3);
                x1 = smooth(x0,smspan/numel(x0),'lowess');
                v1 = diff(x1)*trialSets.fs*-gainDiv; v1(end+1) = v1(end);
                x2 = -(x1-mean(x1((trialSets.imtime-1)*trialSets.fs+1:(trialSets.imtime-1)*trialSets.fs+trialSets.fs/2)));
                trialStruct(jUCO).VRLmov = [decimate(x2,50) decimate(v1,50)];
                
                % Process rewT behaviour
                x0 = trialStruct(jUCO).VRLoff(:,3);
                x1 = smooth(x0,smspan/numel(x0),'lowess');
                v1 = diff(x1)*trialSets.fs*-gainDiv; v1(end+1) = v1(end);
                x2 = -(x1-mean(x1((trialSets.imtime-1)*trialSets.fs+1:(trialSets.imtime-1)*trialSets.fs+trialSets.fs/2)));
                trialStruct(jUCO).VRLoff = [decimate(x2,nDec) decimate(v1,nDec)];
            end
        end
        trialStruct = trialStruct(~emptyFlag);
        
        iTrials = cat(1,trialStruct.trialind);
        [~,iSort] = sort(iTrials,'ascend');
        trialStruct = trialStruct(iSort);
        
        % Take the good bits from the trial struct and put them into seshV4
        seshV4(iGood).iTrial = cat(1,trialStruct.trialind);
        seshV4(iGood).spkMov = single(cat(3,trialStruct.spkLmov));
        VRmov = single(cat(3,trialStruct.VRLmov));
        seshV4(iGood).xMov = single(squeeze(VRmov(:,1,:)));
        seshV4(iGood).vMov = single(squeeze(VRmov(:,2,:)));
        seshV4(iGood).spkRew = single(cat(3,trialStruct.spkLoff));
        VRoff = single(cat(3,trialStruct.VRLoff));
        seshV4(iGood).xRew = single(squeeze(VRoff(:,1,:)));
        seshV4(iGood).vRew = single(squeeze(VRoff(:,2,:)));
        seshV4(iGood).spkMov = single(cat(3,trialStruct.spkLmov));
        seshV4(iGood).iUCO = logical(cat(1,trialStruct.iUCO));
        
        
        % Get correlation with mean rewarded movement
        iUCO = cat(1,trialStruct.iUCO);
        vMovC = cat(3,trialStruct(logical(iUCO(:,2))).VRLmov);
        vMovCMu = single(squeeze(mean(vMovC(:,2,:),3)));
        seshV4(iGood).vMovCMu = vMovCMu;
        movCorrs = corrcoef([seshV4(iGood).vMov(end/2-vrfsDec/2+1:end/2+vrfsDec,:) vMovCMu(end/2-vrfsDec/2+1:end/2+vrfsDec)]);
        seshV4(iGood).RmovC = movCorrs(1:end-1,end);
        seshV4(iGood).lastTrialC = [false; seshV4(iGood).iUCO(1:end-1,2)];
        
        % ID movement and reward neurons
        glmFold = fullfile(dataFold,'GLM','glmRfiles',glmSubDir);
        glmStruct = load(fullfile(glmFold,'GLMsig+VisMovRewLick.mat')); glmStruct = glmStruct.GLMoutput;
        seshV4(iGood).iPCsig = false(numel(glmStruct),1);
        seshV4(iGood).iPCmov = seshV4(iGood).iPCsig; seshV4(iGood).iPCrew = seshV4(iGood).iPCsig;
        for jLocal = 1:numel(glmStruct)
            if ~isempty(glmStruct(jLocal).sigBool)
                seshV4(iGood).sigBool(jLocal,1) = glmStruct(jLocal).sigBool;
            end
            if ~isempty(glmStruct(jLocal).sigBool_shufMov_proj)
                seshV4(iGood).iPCmov(jLocal,1) = glmStruct(jLocal).sigBool_shufMov_proj;
            end
            if ~isempty(glmStruct(jLocal).sigBool_shufRew_proj)
                seshV4(iGood).iPCrew(jLocal,1) = glmStruct(jLocal).sigBool_shufRew_proj;
            end
        end
    end
    
    if savebool
        save(fullfile(figFold,'v4MovRewData.mat'),'seshV4','imfs','vrfsDec');
    end
else
    seshV4 = load(fullfile(figFold,'v4MovRewData.mat'));
    imfs = seshV4.imfs;
    vrfsDec = seshV4.vrfsDec;
    seshV4 = seshV4.seshV4;    
end

%% Now analyze movement by parameters listed below (separate loop for ease)

ilobv = find(cat(1,seshV4.ilobv));
isim2 = find(cat(1,seshV4.isim2));

seshV4 = seshV4([ilobv; isim2]);

% plotbool = false;

for iSesh = 1:numel(seshV4)
    % Pass into another fxn that will process movement
    seshOut = fig3_movAnalysis(seshV4(iSesh),imfs,vrfsDec,figFold,plotbool,savebool);    
    % Pass into another fxn that will process reward history
    seshOut = fig3_rewAnalysis(seshV4(iSesh),seshOut,imfs,vrfsDec,figFold,plotbool,savebool);
    seshSummary(iSesh) = seshOut; %#ok<AGROW>
end


%% Make summary plots

ilobv = [1:numel(ilobv)]';
isim2 = [numel(ilobv)+1:numel(seshSummary)]';

[figsField, statsField] = fig3_summary_plot(seshSummary,ilobv,isim2);
% figsCell = fig3_summaryCell_plot(seshSummary,ilobv,isim2);




if savebool
    figname = 'v4_movSummary_field';
    savefig(figsField(1), fullfile(figFold,[figname,'.fig']));
    saveas(figsField(1),fullfile(figFold,[figname,'.png']));
    print(figsField(1),fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    
    figname = 'v4_rewSummary_field';
    savefig(figsField(2), fullfile(figFold,[figname,'.fig']));
    saveas(figsField(2),fullfile(figFold,[figname,'.png']));
    print(figsField(2),fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');    
    save(fullfile(figFold,'v4_parametrize.mat'),'seshSummary');
    
    figname = 'v4_sigSummary_field';
    savefig(figsField(3), fullfile(figFold,[figname,'.fig']));
    saveas(figsField(3),fullfile(figFold,[figname,'.png']));
    print(figsField(3),fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');    
    
    save(fullfile(figFold,'v4_parametrize.mat'),'seshSummary', 'statsField');
end
       
end