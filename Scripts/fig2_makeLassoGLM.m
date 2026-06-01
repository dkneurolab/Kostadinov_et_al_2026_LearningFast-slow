function GLMoutput = fig2_makeLassoGLM(DMpca,DM, Ydata,glmParams,iFold,iEnds,nShuf)
%% Function information:
% DM = design matrix (should be PCA'd), shape = M (data points) x N predictors
% Ydata is full data matrix, shape = M (data points) x P neurons
% iFold is folds of cross-validation for grouping, shape = M (data point) x 1,
%    with values going from 1:CV folds + 1, CV+1 fold is holdout data
rng(2349);
%% Make CV partitions, etc.
% Make array of trial indices
iTrial = zeros(size(DMpca,1),1);
for i = 1:numel(iEnds)-1
    iTrial(iEnds(i)+1:iEnds(i+1)) = i;
end

% Split into test-train+holdout
iHoldout = iFold == max(iFold);
iTestTrain = iFold ~= max(iFold);

% Parition DM
DMpca = double(DMpca);
DM = double(DM);
DMtt_full = DMpca(iTestTrain,:);
DMholdout_full = DMpca(iHoldout,:);
DMshifts = randi(size(DMtt_full,1),nShuf,1);

Ydata_tt = Ydata(iTestTrain,:);
Ydata_holdout = Ydata(iHoldout,:);

% Make cv partition object
iFold_tt = iFold(iTestTrain);

Ymfr = mean(Ydata)*30;

%% Run GLM
GLMoutput = struct;
optsGLM.alpha = glmParams.alpha; % Elastic net hyperparam - somewhere between 0.5 and 0.9
optsGLM.nlambda = glmParams.nLambda; % Number of penalties to try - probably 100;
predMax = glmParams.maxPred; % Max number of PCA predictors to use
distr = glmParams.distr;

iEx = 0;
tic;
iCell = 0;
for iPC = 1:size(Ydata,2)
    iCell = iCell+1;
    GLMoutput(iCell).iPC = iPC;
    
    if Ymfr(iPC) > 0.1
    % Run GLM with different number of PCA'd predictors to figure out best DM size
    nPCpreds = nan*zeros(size(DMpca,2),1);
    GLMtemp = struct;
    % Initial parameter sweep
    for iPreds = 5:10:predMax
        DMtt = DMtt_full(:,1:iPreds);
        DMholdout = DMholdout_full(:,1:iPreds);
        
        % Run GLM - last boolean signifies parallel or non
        GLMtemp(iPreds).GLMobj = cvglmnet(DMtt,Ydata_tt(:,iPC),distr,optsGLM,[],[],iFold_tt,1); % 1 = Parallel, 0 = Serial
        % Predict on heldout data
        CVpred = cvglmnetPredict(GLMtemp(iPreds).GLMobj, DMholdout, 'lambda_min');
        % Calculate real deviance explained
        if strcmpi(glmParams.distr,'gaussian')
            yHat = CVpred;
            nPCpreds(iPreds) = gaussianDevExp(Ydata_holdout(:,iPC),yHat);
        elseif strcmpi(glmParams.distr,'poisson')
            yHat = exp(CVpred);
            nPCpreds(iPreds) = poissonDevExp(Ydata_holdout(:,iPC),yHat);
        end
        
        % Save just in case useful:
        GLMtemp(iPreds).yData = single(Ydata_holdout(:,iPC));
        GLMtemp(iPreds).yHat = single(yHat);
    end
    
    % Run again over 2 stretches near peaks to find true maximum:
    [~,iMaxPred] = maxk(nPCpreds,2);
    predVals = unique([iMaxPred(1)-9:iMaxPred(1)+9; iMaxPred(2)-9:iMaxPred(2)+9]');
    predVals(predVals > predMax) = []; predVals(predVals < 1) = [];
    for iPreds = 1:numel(predVals)
        if isnan(nPCpreds(predVals(iPreds)))
            DMtt = DMtt_full(:,1:predVals(iPreds));
            DMholdout = DMholdout_full(:,1:predVals(iPreds));
            
            % Run GLM - last boolean signifies parallel or non
            GLMtemp(predVals(iPreds)).GLMobj = cvglmnet(DMtt,Ydata_tt(:,iPC),distr,optsGLM,[],[],iFold_tt,1); % 1 = Parallel, 0 = Serial
            % Predict on heldout data
            CVpred = cvglmnetPredict(GLMtemp(predVals(iPreds)).GLMobj, DMholdout, 'lambda_min');
            % Calculate real deviance explained
            if strcmpi(glmParams.distr,'gaussian')
                yHat = CVpred;
                nPCpreds(predVals(iPreds)) = gaussianDevExp(Ydata_holdout(:,iPC),yHat);
            elseif strcmpi(glmParams.distr,'poisson')
                yHat = exp(CVpred);
                nPCpreds(predVals(iPreds)) = poissonDevExp(Ydata_holdout(:,iPC),yHat);
            end
            
            % Save just in case useful:
            GLMtemp(predVals(iPreds)).yData = single(Ydata_holdout(:,iPC));
            GLMtemp(predVals(iPreds)).yHat = single(yHat);
        end
    end
    
    % Run full GLM - last boolean signifies parallel or non
    GLMobj_full = cvglmnet(DM(iTestTrain,:),Ydata_tt(:,iPC),distr,optsGLM,[],[],iFold_tt,1); % 1 = Parallel, 0 = Serial
    CVpred_full = cvglmnetPredict(GLMobj_full, DM(iHoldout,:), 'lambda_min');
    if strcmpi(glmParams.distr,'gaussian')
        fullPred = gaussianDevExp(Ydata_holdout(:,iPC),CVpred_full);
    elseif strcmpi(glmParams.distr,'poisson')
        fullPred = poissonDevExp(Ydata_holdout(:,iPC),exp(CVpred_full));
    end
    
    if fullPred > 0.1
        iEx = iEx + 1;
        fEx(iEx) = figure;             %#ok<AGROW>
        plot(1:numel(nPCpreds),nPCpreds,'or')
        hold on;
        plot(5:5:numel(nPCpreds),nPCpreds(5:5:numel(nPCpreds)),'ok')
        plot(0,fullPred,'ob');
        set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
        xlim([-5 numel(nPCpreds)]);
        ylim([0 1.2*max([nPCpreds;fullPred])]);
        box off;
        title(sprintf('PC num: %i',iPC));
    end
    
    % Take best GLMobject and data
    [~,iPredBest] = max(nPCpreds);
    GLMoutput(iCell).GLMobj = GLMtemp(iPredBest).GLMobj;
    GLMoutput(iCell).GLMrankPCA = iPredBest;
    GLMoutput(iCell).yData = GLMtemp(iPredBest).yData;
    GLMoutput(iCell).yHat = GLMtemp(iPredBest).yHat;
    GLMoutput(iCell).devExp = nPCpreds(iPredBest);
    GLMoutput(iCell).GLMobj_full = GLMobj_full;
    GLMoutput(iCell).yHat_full = single(exp(CVpred_full));
    GLMoutput(iCell).devExp_full = fullPred;
    
    % Now run shuffles to determine significance of fit:
    DMtt = DMtt_full(:,1:iPredBest);
    devExp_shuf = zeros(nShuf,1);
    GLMobjs = cell(nShuf,1);
    for jShuf = 1:nShuf
        DMshuf = circshift(DMtt,DMshifts(jShuf));
        GLMshuf = cvglmnet(DMshuf,Ydata_tt(:,iPC),distr,optsGLM,[],[],circshift(iFold_tt,DMshifts(jShuf)),1);
        %             GLMshuf = cvglmnet(DMshuf,Ydata_tt(:,iPC),distr,optsGLM,[],[],circshift(iFold_tt,DMshifts(jShuf)),0);
        CVpred_shuf = cvglmnetPredict(GLMshuf, DMholdout_full(:,1:iPredBest), 'lambda_min');
        
        % Calculate real deviance explained
        if strcmpi(glmParams.distr,'gaussian')
            yHat_shuf = CVpred_shuf;
            devExp_shuf(jShuf) = gaussianDevExp(Ydata_holdout(:,iPC),yHat_shuf);
        elseif strcmpi(glmParams.distr,'poisson')
            yHat_shuf = exp(CVpred_shuf);
            devExp_shuf(jShuf) = poissonDevExp(Ydata_holdout(:,iPC),yHat_shuf);
        end
        
        GLMobjs{jShuf,1} = GLMshuf;
    end
    GLMoutput(iCell).GLMobj_shuf = GLMobjs;
    GLMoutput(iCell).devExp_shuf = devExp_shuf;
    GLMoutput(iCell).sigBool = GLMoutput(iCell).devExp > (mean(devExp_shuf)+2*std(devExp_shuf)) & GLMoutput(iCell).devExp > 0;
    
    if mod(iPC,5) == 0 || iPC == size(Ydata,2)
        tNow = clock;
        fprintf('Done %i cells in %i seconds, time is %i:%02.f on %i-%02.f-%02.f.\n',iCell,round(toc),tNow(4),tNow(5),tNow(1),tNow(2),tNow(3));
    end
    end
end

end
