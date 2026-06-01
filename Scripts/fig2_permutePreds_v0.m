function GLMoutput = fig2_permutePreds_v0(DM,DMdata,Ydata,iFold,iEnds,nPCvar,permPreds,nShuf)
% fig2_permutePreds(DM,DMdata,Ydata,iFoldsmall,iEndTrial,nPCvar,permPreds,nShuf);

rng(2349);
distr = 'poisson';

%% Find indices of predictors to permute

fields = DMdata.fields;
nBaseEnds = [0; cumsum(DMdata.nBases)];

permInds = zeros(size(permPreds));
iPerm = [];
a = 0;
for iPred = 1:numel(permPreds)
    for jField = 1:numel(fields)
        if strcmpi(permPreds(iPred),fields(jField))
            a = a + 1;
            permInds(a) = jField;
            iPerm = [iPerm; (nBaseEnds(jField)+1:nBaseEnds(jField+1))']; %#ok<AGROW>
        end
    end
end

%% Permute and run PCA
% [DMpca.coef,DMpca.score,DMpca.latent,~,DMpca.exp] = pca(DM);

DMshifts = randi(size(DM,1),nShuf,1);

DMshuf_pca = zeros(size(DM,1),nPCvar,nShuf);
for jShuf = 1:nShuf
    DMshuf = DM;
    DMshuf(:,iPerm) = circshift(DMshuf(:,iPerm),DMshifts(jShuf));
    [DMpca.coef,DMpca.score,DMpca.latent,~,DMpca.exp] = pca(DMshuf);
    
    DMshuf_pca(:,:,jShuf) = DMpca.score(:,1:nPCvar);
end

% Get real DMpca
[DMpca.coef,DMpca.score,DMpca.latent,~,DMpca.exp] = pca(DM);
DMpca = DMpca.score(:,1:nPCvar);

%% Make CV partitions, etc.
% Make array of trial indices
iTrial = zeros(size(DM,1),1);
for i = 1:numel(iEnds)-1
    iTrial(iEnds(i)+1:iEnds(i+1)) = i;
end
% Split into test-train+holdout
iHoldout = iFold == max(iFold);
iTestTrain = iFold ~= max(iFold);

% Make cv partition object
iFold_tt = iFold(iTestTrain);

Ymfr = mean(Ydata,1)*30;

%% Run GLM
GLMoutput = struct;
optsGLM.alpha = 0.9;
optsGLM.nlambda = 100;

% Partition Ydata
Ydata_tt = Ydata(iTestTrain,:);
Ydata_holdout = Ydata(iHoldout,:);

% Parition DM
DMtt = DMpca(iTestTrain,:);
DMholdout = DMpca(iHoldout,:);

tic;
iCell = 0;
% for iPC = 200
for iPC = 1:size(Ydata,2)
    iCell = iCell+1;
    GLMoutput(iCell).iPC = iPC;
    
    if Ymfr(iPC) > 0.1
    % Run GLM
    GLMoutput(iCell).GLMobj = cvglmnet(DMtt,Ydata_tt(:,iPC),distr,optsGLM,[],[],iFold_tt,1);
    
    % Predict on heldout data
    CVpred = cvglmnetPredict(GLMoutput(iCell).GLMobj, DMholdout, 'lambda_min');
    yHat = exp(CVpred);
    
    GLMoutput(iCell).yData = single(Ydata_holdout(:,iPC));
    GLMoutput(iCell).yHat = single(yHat);
    
    % Calculate real deviance explained
    GLMoutput(iCell).devExp = poissonDevExp(Ydata_holdout(:,iPC),yHat);
    
    % Shuffle and calculate permuted GLMs and dev explained
    devExp_shuf = zeros(nShuf,1);
    yHat_shuf = zeros(size(Ydata_holdout,1),nShuf);
    for jShuf = 1:nShuf
        DMtt_local = DMshuf_pca(iTestTrain,:,jShuf);
        GLMshuf = cvglmnet(DMtt_local,Ydata_tt(:,iPC),distr,optsGLM,[],[],circshift(iFold_tt,DMshifts(jShuf)),1);
        DMholdout_local = DMshuf_pca(iHoldout,:,jShuf);
        CVpred_shuf = cvglmnetPredict(GLMshuf, DMholdout_local, 'lambda_min');
        yHat_shuf(:,jShuf) = exp(CVpred_shuf);
        devExp_shuf(jShuf) = poissonDevExp(Ydata_holdout(:,iPC),yHat_shuf(:,jShuf));
    end
    GLMoutput(iCell).yHat_shuf_mu = single(mean(yHat_shuf,2));
    GLMoutput(iCell).yHat_shuf_sd = single(std(yHat_shuf,[],2));
    GLMoutput(iCell).devExp_shuf = devExp_shuf;
    GLMoutput(iCell).sigBool = GLMoutput(iCell).devExp > (mean(devExp_shuf)+2*std(devExp_shuf)) & GLMoutput(iCell).devExp > 0;
    
    if mod(iCell,10) == 0 || iPC == size(Ydata,2)
        fprintf('Done %i cells in %i seconds.\n',iCell,round(toc));
    end
   
    end
end

end