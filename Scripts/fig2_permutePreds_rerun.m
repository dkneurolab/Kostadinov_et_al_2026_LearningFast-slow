function GLMoutput = fig2_permutePreds_rerun(GLMoutput,DM,DMdata,Ydata,iFold,iEnds,nPCvar,permPreds,nShuf,binBool)
% (DM,DMdata,DMpcadata,Ytrials,iFoldsmall,iEndTrial,nPCvar,permPreds,nShuf);

rng(2349);
if binBool
    Ydata = double(logical(Ydata));
    distr = 'binomial';
    link = 'logit';
else
    Ydata = double(Ydata);
    distr = 'poisson';
    link = 'log';
end

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
[DMpca.coef,DMpca.score,DMpca.latent,~,DMpca.exp] = pca(DM);

DMshifts = randi(size(DM,1),nShuf,1);


DMpca_shuf = cell(size(Ydata,2),nShuf);
DMpca_shuf_project = DMpca_shuf;

for jShuf = 1:nShuf
    DMshuf = DM;
    % Circular permuation of specific predictors
    DMshuf(:,iPerm) = circshift(DMshuf(:,iPerm),DMshifts(jShuf));
    % Recomputing basis PCA
    [DMshuf_pcaObj.coef,DMshuf_pcaObj.score,DMshuf_pcaObj.latent,~,DMshuf_pcaObj.exp] = pca(DMshuf);    
    % Projecting shuffled DM into original PCA subspace
    DMproject = (DMshuf-mean(DMshuf))*DMpca.coef;
    
    % Taking top n components based on best full model
    for iPC = 1:numel(GLMoutput)
        if ~isempty(GLMoutput(iPC).GLMrankPCA)
            nPCs = GLMoutput(iPC).GLMrankPCA;
            DMpca_shuf{iPC,jShuf} = DMshuf_pcaObj.score(:,1:nPCs);
            DMpca_shuf_project{iPC,jShuf} = DMproject(:,1:nPCs);
        end
    
    end
end

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
% GLMoutput = struct;
optsGLM.alpha = 0.9;
optsGLM.nlambda = 100;

% Partition Ydata
Ydata_tt = Ydata(iTestTrain,:);
Ydata_holdout = Ydata(iHoldout,:);

% % Parition DM
% DMtt = DMpca(iTestTrain,:);
% DMholdout = DMpca(iHoldout,:);

tic;
iCell = 0;
% for iPC = 200
for iPC = 1:size(Ydata,2)
    iCell = iCell+1;
%     GLMoutput(iCell).iPC = iPC;
    
%     if Ymfr(iPC) > 0.1
    if ~isempty(GLMoutput(iCell).devExp) && GLMoutput(iCell).devExp > 0 && GLMoutput(iCell).sigBool > 0 %&& iPC ~= 149
    
    % Shuffle and calculate permuted GLMs and dev explained
    GLMobj_shuf = cell(nShuf,1);
    GLMobj_shuf_proj = cell(nShuf,1);
    devExp_shuf = zeros(nShuf,1);
    devExp_shuf_project = zeros(nShuf,1);
    yHat_shuf = zeros(size(Ydata_holdout,1),nShuf);
    yHat_shuf_proj = zeros(size(Ydata_holdout,1),nShuf);
    for jShuf = 1:nShuf
        % Full shuffle
        DMtt_local = double(DMpca_shuf{iCell,jShuf}(iTestTrain,:));
        GLMshuf = cvglmnet(DMtt_local,Ydata_tt(:,iPC),distr,optsGLM,[],[],circshift(iFold_tt,DMshifts(jShuf)),1);
        DMholdout_local = double(DMpca_shuf{iCell,jShuf}(iHoldout,:));
        CVpred_shuf = cvglmnetPredict(GLMshuf, DMholdout_local, 'lambda_min');
        yHat_shuf(:,jShuf) = exp(CVpred_shuf);
        devExp_shuf(jShuf) = poissonDevExp(Ydata_holdout(:,iPC),yHat_shuf(:,jShuf));
        GLMobj_shuf{jShuf,1} = GLMshuf;
        
        % Projection shuffle
        DMtt_local_proj = double(DMpca_shuf_project{iCell,jShuf}(iTestTrain,:));
        GLMshuf_proj = cvglmnet(DMtt_local_proj,Ydata_tt(:,iPC),distr,optsGLM,[],[],circshift(iFold_tt,DMshifts(jShuf)),1);
        DMholdout_local = double(DMpca_shuf_project{iCell,jShuf}(iHoldout,:));
        CVpred_shuf_proj = cvglmnetPredict(GLMshuf_proj, DMholdout_local, 'lambda_min');
        yHat_shuf_proj(:,jShuf) = exp(CVpred_shuf_proj);
        devExp_shuf_project(jShuf) = poissonDevExp(Ydata_holdout(:,iPC),yHat_shuf_proj(:,jShuf));
        GLMobj_shuf_proj{jShuf,1} = GLMshuf_proj;
        
    end
    % Full shuffle
    GLMoutput(iCell).GLMobj_shuf = GLMobj_shuf;
    GLMoutput(iCell).yHat_shuf_mu = single(mean(yHat_shuf,2));
    GLMoutput(iCell).yHat_shuf_sd = single(std(yHat_shuf,[],2));
    GLMoutput(iCell).devExp_shuf = devExp_shuf;
    % Projection shuffle
    GLMoutput(iCell).GLMobj_shuf_proj = GLMobj_shuf_proj;
    GLMoutput(iCell).yHat_shuf_proj_mu = single(mean(yHat_shuf_proj,2));
    GLMoutput(iCell).yHat_shuf_proj_sd = single(std(yHat_shuf_proj,[],2));
    GLMoutput(iCell).devExp_shuf_proj = devExp_shuf_project;
    
    GLMoutput(iCell).sigBool_shuf = GLMoutput(iCell).devExp > (mean(devExp_shuf)+2*std(devExp_shuf)) & GLMoutput(iCell).devExp > 0;
    GLMoutput(iCell).sigBool_proj = GLMoutput(iCell).devExp > (mean(devExp_shuf_project)+2*std(devExp_shuf_project)) & GLMoutput(iCell).devExp > 0;
    end    
    if mod(iCell,10) == 0 || iPC == size(Ydata,2)
        tNow = clock;
        fprintf('Done %i cells in %i seconds, time is %i:%02.f on %i-%02.f-%02.f.\n',iCell,round(toc),tNow(4),tNow(5),tNow(1),tNow(2),tNow(3));    
    end

end

end