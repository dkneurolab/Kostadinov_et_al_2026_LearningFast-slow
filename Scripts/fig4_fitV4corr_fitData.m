function GLMoutput = fig4_fitV4corr_fitData(GLMinput,DM,DMdata,Ydata,glmParams,permGrp)
%% Start output structure that will have the relevant information
nShuf = numel(GLMinput(1).GLMobj_shuf);
glmFields = fieldnames(GLMinput);
GLMoutput = rmfield(GLMinput,glmFields([2,4:11,13:end]));
xvalFold = 5;

rng(2349);

%% Find indices of predictors to that were permuted and will now be tested
switch permGrp
    case 'Vis'
%         permPreds = {'TrialON','BSwheelpos'};
        permPreds = {'TrialON'};
    case 'Mov'
        permPreds = {'MovT'};
%         permPreds = {'MovT','BSwheelvelL','BSwheelvelR'};
    case 'Rew'
        permPreds = {'RewTcorrect'};
    case 'Lick'
        permPreds = {'LickT'};
end

yHat = sprintf('yHat_shuf%s_proj_mu',permGrp);
sigBoolLocal = sprintf('sigBool_shuf%s_proj',permGrp);
GLMobj = sprintf('GLMobj_shuf%s_proj',permGrp);


%% Permute and run PCA
% Figure out which indices to shift
dmFields = DMdata.fields;
nBaseEnds = [0; cumsum(DMdata.nBases)];
permInds = zeros(size(permPreds));
iPerm = [];
a = 0;
for iPred = 1:numel(permPreds)
    for jField = 1:numel(dmFields)
        if strcmpi(permPreds(iPred),dmFields(jField))
            a = a + 1;
            permInds(a) = jField;
            iPerm = [iPerm; (nBaseEnds(jField)+1:nBaseEnds(jField+1))']; %#ok<AGROW>
        end
    end
end

% Run PCA
[DMpca.coef,DMpca.score,DMpca.latent,~,DMpca.exp] = pca(DM);
% Pre-allocate shifts
DMshifts = randi(size(DM,1),nShuf,1);
DMpca_shuf_project = cell(nShuf,1);
for jShuf = 1:nShuf
    DMshuf = DM;
    % Circular permuation of specific predictors
    DMshuf(:,iPerm) = circshift(DMshuf(:,iPerm),DMshifts(jShuf));
    % Recomputing basis PCA
    [DMshuf_pcaObj.coef,DMshuf_pcaObj.score,DMshuf_pcaObj.latent,~,DMshuf_pcaObj.exp] = pca(DMshuf);    
    % Projecting shuffled DM into original PCA subspace
    DMpca_shuf_project{jShuf,1} = (DMshuf-mean(DMshuf))*DMpca.coef;
    
%     % Taking top n components based on best full model
%     for iPC = 1:numel(GLMoutput)
%         if ~isempty(GLMoutput(iPC).GLMrankPCA)
%             nPCs = GLMoutput(iPC).GLMrankPCA;
%             DMpca_shuf_project{iPC,jShuf} = DMproject(:,1:nPCs);
%         end
%     end
end

% Make DM for residual predictors
DMres = double(DM(:,iPerm));
iDMresMid = floor(size(DMres,2));


%% Make CV partitions, etc.
iEnds = DMdata.endTrialIndices;
% Make array of trial indices
iTrial = zeros(size(DM,1),1);
for i = 1:numel(iEnds)-1
    iTrial(iEnds(i)+1:iEnds(i+1)) = i;
end
iFold = mod(iTrial,xvalFold);
iFold(iFold == 0) = xvalFold;

%% Iterate through PCs and run analysis
% Set GLM parameters
optsGLM.alpha = glmParams.alpha; % Elastic net hyperparam - somewhere between 0.5 and 0.9
optsGLM.nlambda = glmParams.nLambda; % Number of penalties to try - probably 100;
distr = glmParams.distr;

tic;
for iPC = 1:size(Ydata,2)
    if ~isempty(GLMinput(iPC).(yHat)) && sum(DMres(:,iDMresMid)) >= xvalFold
    % Shuffle and calculate permuted GLMs and dev explained
    GLMoutput(iPC).yData = single(Ydata(:,iPC));    
    GLMobjToFit = GLMinput(iPC).(GLMobj);        
    GLMoutput(iPC).sigBoolLocal = GLMinput(iPC).(sigBoolLocal);
    yHat_shuf_proj = zeros(size(Ydata,1),nShuf);
    for jShuf = 1:numel(nShuf)
        CVpred_shuf_proj = cvglmnetPredict(GLMobjToFit{jShuf}, DMpca_shuf_project{jShuf}(:,1:GLMoutput(iPC).GLMrankPCA), 'lambda_min');        
        if strcmpi(glmParams.distr,'gaussian')
            yHat_shuf_proj(:,jShuf) = CVpred_shuf_proj;                                            
        elseif strcmpi(glmParams.distr,'poisson')
            yHat_shuf_proj(:,jShuf) = exp(CVpred_shuf);            
        end
    end
    
    % Fit GLM on residual with appropriate predictors
    yRes = double(GLMoutput(iPC).yData-mean(yHat_shuf_proj,2));
    GLMoutput(iPC).GLMobjLocal = cvglmnet(DMres,yRes,distr,optsGLM,[],[],iFold,1); % 1 = Parallel, 0 = Serial
    iBest = find(GLMoutput(iPC).GLMobjLocal.lambda == GLMoutput(iPC).GLMobjLocal.lambda_min);
    GLMoutput(iPC).kernelLocal = GLMoutput(iPC).GLMobjLocal.glmnet_fit.beta(:,iBest);
    
    else
        GLMoutput(iPC).sigBoolLocal = false;
    end
    
    if mod(iPC,50) == 0 || iPC == numel(GLMinput)
        tNow = clock;
        fprintf('Done %i cells in %i seconds, time is %i:%02.f on %i-%02.f-%02.f.\n',iPC,round(toc),tNow(4),tNow(5),tNow(1),tNow(2),tNow(3));        
    end
end

end


%% Make CV partitions, etc.
% iEnds = DMdata.endTrialIndices;
% % Make array of trial indices
% iTrial = zeros(size(DM,1),1);
% for i = 1:numel(iEnds)-1
%     iTrial(iEnds(i)+1:iEnds(i+1)) = i;
% end
% iFold = mod(iTrial,xvalFold);
% iFold(iFold == 0) = xvalFold;
% 
% 
% 
% yHat = sprintf('yHat_shuf%s_proj_mu',permGrp);
% sigBoolLocal = sprintf('sigBool_shuf%s_proj',permGrp);
% 
% 
% 
% permInds = zeros(size(permPreds));
% iPerm = [];
% a = 0;
% for iPred = 1:numel(permPreds)
%     for jField = 1:numel(dmFields)
%         if strcmpi(permPreds(iPred),dmFields(jField))
%             a = a + 1;
%             permInds(a) = jField;
%             iPerm = [iPerm; (nBaseEnds(jField)+1:nBaseEnds(jField+1))']; %#ok<AGROW>
%         end
%     end
% end
% 
% % Make DM with only heldout data
% iHeldout = iFold == max(iFold);
% DMres = double(DM(iHeldout,iPerm));
% iDMres = logical(sum(DMres,2));
% 
% 
% % Find ends of each holdout trial and use them to make 5-folds
% iHeldoutEnds = diff(iFold) < 0; iHeldoutEnds(end+1) = true;
% iHeldoutEnds = [0; find(iHeldoutEnds(iHeldout))];
% iFoldInds = mod([1:numel(iHeldoutEnds)-1]',xvalFold); iFoldInds(iFoldInds == 0) = xvalFold;
% iFoldHoldout = zeros(sum(iHeldout),1);
% for i = 1:numel(iHeldoutEnds)-1
%     iFoldHoldout(iHeldoutEnds(i)+1:iHeldoutEnds(i+1)) = iFoldInds(i);
% end
% 
% clear glmFields dmFields nBaseEnds a iPred jField
% 
% %% Iterate through PCs and run analysis
% % Set GLM parameters
% optsGLM.alpha = glmParams.alpha; % Elastic net hyperparam - somewhere between 0.5 and 0.9
% optsGLM.nlambda = glmParams.nLambda; % Number of penalties to try - probably 100;
% distr = glmParams.distr;
% 
% tic;
% for iPC = 1:numel(GLMinput)
%     if ~isempty(GLMinput(iPC).(yHat)) && sum(DMres(:,1)) >= xvalFold
% %         yRes = double(GLMinput(iPC).yData(iDMres)-GLMinput(iPC).(yHat)(iDMres));
% %         GLMoutput(iPC).sigBoolLocal = GLMinput(iPC).(sigBoolLocal);
% %         
% %         % Run GLM - last boolean signifies parallel or non
% %         GLMoutput(iPC).GLMobjLocal = cvglmnet(DMres(iDMres,:),yRes,distr,optsGLM,[],[],iFoldHoldout(iDMres),1); % 1 = Parallel, 0 = Serial
% %         iBest = find(GLMoutput(iPC).GLMobjLocal.lambda == GLMoutput(iPC).GLMobjLocal.lambda_min);
% %         GLMoutput(iPC).kernelLocal = GLMoutput(iPC).GLMobjLocal.glmnet_fit.beta(:,iBest);
%                 
%         yRes = double(GLMinput(iPC).yData-GLMinput(iPC).(yHat));
%         GLMoutput(iPC).sigBoolLocal = GLMinput(iPC).(sigBoolLocal);
%         % Run GLM - last boolean signifies parallel or non
%         GLMoutput(iPC).GLMobjLocal = cvglmnet(DMres,yRes,distr,optsGLM,[],[],iFoldHoldout,1); % 1 = Parallel, 0 = Serial
%         iBest = find(GLMoutput(iPC).GLMobjLocal.lambda == GLMoutput(iPC).GLMobjLocal.lambda_min);
%         GLMoutput(iPC).kernelLocal = GLMoutput(iPC).GLMobjLocal.glmnet_fit.beta(:,iBest);
%         
%         aa =4;
%     else
%         GLMoutput(iPC).sigBoolLocal = false;
%     end
%     
%     if mod(iPC,50) == 0 || iPC == numel(GLMinput)
%         tNow = clock;
%         fprintf('Done %i cells in %i seconds, time is %i:%02.f on %i-%02.f-%02.f.\n',iPC,round(toc),tNow(4),tNow(5),tNow(1),tNow(2),tNow(3));        
%     end
% end
% 
% end