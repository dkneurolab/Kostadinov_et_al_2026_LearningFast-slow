function GLMlocal = fig4_fitKernels_crossDay(GLMmatched, glmDataTest, glmParams)
%% Set up and unpack data struct
glmFields = fieldnames(GLMmatched.GLMstruct);
GLMlocal = rmfield(GLMmatched.GLMstruct,glmFields([3:4,6:17,20,23:end]));
xvalFold = 5;
v2struct(glmDataTest);

% Find ends of each holdout trial and use them to make 5-folds
iHeldoutEnds = DMdata.endTrialIndices;
iFoldInds = mod([1:numel(iHeldoutEnds)-1]',xvalFold); iFoldInds(iFoldInds == 0) = xvalFold;
iFold = zeros(size(Ytrials,1),1);
for i = 1:numel(iHeldoutEnds)-1
    iFold(iHeldoutEnds(i)+1:iHeldoutEnds(i+1)) = iFoldInds(i);
end

%% Iterate through PCs and run analysis
% Set GLM parameters
optsGLM.alpha = glmParams.alpha; % Elastic net hyperparam - somewhere between 0.5 and 0.9
optsGLM.nlambda = glmParams.nLambda; % Number of penalties to try - probably 100;
distr = glmParams.distr;
dmFields = DMdata.fields;
nBaseEnds = [0; cumsum(DMdata.nBases)];

% Movement first:
permPreds = {'MovT'};
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
DMresMov = double(DM(:,iPerm));
iDMresMov = logical(sum(DMresMov,2));

% Reward second:
permPreds = {'RewTcorrect'};
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
DMresRew = double(DM(:,iPerm));
iDMresRew = logical(sum(DMresRew,2));

% Now run it in a for loop:
for iPC = 1:numel(GLMlocal)
    if ~isnan(GLMlocal(iPC).sigBool_shufMov) && sum(DMresMov(:,1)) >= xvalFold
%         yRes = double(GLMlocal(iPC).yData(iDMresMov,:)-GLMlocal(iPC).yHat_shufMov(iDMresMov));
%         % Run GLM - last boolean signifies parallel or non
%         GLMlocal(iPC).GLMobjMov = cvglmnet(DMresMov(iDMresMov,:),yRes,distr,optsGLM,[],[],iFold(iDMresMov),1); % 1 = Parallel, 0 = Serial
%         iBest = find(GLMlocal(iPC).GLMobjMov.lambda == GLMlocal(iPC).GLMobjMov.lambda_min);
%         GLMlocal(iPC).kernelMov = GLMlocal(iPC).GLMobjMov.glmnet_fit.beta(:,iBest);
        
        yRes = double(GLMlocal(iPC).yData-GLMlocal(iPC).yHat_shufMov);
        % Run GLM - last boolean signifies parallel or non
        GLMlocal(iPC).GLMobjMov = cvglmnet(DMresMov,yRes,distr,optsGLM,[],[],iFold,1); % 1 = Parallel, 0 = Serial
        iBest = find(GLMlocal(iPC).GLMobjMov.lambda == GLMlocal(iPC).GLMobjMov.lambda_min);
        GLMlocal(iPC).kernelMov = GLMlocal(iPC).GLMobjMov.glmnet_fit.beta(:,iBest);
        
    end
    
    if ~isnan(GLMlocal(iPC).sigBool_shufRew) && sum(DMresRew(:,1)) >= xvalFold
%         yRes = double(GLMlocal(iPC).yData(iDMresRew)-GLMlocal(iPC).yHat_shufRew(iDMresRew));
%         % Run GLM - last boolean signifies parallel or non
%         GLMlocal(iPC).GLMobjRew = cvglmnet(DMresRew(iDMresRew,:),yRes,distr,optsGLM,[],[],iFold(iDMresRew),1); % 1 = Parallel, 0 = Serial
%         iBest = find(GLMlocal(iPC).GLMobjRew.lambda == GLMlocal(iPC).GLMobjRew.lambda_min);
%         GLMlocal(iPC).kernelRew = GLMlocal(iPC).GLMobjRew.glmnet_fit.beta(:,iBest);
        
        yRes = double(GLMlocal(iPC).yData-GLMlocal(iPC).yHat_shufRew);
        % Run GLM - last boolean signifies parallel or non
        GLMlocal(iPC).GLMobjRew = cvglmnet(DMresRew,yRes,distr,optsGLM,[],[],iFold,1); % 1 = Parallel, 0 = Serial
        iBest = find(GLMlocal(iPC).GLMobjRew.lambda == GLMlocal(iPC).GLMobjRew.lambda_min);
        GLMlocal(iPC).kernelRew = GLMlocal(iPC).GLMobjRew.glmnet_fit.beta(:,iBest);
        
    end
    if mod(iPC,50) == 0 || iPC == numel(GLMlocal)
        tNow = clock;
        fprintf('Done %i cells in %i seconds, time is %i:%02.f on %i-%02.f-%02.f.\n',iPC,round(toc),tNow(4),tNow(5),tNow(1),tNow(2),tNow(3));        
    end

end


end