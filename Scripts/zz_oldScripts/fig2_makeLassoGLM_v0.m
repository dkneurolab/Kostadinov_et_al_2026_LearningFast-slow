function GLMoutput = fig2_makeLassoGLM_v0(DM,Ydata,iFold,iEnds,nShuf,binBool)
%% Function information:
% DM = design matrix (can be PCA'd or not), shape = M (data points) x N predictors
% Ydata is full data matrix, shape = M (data points) x P neurons
% iFold is folds of cross-validation for grouping, shape = M (data point) x 1,
%    with values going from 1:CV folds + 1, CV+1 fold is holdout data
rng(2349);
%% Make CV partitions, etc.
% Make array of trial indices
iTrial = zeros(size(DM,1),1);
for i = 1:numel(iEnds)-1
    iTrial(iEnds(i)+1:iEnds(i+1)) = i;
end
% Split into test-train+holdout
iHoldout = iFold == max(iFold);
iTestTrain = iFold ~= max(iFold);

% Parition DM
DM = double(DM);
DMtt = DM(iTestTrain,:);
DMholdout = DM(iHoldout,:);
DMshifts = randi(size(DMtt,1),nShuf,1);

% Partition Ydata, use bernoulli or exponential glm
Ydata = double(Ydata);
if binBool
    Ydata = double(logical(Ydata));
    distr = 'binomial';
    link = 'logit';
else
    distr = 'poisson';
    link = 'log';
end
Ydata_tt = Ydata(iTestTrain,:);
Ydata_holdout = Ydata(iHoldout,:);

% Make cv partition object
iFold_tt = iFold(iTestTrain);
nCV = max(iFold_tt);
% cvPart = cvpartition(iFold_tt,'KFold',nCV);

Ymfr = mean(Ydata,1)*30;

%% Run GLM
GLMoutput = struct;
optsGLM.alpha = 0.9;
optsGLM.nlambda = 100;

tic;
iCell = 0;
% for iPC = 200
for iPC = 1:size(Ydata,2)
% for iPC = 1:20
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
    
    devExp_shuf = zeros(nShuf,1);
    for jShuf = 1:nShuf
        DMshuf = circshift(DMtt,DMshifts(jShuf));
%         GLMshuf = cvglmnet(DMshuf,Ydata_tt(:,iPC),distr,optsGLM,[],[],iFold_tt,0);
        GLMshuf = cvglmnet(DMshuf,Ydata_tt(:,iPC),distr,optsGLM,[],[],circshift(iFold_tt,DMshifts(jShuf)),1);
        CVpred_shuf = cvglmnetPredict(GLMshuf, DMholdout, 'lambda_min');
        yHat_shuf = exp(CVpred_shuf);
        devExp_shuf(jShuf) = poissonDevExp(Ydata_holdout(:,iPC),yHat_shuf);
    end
    GLMoutput(iCell).devExp_shuf = devExp_shuf;
    GLMoutput(iCell).sigBool = GLMoutput(iCell).devExp > (mean(devExp_shuf)+2*std(devExp_shuf)) & GLMoutput(iCell).devExp > 0; 
    
    if mod(iPC,10) == 0 || iPC == size(Ydata,2)
        fprintf('Done %i cells in %i seconds.\n',iPC,round(toc));
    end
    
    end
end

end
    
    
    
    
    
    
    
    
    
%     % Use bernoulli or exponential glm
%     [b, binfo] = lassoglm(DMtt,Ydata_tt(:,iPC),distr,'Alpha',0.9,'CV',cvPart);
%     
%     
%     GLMoutput(iCell).fitstats.distr = distr;
%     GLMoutput(iCell).fitstats.link = link;
%     GLMoutput(iCell).fitstats.b = b;
%     GLMoutput(iCell).fitstats.binfo = binfo;
%     b0 = binfo.Intercept(binfo.IndexMinDeviance);    
%     GLMoutput(iCell).Bcoefs = [b0; b(:,binfo.IndexMinDeviance)]; clear b0;
%     
%     yhat_holdout = glmval(GLMoutput(iCell).Bcoefs,DMholdout,link);
%     
%     % % Visualize if you want to:
%     % lassoPlot(b,binfo,'plottype','CV');
%     % legend('show') % Show legend
% 
%     
%     % Make permuted data for shuffling and computing significance of full model   
%     Bshuf = cell(nShuf,1);
%     for jShuf = 1:nShuf
%         DMshuf = circshift(DMtt,DMshifts(jShuf));
%         if binBool
%             [b, binfo] = lassoglm(DMshuf,Ydata_tt(:,iPC),'binomial','Alpha',0.9,'CV',cvPart);
%             GLMoutput(iCell).link = 'binomial';
%         else
%             [b, binfo] = lassoglm(DMshuf,Ydata_tt(:,iPC),'poisson','Alpha',0.9,'CV',cvPart);
%             GLMoutput(iCell).link = 'poisson';
%         end
%         b0 = binfo.Intercept(binfo.IndexMinDeviance);    
%         Bshuf{jShuf} = [b0; b(:,binfo.IndexMinDeviance)]; clear b0;
%     end


% end
% 
% 
% %% Test GLM
% % Use bernoulli or exponential glm
% if binBool
%     yhat_holdout = glmval(Bcoefs,DMholdout,'logit');
% else
%     yhat_holdout = glmval(coef,DMholdout,'log');
% end
% 
% 
% 
% tTrials = diff(iEnds);
% tTrials_holdout = tTrials(unique(iTrial(iHoldout)));
% iEnd_holdout = [0; cumsum(tTrials_holdout)];
% 
% Yon_holdout = zeros(30,numel(tTrials_holdout));
% Yon_hat = Yon_holdout;
% for iHold = 1:numel(tTrials_holdout)
%     Yon_holdout(:,iHold) = Ydata_holdout(iEnd_holdout(iHold)+1:iEnd_holdout(iHold)+30, iPC);
%     Yon_hat(:,iHold) = yhat_holdout(iEnd_holdout(iHold)+1:iEnd_holdout(iHold)+30);
% end
% 
% figure;
% subplot(2,1,1); hold on;
% plot(Ydata_holdout(:,iPC),'-k','LineWidth',1.5);
% plot(yhat_holdout,'-r','LineWidth',1.5);
% subplot(2,1,2); hold on;
% plot(mean(Yon_holdout,2),'-k','LineWidth',1.5);
% plot(mean(Yon_hat,2),'-r','LineWidth',1.5);
% 
% end
