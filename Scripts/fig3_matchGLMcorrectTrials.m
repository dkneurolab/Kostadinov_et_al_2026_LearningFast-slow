function GLMoutput_v1 = fig3_matchGLMcorrectTrials(v4Sesh,v1Sesh,glmParams,dataFold)

v1Date = v1Sesh.date(end-4:end);
v4Date = v4Sesh.date(end-4:end);
strDash = strfind(v1Date,'-'); v1Date(strDash) = [];
strDash = strfind(v4Date,'-'); v4Date(strDash) = [];
regPath = sprintf('REG_%s_%s',v4Date,v1Date);

imFold = fullfile(dataFold,'02_Imaging',v4Sesh.name,regPath);
matchFile = dir(fullfile(imFold,'*analysis.mat'));
matchFile = load(fullfile(imFold,matchFile.name));
iMatch0 = matchFile.regi.rois.idcs;
[~,iSort] = sort(iMatch0(:,1),'ascend');
iMatch = iMatch0(iSort,:);

glmFold1 = fullfile(dataFold,'01_Behav+imaging',sprintf('Version_%i',v1Sesh.version),v1Sesh.name,...
    sprintf('%s_%s',v1Sesh.date,v1Sesh.fov),'GLM','glmRfiles');
glmFold4 = fullfile(dataFold,'01_Behav+imaging',sprintf('Version_%i',v4Sesh.version),v4Sesh.name,...
    sprintf('%s_%s',v4Sesh.date,v4Sesh.fov),'GLM','glmRfiles');
glmFold4_model = fullfile(glmFold4,sprintf('GLMfull_%s_%s_alpha%i',glmParams.dataType,glmParams.distr,glmParams.alpha*100));

%%
glmV1 = load(fullfile(glmFold1,'GLMcor_small.mat')); glmV1 = glmV1.GLMcor_small;
glmV4 = load(fullfile(glmFold4,'GLMcor_small.mat')); glmV4 = glmV4.GLMcor_small;
GLMoutput_v4 = load(fullfile(glmFold4_model,'GLMsigCor+VisMovRewLick.mat')); GLMoutput_v4 = GLMoutput_v4.GLMoutput;
GLMoutput_v4 = GLMoutput_v4(iMatch(:,1));

DMv1 = glmV1.DM;
yDatav1 = glmV1.Ytrials(:,iMatch(:,2));

DMv1_project = (DMv1-mean(DMv1))*glmV4.DMpca.coef;

nShuf = numel(GLMoutput_v4(1).GLMobj_shuf);

%%
GLMoutput_v1 = struct;
for iRoi = 1:size(iMatch,1)
    GLMoutput_v1(iRoi).iPC = iMatch(iRoi,2);
    GLMoutput_v1(iRoi).iPC_v4 = iMatch(iRoi,1);
    GLMoutput_v1(iRoi).GLMobj_v4 = GLMoutput_v4(iRoi).GLMobj;
    GLMoutput_v1(iRoi).GLMrankPCA = GLMoutput_v4(iRoi).GLMrankPCA;
    GLMoutput_v1(iRoi).yData = yDatav1(:,iRoi);
    GLMoutput_v1(iRoi).sigV4 = false;
    DMv1_local = DMv1_project(:,1:GLMoutput_v1(iRoi).GLMrankPCA);
    
    if ~isempty(GLMoutput_v1(iRoi).GLMobj_v4)
        GLMoutput_v1(iRoi).sigV4 = GLMoutput_v4(iRoi).sigBool;
        
        % Reduced model
        CVpred = cvglmnetPredict(GLMoutput_v1(iRoi).GLMobj_v4, DMv1_local, 'lambda_min');
        if strcmpi(glmParams.distr,'gaussian')
            GLMoutput_v1(iRoi).yHat = CVpred;
            GLMoutput_v1(iRoi).devExp = gaussianDevExp(GLMoutput_v1(iRoi).yData,GLMoutput_v1(iRoi).yHat);
        elseif strcmpi(glmParams.distr,'poisson')
            GLMoutput_v1(iRoi).yHat = exp(CVpred);
            GLMoutput_v1(iRoi).devExp = poissonDevExp(GLMoutput_v1(iRoi).yData,GLMoutput_v1(iRoi).yHat);
        end
        
        % Full model
        CVpred_full = cvglmnetPredict(GLMoutput_v4(iRoi).GLMobj_full, DMv1, 'lambda_min');
        if strcmpi(glmParams.distr,'gaussian')
            GLMoutput_v1(iRoi).yHat_full = CVpred_full;
            GLMoutput_v1(iRoi).devExp_full = gaussianDevExp(GLMoutput_v1(iRoi).yData,GLMoutput_v1(iRoi).yHat_full);
        elseif strcmpi(glmParams.distr,'poisson')
            GLMoutput_v1(iRoi).yHat_full = exp(CVpred_full);
            GLMoutput_v1(iRoi).devExp_full = poissonDevExp(GLMoutput_v1(iRoi).yData,GLMoutput_v1(iRoi).yHat_full);
        end
        
        % Shuf projection - significant fit?
        if ~isempty(GLMoutput_v4(iRoi).GLMobj_shuf)
            devExp_shuf = zeros(nShuf,1);
            yHat_shuf = zeros(numel(GLMoutput_v1(iRoi).yData),nShuf);
            for jShuf = 1:nShuf
                CVpred_shuf = cvglmnetPredict(GLMoutput_v4(iRoi).GLMobj_shuf{jShuf,1}, DMv1_local, 'lambda_min');
                if strcmpi(glmParams.distr,'gaussian')
                    yHat_shuf(:,jShuf) = CVpred_shuf;
                    devExp_shuf(jShuf) = gaussianDevExp(GLMoutput_v1(iRoi).yData,yHat_shuf(:,jShuf));
                elseif strcmpi(glmParams.distr,'poisson')
                    yHat_shuf(:,jShuf) = exp(CVpred_shuf);
                    % Calculate real deviance explained
                    devExp_shuf(jShuf) = poissonDevExp(GLMoutput_v1(iRoi).yData,yHat_shuf(:,jShuf));
                end
            end
            GLMoutput_v1(iRoi).devExp_shuf = devExp_shuf;
            GLMoutput_v1(iRoi).yHat_shuf = single(mean(yHat_shuf,2));
            GLMoutput_v1(iRoi).sigBool = GLMoutput_v1(iRoi).devExp > (mean(devExp_shuf)+2*std(devExp_shuf)) & GLMoutput_v1(iRoi).devExp > 0;
            
            % Vis projection - significant fit?
            if GLMoutput_v4(iRoi).sigBool_shufVis_proj
                devExp_shuf = zeros(nShuf,1);
                yHat_shuf = zeros(numel(GLMoutput_v1(iRoi).yData),nShuf);
                for jShuf = 1:nShuf
                    CVpred_shuf = cvglmnetPredict(GLMoutput_v4(iRoi).GLMobj_shufVis_proj{jShuf,1}, DMv1_local, 'lambda_min');
                    if strcmpi(glmParams.distr,'gaussian')
                        yHat_shuf(:,jShuf) = CVpred_shuf;
                        devExp_shuf(jShuf) = gaussianDevExp(GLMoutput_v1(iRoi).yData,yHat_shuf(:,jShuf));
                    elseif strcmpi(glmParams.distr,'poisson')
                        yHat_shuf(:,jShuf) = exp(CVpred_shuf);
                        % Calculate real deviance explained
                        devExp_shuf(jShuf) = poissonDevExp(GLMoutput_v1(iRoi).yData,yHat_shuf(:,jShuf));
                    end
                end
                GLMoutput_v1(iRoi).devExp_shufVis = devExp_shuf;
                GLMoutput_v1(iRoi).yHat_shufVis = single(mean(yHat_shuf,2));
                GLMoutput_v1(iRoi).sigBool_shufVis = GLMoutput_v1(iRoi).devExp > (mean(devExp_shuf)+2*std(devExp_shuf)) & GLMoutput_v1(iRoi).devExp > 0;
            else
                GLMoutput_v1(iRoi).devExp_shufVis = [];
                GLMoutput_v1(iRoi).yHat_shufVis = [];
                GLMoutput_v1(iRoi).sigBool_shufVis = nan;
            end
            
            % Mov projection - significant fit?
            if GLMoutput_v4(iRoi).sigBool_shufMov_proj
                devExp_shuf = zeros(nShuf,1);
                yHat_shuf = zeros(numel(GLMoutput_v1(iRoi).yData),nShuf);
                for jShuf = 1:nShuf
                    CVpred_shuf = cvglmnetPredict(GLMoutput_v4(iRoi).GLMobj_shufMov_proj{jShuf,1}, DMv1_local, 'lambda_min');
                    if strcmpi(glmParams.distr,'gaussian')
                        yHat_shuf(:,jShuf) = CVpred_shuf;
                        devExp_shuf(jShuf) = gaussianDevExp(GLMoutput_v1(iRoi).yData,yHat_shuf(:,jShuf));
                    elseif strcmpi(glmParams.distr,'poisson')
                        yHat_shuf(:,jShuf) = exp(CVpred_shuf);
                        % Calculate real deviance explained
                        devExp_shuf(jShuf) = poissonDevExp(GLMoutput_v1(iRoi).yData,yHat_shuf(:,jShuf));
                    end
                end
                GLMoutput_v1(iRoi).devExp_shuMov = devExp_shuf;
                GLMoutput_v1(iRoi).yHat_shufMov = single(mean(yHat_shuf,2));
                GLMoutput_v1(iRoi).sigBool_shufMov = GLMoutput_v1(iRoi).devExp > (mean(devExp_shuf)+2*std(devExp_shuf)) & GLMoutput_v1(iRoi).devExp > 0;
            else
                GLMoutput_v1(iRoi).devExp_shuMov = [];
                GLMoutput_v1(iRoi).yHat_shufMov = [];
                GLMoutput_v1(iRoi).sigBool_shufMov = nan;
            end
            
            % Rew projection - significant fit?
            if GLMoutput_v4(iRoi).sigBool_shufRew_proj
                devExp_shuf = zeros(nShuf,1);
                yHat_shuf = zeros(numel(GLMoutput_v1(iRoi).yData),nShuf);
                for jShuf = 1:nShuf
                    CVpred_shuf = cvglmnetPredict(GLMoutput_v4(iRoi).GLMobj_shufRew_proj{jShuf,1}, DMv1_local, 'lambda_min');
                    if strcmpi(glmParams.distr,'gaussian')
                        yHat_shuf(:,jShuf) = CVpred_shuf;
                        devExp_shuf(jShuf) = gaussianDevExp(GLMoutput_v1(iRoi).yData,yHat_shuf(:,jShuf));
                    elseif strcmpi(glmParams.distr,'poisson')
                        yHat_shuf(:,jShuf) = exp(CVpred_shuf);
                        % Calculate real deviance explained
                        devExp_shuf(jShuf) = poissonDevExp(GLMoutput_v1(iRoi).yData,yHat_shuf(:,jShuf));
                    end
                end
                GLMoutput_v1(iRoi).devExp_shufRew = devExp_shuf;
                GLMoutput_v1(iRoi).yHat_shufRew = single(mean(yHat_shuf,2));
                GLMoutput_v1(iRoi).sigBool_shufRew = GLMoutput_v1(iRoi).devExp > (mean(devExp_shuf)+2*std(devExp_shuf)) & GLMoutput_v1(iRoi).devExp > 0;
            else
                GLMoutput_v1(iRoi).devExp_shufRew = [];
                GLMoutput_v1(iRoi).yHat_shufRew = [];
                GLMoutput_v1(iRoi).sigBool_shufRew = nan;
            end
            
            % Lick projection - significant fit?
            if GLMoutput_v4(iRoi).sigBool_shufLick_proj
                devExp_shuf = zeros(nShuf,1);
                yHat_shuf = zeros(numel(GLMoutput_v1(iRoi).yData),nShuf);
                for jShuf = 1:nShuf
                    CVpred_shuf = cvglmnetPredict(GLMoutput_v4(iRoi).GLMobj_shufLick_proj{jShuf,1}, DMv1_local, 'lambda_min');
                    if strcmpi(glmParams.distr,'gaussian')
                        yHat_shuf(:,jShuf) = CVpred_shuf;
                        devExp_shuf(jShuf) = gaussianDevExp(GLMoutput_v1(iRoi).yData,yHat_shuf(:,jShuf));
                    elseif strcmpi(glmParams.distr,'poisson')
                        yHat_shuf(:,jShuf) = exp(CVpred_shuf);
                        % Calculate real deviance explained
                        devExp_shuf(jShuf) = poissonDevExp(GLMoutput_v1(iRoi).yData,yHat_shuf(:,jShuf));
                    end
                end
                GLMoutput_v1(iRoi).devExp_shufLick = devExp_shuf;
                GLMoutput_v1(iRoi).yHat_shufLick = single(mean(yHat_shuf,2));
                GLMoutput_v1(iRoi).sigBool_shufLick = GLMoutput_v1(iRoi).devExp > (mean(devExp_shuf)+2*std(devExp_shuf)) & GLMoutput_v1(iRoi).devExp > 0;
            else
                GLMoutput_v1(iRoi).devExp_shufLick = [];
                GLMoutput_v1(iRoi).yHat_shufLick = [];
                GLMoutput_v1(iRoi).sigBool_shufLick = nan;
            end
        end
    end
end


end