function fig2_extractSigs(sessions,dataFold,paperFold,glmParams,savebool)
%%
glmFold = sprintf('GLMfull_%s_%s_alpha%i',glmParams.dataType,glmParams.distr,glmParams.alpha*100);
figFold0 = fullfile(paperFold,'Figures','Fig2','Sub3_SummaryBars');
figFold = fullfile(figFold0,sprintf('%s_%s_alpha%i',glmParams.dataType,glmParams.distr,glmParams.alpha*100));
if ~exist(figFold,'dir'); mkdir(figFold); end

%%

glmSummary = struct;

for iSesh = 1:numel(sessions)
    dataPath = fullfile(dataFold,sprintf('Version_%i',sessions(iSesh).version),sessions(iSesh).name,...
        sprintf('%s_%s',sessions(iSesh).date,sessions(iSesh).fov),'GLM','glmRfiles');
    glmLocal = load(fullfile(dataPath,glmFold,'GLMsig+VisMovRewLick.mat')); glmLocal = glmLocal.GLMoutput;
    glmSummary(iSesh).name = sessions(iSesh).name;
    glmSummary(iSesh).date = sessions(iSesh).date;
    glmSummary(iSesh).fov = sessions(iSesh).fov;
    glmSummary(iSesh).iLobV = strcmpi(glmSummary(iSesh).fov,'lobv');
    glmSummary(iSesh).iSim2 = strcmpi(glmSummary(iSesh).fov,'sim2');
    
    for jLocal = 1:numel(glmLocal)
        if ~isempty(glmLocal(jLocal).sigBool)
            glmSummary(iSesh).sigBool(jLocal,1) = glmLocal(jLocal).sigBool;
        end
        if ~isempty(glmLocal(jLocal).sigBool_shufVis_proj)
            glmSummary(iSesh).sigBoolVis(jLocal,1) = glmLocal(jLocal).sigBool_shufVis_proj;
        end
        if ~isempty(glmLocal(jLocal).sigBool_shufMov_proj)
            glmSummary(iSesh).sigBoolMov(jLocal,1) = glmLocal(jLocal).sigBool_shufMov_proj;
        end
        if ~isempty(glmLocal(jLocal).sigBool_shufRew_proj)
            glmSummary(iSesh).sigBoolRew(jLocal,1) = glmLocal(jLocal).sigBool_shufRew_proj;
        end
        if ~isempty(glmLocal(jLocal).sigBool_shufLick_proj)
            glmSummary(iSesh).sigBoolLick(jLocal,1) = glmLocal(jLocal).sigBool_shufLick_proj;
        end
    end
end

if savebool
    save(fullfile(figFold,'sigSummary.mat'),'glmSummary');
end

end