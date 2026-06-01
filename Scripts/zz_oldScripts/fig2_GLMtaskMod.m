function fig2_GLMtaskMod(sessions,dataFold,paperFold,savebool)

for i = 1:numel(sessions)
    glmFold = fullfile(dataFold,'01_Behav+imaging',sprintf('Version_%i',sessions(i).version),sessions(i).name,sprintf('%s_%s',sessions(i).date,sessions(i).fov),'GLM');
%     glmData = load(fullfile(glmFold,'GLMdata_only_cells.mat'));
    DM_spikes = load(fullfile(glmFold,'DM+spikes_cells.mat'));
    setFiles = load(fullfile(glmFold,'setfiles_cells.mat'));
    basisStructs = load(fullfile(glmFold,'basisstructs_cells.mat'));
    
    DMboxdata = DM_spikes.DMboxdata;
    trialBasis = basisStructs.trialbasis;
%     rewBasis = basisStructs.rewbasis;
    xvalFold = setFiles.inits.xvalfold;
    trialOutcomes = setFiles.trialout;
    roiNames = setFiles.mzones;
    
    clear DM_spikes basisStructs
    
    % Make design matrices for test, train, and withold data
    figFold = fullfile(dataFold,'paperFigs','Fig2','Sub1_GLM',sprintf('%s_%s_v%i_%s',sessions(i).name,sessions(i).fov,sessions(i).version,sessions(i).date(end-4:end)));
    figFold2 = fullfile(paperFold,'Figures','Fig2','Sub1_GLM',sprintf('%s_%s_v%i_%s',sessions(i).name,sessions(i).fov,sessions(i).version,sessions(i).date(end-4:end)));
    if ~exist(figFold,'dir'); mkdir (figFold); end
    if ~exist(fullfile(figFold,sprintf('GLMdata_only_%s.mat',setFiles.inits.cellgrping)),'file')
        [GLMdata, GLMholdout,DMdataBoxSmall,trialBasisC] = fig2_make_GLMdata_dk(DMboxdata, trialBasis, roiNames, xvalFold,trialOutcomes);
        if savebool; save(fullfile(figFold,sprintf('GLMdata_only_%s.mat',setFiles.inits.cellgrping)),'GLMdata', 'GLMholdout','DMdataBoxSmall','trialBasisC'); end
    else
        GLM_validation = load(fullfile(figFold,sprintf('GLMdata_only_%s.mat',setFiles.inits.cellgrping)));
        GLMdata = GLM_validation.GLMdata;
        GLMholdout = GLM_validation.GLMholdout;
        DMdataBoxSmall = GLM_validation.DMdataBoxSmall;
        trialBasisC = GLM_validation.trialBasisC;
        clear GLM_validation;
    end
    
    % Get best lambda for each PC dend (tuning lambda)
    [EVbox_bestmu, bestBoxLi] = fig2_bestLam(GLMdata, DMdataBoxSmall, setFiles, figFold, figFold2, savebool);
    clear GLMdata
    
    % Plot neurons and calculate EV
    setFiles.inits.resptime = setFiles.setsBehav.responsetime;
    sessions(i).statStruct = fig2_plotHoldout(trialBasisC,roiNames,lambda(bestBoxLi),GLMholdout,DMdataBoxSmall,setFiles.inits,figFold,figFold2,savebool);
    close all
    
    
    
    if ~exist(fullfile(figFold,sprintf('GLMdata_1perc_%s.mat',setFiles.inits.cellgrping)),'file')
        [GLMdata_1perc, GLMholdout_1perc,DMboxPCA,trialBasisC_1perc] = fig2_make_GLMdata_1perc_dk(DMboxdata, trialBasis, roiNames, xvalFold,trialOutcomes);
        if savebool; save(fullfile(figFold,sprintf('GLMdata_1perc_%s.mat',setFiles.inits.cellgrping)),'GLMdata_1perc', 'GLMholdout_1perc','DMboxPCA','trialBasisC_1perc'); end
    else
        GLM_validation = load(fullfile(figFold,sprintf('GLMdata_1perc_%s.mat',setFiles.inits.cellgrping)));
        GLMdata_1perc = GLM_validation.GLMdata_1perc;
        GLMholdout_1perc = GLM_validation.GLMholdout_1perc;
        DMboxPCA = GLM_validation.DMboxPCA;
        trialBasisC_1perc = GLM_validation.trialBasisC_1perc;
        clear GLM_validation;
    end
    [EVbox_bestmu_1perc, bestBoxLi_1perc] = fig2_bestLam(GLMdata_1perc, DMdataBoxSmall, setFiles, figFold, figFold2, savebool);
    clear GLMdata
    
    
    
    close all;
end

glmSessions = sessions;
figFoldAll = fullfile(dataFold,'paperFigs','Fig2','Sub1_GLM');
figFold2All = fullfile(paperFold,'Figures','Fig2','Sub1_GLM');

fig2_sigFitPerFOV(glmSessions,figFoldAll,figFold2All,savebool)

if savebool
    save(fullfile(figFoldAll,'seshStats.mat'),'glmSessions');
    save(fullfile(figFold2All,'seshStats.mat'),'glmSessions');
end

end