function fig3_GLMtaskMod_v4n(sessions,pcaVar, dataFold, paperFold, savebool)

for i = 1:numel(sessions)
    % Load in data needed to make new smaller design matrix
    glmFold = fullfile(dataFold,'01_Behav+imaging',sprintf('Version_%i',sessions(i).version),sessions(i).name,sprintf('%s_%s',sessions(i).date,sessions(i).fov),'GLM');
%     DM_spikes = load(fullfile(glmFold,'DM+spikes_cells.mat'));
    setFiles = load(fullfile(glmFold,'setfiles_cells.mat'));
    basisStructs = load(fullfile(glmFold,'basisstructs_cells.mat'));
    
    % Extract the data that will be fed into script to make small design matrix
%     DMboxdata = DM_spikes.DMboxdata;
    trialBasis = basisStructs.trialbasis;
%     xvalFold = setFiles.inits.xvalfold;
    trialOutcomes = setFiles.trialout;    
    trialBasis = rmfield(trialBasis,'bs');
    setFiles.inits.imfs = setFiles.setsImg.imfs;
    setFiles.inits.vrfs = setFiles.setsImg.vrfs;
    setFiles.inits.paqfsfs = setFiles.setsImg.fs;
    setFiles.inits.pcaVar = pcaVar;
    clear DM_spikes basisStructs
    
    % Make small design matrix and save 
    [GLMfull, GLMfull_small, GLMcor, GLMcor_small] = fig2_GLMtaskMod_v4n_makeDM(trialBasis,setFiles.inits,trialOutcomes,glmFold,savebool);    
end
   