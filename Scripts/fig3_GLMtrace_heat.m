function fig3_GLMtrace_heat(sessions,dataFold,paperFold,glmParams,tracebool,savebool)
%% Make save dirs
figFold0 = fullfile(paperFold,'Figures','Fig2','Sub2_Traces+Heatmaps');
figFold = fullfile(figFold0,sprintf('%s_%s_alpha%i',glmParams.dataType,glmParams.distr,glmParams.alpha*100));
figFold1 = fullfile(figFold,'Heatmaps');
figFold2 = fullfile(figFold,'Traces');
% if ~exist(figFold,'dir'); mkdir(figFold); end
if ~exist(figFold1,'dir'); mkdir(figFold1); end
if ~exist(figFold2,'dir'); mkdir(figFold2); end

%% Loop through FOVs and process them
for iSesh = 1:numel(sessions)
    % Load in GLM data structures
    dataPath0 = fullfile(dataFold,sprintf('Version_%i',sessions(iSesh).version),sessions(iSesh).name,...
        sprintf('%s_%s',sessions(iSesh).date,sessions(iSesh).fov));
    dataPath = fullfile(dataPath0,'GLM');
    glmFold = fullfile(dataPath,'glmRfiles',sprintf('GLMfull_%s_%s_alpha%i',glmParams.dataType,glmParams.distr,glmParams.alpha*100));
    GLMsmall = load(fullfile(dataPath,'glmRfiles','GLMfull_small.mat')); GLMsmall = GLMsmall.GLMfull_small;
    GLMbasis = load(fullfile(dataPath,'basisstructs_cells.mat'));
    expt = GLMbasis.trialbasis.expt;
    GLMsig0 = load(fullfile(glmFold,'GLMsig+VisMovRewLick.mat')); GLMsig0 = GLMsig0.GLMoutput;
    sigBool = false(numel(GLMsig0),1);
    for iMod = 1:numel(GLMsig0)
        if ~isempty(GLMsig0(iMod).GLMobj) && GLMsig0(iMod).devExp > 0
            sigBool(iMod) = true;
        end
    end
    GLMsig = GLMsig0(sigBool);
    GLMsets = load(fullfile(dataPath,'setfiles_cells.mat'));
    clear trialsAll trialstruct GLMbasis GLMsig0;
    
    % Trim out bad trials from expt!
    trialstruct = load(fullfile(dataPath0,'trialsonly.mat'));
    trialsUCO = sort(cat(1,trialstruct.trialstructs.LUstruct.trialind,trialstruct.trialstructs.LCstruct.trialind,trialstruct.trialstructs.LOstruct.trialind));
    trialsALL = cat(1,expt.trialnum);    
    expt = expt(ismember(trialsALL,trialsUCO));
    
    % Make matrices and align to events
    eventsOut = fig2_findEvents(GLMsmall,GLMsig,expt,GLMsets.inits);
    
    % Get PCs with significant fit
    iGood = eventsOut(1).iCells(logical(cat(1,GLMsig.sigBool)));
    
    eventsOutSig = fig2_findEvents(GLMsmall,GLMsig(iGood),expt,GLMsets.inits);
    
    
    if tracebool
        % Plot traces and maybe save
        fig2_GLM_plotTraces(sessions(iSesh),GLMsig,eventsOut,figFold2,savebool);
    end
    
    % Plot heatmaps and maybe save    
    nCells = numel(eventsOut(1).iCells);
    nCellsSig = numel(eventsOutSig(1).iCells);
    fs = GLMsets.inits.imfs;
    t0 = [fs/5,fs/2,fs/2,fs/2]+.5;
    
    fVis = fig2_GLM_plotHeatmap(sessions(iSesh),eventsOut(1),nCells,t0(1),'Trial onset');
    fVisSig = fig2_GLM_plotHeatmap(sessions(iSesh),eventsOutSig(1),nCellsSig,t0(1),'Trial onset');
    fMov = fig2_GLM_plotHeatmap(sessions(iSesh),eventsOut(2),nCells,t0(2),'Movement onset');
    fMovSig = fig2_GLM_plotHeatmap(sessions(iSesh),eventsOutSig(2),nCellsSig,t0(2),'Movement onset');
    fRew = fig2_GLM_plotHeatmap(sessions(iSesh),eventsOut(3),nCells,t0(3),'Reward');
    fRewSig = fig2_GLM_plotHeatmap(sessions(iSesh),eventsOutSig(3),nCellsSig,t0(3),'Reward');
    fLick = fig2_GLM_plotHeatmap(sessions(iSesh),eventsOut(4),nCells,t0(4),'Lick');
    fLickSig = fig2_GLM_plotHeatmap(sessions(iSesh),eventsOutSig(4),nCellsSig,t0(4),'Lick');
    
    if savebool
        fname = sprintf('%s_%s_eventsOut.mat',sessions(iSesh).name,sessions(iSesh).fov);
        save(fullfile(figFold,fname),'eventsOut','eventsOutSig');
        figName = sprintf('%s_%s_Heatmap',sessions(iSesh).name,sessions(iSesh).fov);

        savefig(fVis,fullfile(figFold1,[figName,'_On.fig']));
        saveas(fVis,fullfile(figFold1,[figName,'_On.png']));
        print(fVis,fullfile(figFold1,[figName,'_On.eps']), '-depsc', '-painters');
        savefig(fVisSig,fullfile(figFold1,[figName,'_OnSig.fig']));
        saveas(fVisSig,fullfile(figFold1,[figName,'_OnSig.png']));
        print(fVisSig,fullfile(figFold1,[figName,'_OnSig.eps']), '-depsc', '-painters');
        
        savefig(fMov,fullfile(figFold1,[figName,'_Mov.fig']));
        saveas(fMov,fullfile(figFold1,[figName,'_Mov.png']));
        print(fMov,fullfile(figFold1,[figName,'_Mov.eps']), '-depsc', '-painters');
        savefig(fMovSig,fullfile(figFold1,[figName,'_MovSig.fig']));
        saveas(fMovSig,fullfile(figFold1,[figName,'_MovSig.png']));
        print(fMovSig,fullfile(figFold1,[figName,'_MovSig.eps']), '-depsc', '-painters');
        
        savefig(fRew,fullfile(figFold1,[figName,'_Rew.fig']));
        saveas(fRew,fullfile(figFold1,[figName,'_Rew.png']));
        print(fRew,fullfile(figFold1,[figName,'_Rew.eps']), '-depsc', '-painters');
        savefig(fRewSig,fullfile(figFold1,[figName,'_RewSig.fig']));
        saveas(fRewSig,fullfile(figFold1,[figName,'_RewSig.png']));
        print(fRewSig,fullfile(figFold1,[figName,'_RewSig.eps']), '-depsc', '-painters');
        
        savefig(fLick,fullfile(figFold1,[figName,'_Lick.fig']));
        saveas(fLick,fullfile(figFold1,[figName,'_Lick.png']));
        print(fLick,fullfile(figFold1,[figName,'_Lick.eps']), '-depsc', '-painters');
        savefig(fLickSig,fullfile(figFold1,[figName,'_LickSig.fig']));
        saveas(fLickSig,fullfile(figFold1,[figName,'_LickSig.png']));
        print(fLickSig,fullfile(figFold1,[figName,'_LickSig.eps']), '-depsc', '-painters');
    end
    
    close all
end

end