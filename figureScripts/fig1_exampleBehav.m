function fig1_exampleBehav(dataPaths,dataFold,paperFold,savebool)
v1Path = fullfile(dataFold,'01_Behav+imaging',sprintf('Version_%i',dataPaths(1).version),dataPaths(1).name,...
    sprintf('%s_%s',dataPaths(1).date,dataPaths(1).fov));
v41Path = fullfile(dataFold,'01_Behav+imaging',sprintf('Version_%i',dataPaths(2).version),dataPaths(2).name,...
    sprintf('%s_%s',dataPaths(2).date,dataPaths(2).fov));
v4nPath = fullfile(dataFold,'01_Behav+imaging',sprintf('Version_%i',dataPaths(3).version),dataPaths(3).name,...
    sprintf('%s_%s',dataPaths(3).date,dataPaths(3).fov));
v5Path = fullfile(dataFold,'01_Behav+imaging',sprintf('Version_%i',dataPaths(4).version),dataPaths(4).name,...
    sprintf('%s_%s',dataPaths(4).date,dataPaths(4).fov));

%% DK169 - v1
v1glm = load(fullfile(v1Path,'GLM','glmRfiles','GLMcor_small.mat')); v1glm = v1glm.GLMcor_small;
v1basis = load(fullfile(v1Path,'GLM','basisstructs_cells.mat')); v1expt = v1basis.trialbasis.expt;
exTrials1 = [2,5]; % Real trial num 89-98
figv1 = fig1_exampleBehav_plot(v1glm,v1expt,exTrials1,dataPaths(1));

%% DK169 - v4_1
v41glm = load(fullfile(v41Path,'GLM','glmRfiles','GLMcor_small.mat')); v41glm = v41glm.GLMcor_small;
v41basis = load(fullfile(v41Path,'GLM','basisstructs_cells.mat')); v41expt = v41basis.trialbasis.expt;
exTrials41 = [1:10]; % Real trial num 89-98
figv41 = fig1_exampleBehav_plot(v41glm,v41expt,exTrials41,dataPaths(2));

%% DK169 - v4_n
v4nglm = load(fullfile(v4nPath,'GLM','glmRfiles','GLMcor_small.mat')); v4nglm = v4nglm.GLMcor_small;
v4nbasis = load(fullfile(v4nPath,'GLM','basisstructs_cells.mat')); v4nexpt = v4nbasis.trialbasis.expt;
exTrials4n = [48:57]; % Real trial num 89-98
figv4n = fig1_exampleBehav_plot(v4nglm,v4nexpt,exTrials4n,dataPaths(3));

%% DK169 - v5
v5glm = load(fullfile(v5Path,'GLM','glmRfiles','GLMfull_small.mat')); v5glm = v5glm.GLMfull_small;
v5basis = load(fullfile(v5Path,'GLM','basisstructs_cells.mat')); v5expt = v5basis.trialbasis.expt;
exTrials5 = [71:80]; % Real trial num 89-98
figv5 = fig1_exampleBehav_plot(v5glm,v5expt,exTrials5,dataPaths(4));

%% Save maybe

if savebool
    figFold = fullfile(paperFold,'Figures','Fig1','Sub1_BehavExample');
    if ~exist(figFold,'dir'); mkdir(figFold); end
    
    figname = sprintf('%s_v1_Trials_%i_to_%i',dataPaths(1).name,exTrials1(1),exTrials1(end));
    savefig(figv1, fullfile(figFold,[figname,'.fig']));
    saveas(figv1,fullfile(figFold,[figname,'.png']));
    print(figv1,fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    
    figname = sprintf('%s_v4_1_Trials_%i_to_%i',dataPaths(1).name,exTrials41(1),exTrials41(end));
    savefig(figv41, fullfile(figFold,[figname,'.fig']));
    saveas(figv41,fullfile(figFold,[figname,'.png']));
    print(figv41,fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    
    figname = sprintf('%s_v4_n_Trials_%i_to_%i',dataPaths(1).name,exTrials4n(1),exTrials4n(end));
    savefig(figv4n, fullfile(figFold,[figname,'.fig']));
    saveas(figv4n,fullfile(figFold,[figname,'.png']));
    print(figv4n,fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    
    figname = sprintf('%s_v5_Trials_%i_to_%i',dataPaths(1).name,exTrials5(1),exTrials5(end));
    savefig(figv5, fullfile(figFold,[figname,'.fig']));
    saveas(figv5,fullfile(figFold,[figname,'.png']));
    print(figv5,fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
end

end