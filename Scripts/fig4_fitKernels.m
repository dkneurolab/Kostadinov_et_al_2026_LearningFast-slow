function fig4_fitKernels(sessionsv4_n,sessionsv1,sessionsv4_1,glmParams,dataFold,paperFold,savebool)
%% Load in GLM files for training


%% Fit kernels on V4 sessions if it does not exist
figFold = fullfile(paperFold,'Figures','Fig4','Sub5_fitKernels');

if ~exist(fullfile(figFold,'kernelsGLMv144.mat'),'file')

matchFold = fullfile(paperFold,'Figures','Fig4','Sub2_CrossDayPredictions');
GLMv1v4n = load(fullfile(matchFold,sprintf('GLM_v1_v4n_corrects.mat'))); GLMv1v4n = GLMv1v4n.GLMout;
GLMv41v4n = load(fullfile(matchFold,sprintf('GLM_v41_v4n_corrects.mat'))); GLMv41v4n = GLMv41v4n.GLMout;

for iSesh = 1:numel(sessionsv4_n)
    tSesh = tic;
    % Process GLM correct trials only in V4_n sessions using cvglmnet, etc.
    glmFold = fullfile(dataFold,'01_Behav+imaging',sprintf('Version_%i',sessionsv4_n(iSesh).version),sessionsv4_n(iSesh).name,...
        sprintf('%s_%s',sessionsv4_n(iSesh).date,sessionsv4_n(iSesh).fov),'GLM','glmRfiles');
    GLMoutput = fig4_fitV4corr_kernels(glmFold,glmParams);
    fprintf('Session took %i seconds.\n\n',round(toc(tSesh)));
    if savebool
        savFold = fullfile(figFold,'v4Fits');
        if ~exist(savFold,'dir'); mkdir(savFold); end
        save(fullfile(savFold,sprintf('%s_%s.mat',sessionsv4_n(iSesh).fov,sessionsv4_n(iSesh).name)),'GLMoutput');
    end
    
    % Find matching dataset in v1 data
    for j = 1:numel(GLMv1v4n)
        if strcmpi(GLMv1v4n(j).name,sessionsv4_n(iSesh).name)
            if (GLMv1v4n(j).ilobv && strcmpi(sessionsv4_n(iSesh).fov,'lobv')) || ...
                   (GLMv1v4n(j).isim2 && strcmpi(sessionsv4_n(iSesh).fov,'sim2')) 
                GLMv1v4n(j).GLMv4 = GLMoutput(cat(1,GLMv1v4n(j).GLMstruct.iPC_v4));
%             elseif GLMv1v4n(j).isim2 && strcmpi(sessionsv4_n(iSesh).fov,'sim2')
%                 GLMv1v4n(j).GLMv4 = GLMoutput(cat(1,GLMv1v4n(j).GLMstruct.iPC_v4));
            end
        end
    end
    % Find matching dataset in v41 data
    for j = 1:numel(GLMv41v4n)
        if strcmpi(GLMv41v4n(j).name,sessionsv4_n(iSesh).name)
            if (GLMv41v4n(j).ilobv && strcmpi(sessionsv4_n(iSesh).fov,'lobv')) || ...
                    (GLMv41v4n(j).isim2 && strcmpi(sessionsv4_n(iSesh).fov,'sim2'))
                GLMv41v4n(j).GLMv4 = GLMoutput(cat(1,GLMv41v4n(j).GLMstruct.iPC_v4));
%             elseif (GLMv41v4n(j).isim2 && strcmpi(sessionsv4_n(iSesh).fov,'sim2'))
%                 GLMv41v4n(j).GLMv4 = GLMoutput(cat(1,GLMv41v4n(j).GLMstruct.iPC_v4));
            end
        end
    end
end

%% Process v1 sessions
for iSesh = 1:numel(GLMv1v4n)
    glmFold = fullfile(dataFold,'01_Behav+imaging','Version_1',sessionsv1(iSesh).name,...
        sprintf('%s_%s',sessionsv1(iSesh).date,sessionsv1(iSesh).fov),'GLM','glmRfiles');
    glmDatav1 = load(fullfile(glmFold,'GLMcor_small.mat')); glmDatav1 = glmDatav1.GLMcor_small;
    GLMout = fig4_fitKernels_crossDay(GLMv1v4n(iSesh),glmDatav1,glmParams);
    GLMv1v4n(iSesh).GLMv1 = GLMout;
end

for iSesh = 1:numel(GLMv41v4n)
    glmFold = fullfile(dataFold,'01_Behav+imaging','Version_4',sessionsv1(iSesh).name,...
        sprintf('%s_%s',sessionsv4_1(iSesh).date,sessionsv4_1(iSesh).fov),'GLM','glmRfiles');
    glmDatav4 = load(fullfile(glmFold,'GLMcor_small.mat')); glmDatav4 = glmDatav4.GLMcor_small;
    GLMout = fig4_fitKernels_crossDay(GLMv41v4n(iSesh),glmDatav4,glmParams);
    GLMv41v4n(iSesh).GLMv41 = GLMout;
end

%% Combine in a single data structure
glmFields = fieldnames(GLMv41v4n);
GLMv144 = rmfield(GLMv41v4n,glmFields(4:end));

for iSesh = 1:numel(GLMv41v4n)
    [iPCv1,iPCv41] = ismember(cat(1,GLMv1v4n(iSesh).GLMstruct.iPC_v4),cat(1,GLMv41v4n(iSesh).GLMstruct.iPC_v4));
    
    iPCv1 = find(iPCv1);
    iPCv41 = iPCv41(iPCv41>0);
    GLMv144(iSesh).GLMmatchv4n_v1 = GLMv1v4n(iSesh).GLMstruct(iPCv1);
    GLMv144(iSesh).GLMmatchv4n_v41 = GLMv41v4n(iSesh).GLMstruct(iPCv41);
    GLMv144(iSesh).GLMkernelsv1 = GLMv1v4n(iSesh).GLMv1(iPCv1);
    GLMv144(iSesh).GLMkernelsv41 = GLMv41v4n(iSesh).GLMv41(iPCv41);
    GLMv144(iSesh).GLMkernelsv4n = GLMv1v4n(iSesh).GLMv4(iPCv1);
    
    % Get the relevant cell indices and extract movement kernels
    iMovV4 = cat(1,GLMv144(iSesh).GLMkernelsv4n.sigBoolMov);
    v4ntemp = GLMv144(iSesh).GLMkernelsv4n(iMovV4);
    v41temp = GLMv144(iSesh).GLMkernelsv41(iMovV4);
    v1temp = GLMv144(iSesh).GLMkernelsv1(iMovV4);
    GLMv144(iSesh).movKernel_v4n = cat(2,v4ntemp.kernelMov);
    if isfield(v41temp,'kernelMov')
        GLMv144(iSesh).movKernel_v41 = cat(2,v41temp.kernelMov);
    else
        GLMv144(iSesh).movKernel_v41 = nan(size(GLMv144(iSesh).movKernel_v4n));
    end
    if isfield(v1temp,'kernelMov')
        GLMv144(iSesh).movKernel_v1 = cat(2,v1temp.kernelMov);
    else
        GLMv144(iSesh).movKernel_v1 = nan(size(GLMv144(iSesh).movKernel_v4n));
    end
    
    % Get the relevant cell indices and extract movement kernels
    iRewV4 = cat(1,GLMv144(iSesh).GLMkernelsv4n.sigBoolRew);
    v4ntemp = GLMv144(iSesh).GLMkernelsv4n(iRewV4);
    v41temp = GLMv144(iSesh).GLMkernelsv41(iRewV4);
    v1temp = GLMv144(iSesh).GLMkernelsv1(iRewV4);
    GLMv144(iSesh).rewKernel_v4n = cat(2,v4ntemp.kernelRew);
    if isfield(v41temp,'kernelRew')
        GLMv144(iSesh).rewKernel_v41 = cat(2,v41temp.kernelRew);
    else
        GLMv144(iSesh).rewKernel_v41 = nan(size(GLMv144(iSesh).rewKernel_v4n));
    end
    if isfield(v1temp,'kernelRew')
        GLMv144(iSesh).rewKernel_v1 = cat(2,v1temp.kernelRew);
    else
        GLMv144(iSesh).rewKernel_v1 = nan(size(GLMv144(iSesh).rewKernel_v4n));
    end

end

else
    GLMv144 = load(fullfile(figFold,'kernelsGLMv144.mat'));
    GLMv144 = GLMv144.GLMv144;

end



%% Make plots
%split by FOV
GLMv144_lobv = GLMv144(cat(1,GLMv144.ilobv));
GLMv144_sim2 = GLMv144(cat(1,GLMv144.isim2));

[figMovLobv, statMovLobv] = fig4_fitKernels_plot(GLMv144_lobv,'mov');
[figMovSim2, statMovSim2] = fig4_fitKernels_plot(GLMv144_sim2,'mov');
[figRewLobv, statRewLobv] = fig4_fitKernels_plot(GLMv144_lobv,'rew');
[figRewSim2, statRewSim2] = fig4_fitKernels_plot(GLMv144_sim2,'rew');

statsFriedman = v2struct(statMovLobv,statMovSim2,statRewLobv,statRewSim2);

if savebool
    if ~exist(figFold,'dir'); mkdir(figFold); end

    figname = 'kernelMatch_lobv_mov';
    savefig(figMovLobv, fullfile(figFold,[figname,'.fig']));
    saveas(figMovLobv,fullfile(figFold,[figname,'.png']));
    print(figMovLobv,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');
    
    figname = 'kernelMatch_sim2_mov';
    savefig(figMovSim2, fullfile(figFold,[figname,'.fig']));
    saveas(figMovSim2,fullfile(figFold,[figname,'.png']));
    print(figMovSim2,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');
    
    figname = 'kernelMatch_lobv_rew';
    savefig(figRewLobv, fullfile(figFold,[figname,'.fig']));
    saveas(figRewLobv,fullfile(figFold,[figname,'.png']));
    print(figRewLobv,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');
    
    figname = 'kernelMatch_sim2_rew';
    savefig(figRewSim2, fullfile(figFold,[figname,'.fig']));
    saveas(figRewSim2,fullfile(figFold,[figname,'.png']));
    print(figRewSim2,fullfile(figFold,[figname,'.eps']), '-depsc', '-vector');
        
    save(fullfile(figFold,'kernelsGLMv144.mat'),'GLMv144','statsFriedman');
    
end

end