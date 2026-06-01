function fig6_part1

%% Load in stuff to correlate v1-4 learning and v4 trial over trial stuff.
%v4
v4Behav = load(fullfile('C:\Users\Dimitar\Dropbox\DK papers\2022_Kostadinov_FastSlowLearning\Figures\Fig3\Sub5_behavCorrelates\spk_gaussian_alpha50\','v4MovRewData.mat'));
v4Behav = v4Behav.seshV4;
v4Param = load(fullfile('C:\Users\Dimitar\Dropbox\DK papers\2022_Kostadinov_FastSlowLearning\Figures\Fig3\Sub5_behavCorrelates\spk_gaussian_alpha50','v4_parametrize.mat'));
v4Param = v4Param.seshSummary;
v4Behav = v4Behav([1,4,6,7,10,11,2,3,5,8,9,12]);

%v144
v144kernels = load(fullfile('C:\Users\Dimitar\Dropbox\DK papers\2022_Kostadinov_FastSlowLearning\Figures\Fig4\Sub5_fitKernels','kernelsGLMv144.mat'));
v144kernels = v144kernels.GLMv144;


%% Run through and extract relevant data pairs

matchData = struct;
for i144 = 1:numel(v144kernels)
   for jv4 = 1:numel(v4Behav)
       if strcmpi(v144kernels(i144).name, v4Behav(jv4).name) && ...
               v144kernels(i144).ilobv == v4Behav(jv4).ilobv
            matchData(i144).name = v4Behav(jv4).name;
            matchData(i144).fov = v4Behav(jv4).fov;
            if v4Behav(jv4).ilobv
                matchData(i144).ilobv = true;
                matchData(i144).isim2 = false; 
            else
                matchData(i144).ilobv = false;
                matchData(i144).isim2 = true; 
            end
            matchData(i144).Rmov14 = [];
            matchData(i144).v4MovRatio = [];
            matchData(i144).Rrew14 = [];
            matchData(i144).v4RewRatio = [];
            
            v4iMov = find(v4Behav(jv4).iPCmov);
            v4iRew = find(v4Behav(jv4).iPCrew);
            
            % Extract data from v1 GLM analysis
            GLMv1 = v144kernels(i144).GLMkernelsv1;
            if isfield(GLMv1,'kernelMov')
                v1iData = find(~isnan(cat(1,GLMv1.sigBool_shufMov)));
                iPCv4 = cat(1,GLMv1.iPC_v4);
                iPCv4 = iPCv4(v1iData);
                a = ismember(iPCv4,v4iMov);
                
                kernelMov = [v144kernels(i144).movKernel_v1(:,a) v144kernels(i144).movKernel_v4n(:,a)];
                Rmov = corrcoef(kernelMov); Rmov = Rmov(end/2+1:end,1:end/2);
                matchData(i144).Rmov14 = Rmov(1:size(Rmov,1)+1:end)';
                
                iPCv4 = iPCv4(a);
                v4iMovMatched = find(ismember(v4iMov,iPCv4));
                v4MovData = v4Param(jv4).RcorrSpkBCell(v4iMovMatched,:);
                matchData(i144).v4MovRatio = v4MovData(:,3)./v4MovData(:,1);
            end
            if isfield(GLMv1,'kernelRew')
                v1iData = find(~isnan(cat(1,GLMv1.sigBool_shufRew)));
                iPCv4 = cat(1,GLMv1.iPC_v4);
                iPCv4 = iPCv4(v1iData);
                a = ismember(iPCv4,v4iRew);
                
                kernelRew = [v144kernels(i144).rewKernel_v1(:,a) v144kernels(i144).rewKernel_v4n(:,a)];
                Rrew = corrcoef(kernelRew); Rrew = Rrew(end/2+1:end,1:end/2);
                matchData(i144).Rrew14 = Rrew(1:size(Rrew,1)+1:end)';
                
                iPCv4 = iPCv4(a);
                v4iRewMatched = find(ismember(v4iRew,iPCv4));
                v4RewData = v4Param(jv4).lastOutcomeSpkCell(v4iRewMatched,:);
                matchData(i144).v4RewRatio = v4RewData(:,3)./v4RewData(:,2);
            end
%             
%             v4iMatched = cat(1,v144kernels(i144).GLMmatchv4n_v1.iPC_v4);
%             matchData(i144).v4iMovMatched = find(ismember(v4iMov,v4iMatched));
%             matchData(i144).v1iMovMatched = find(ismember(v4iMatched,v4iMov));
%             matchData(i144).v4iRewMatched = find(ismember(v4iRew,v4iMatched));
%             matchData(i144).v1iRewMatched = find(ismember(v4iMatched,v4iRew));
%             
%             % v4 movement ratio
%             v4MovData = v4Param(jv4).RcorrSpkBCell(matchData(i144).v4iMovMatched,:);
%             matchData(i144).v4MovRatio = v4MovData(:,3)./v4MovData(:,1);
%             
%             % v1 movement kernel
%             matchData(i144).v1MovKernelCorr = v4MovData(:,3)./v4MovData(:,1);
%             
%             abcd = 4;
       end
   end
end


%% Figure
lobvMov = [cat(1,matchData(cat(1,matchData.ilobv)).Rmov14) cat(1,matchData(cat(1,matchData.ilobv)).v4MovRatio)]; lobvMov = lobvMov(~isnan(lobvMov(:,1)),:);
oLobvMov = regress_hd(lobvMov(:,1),lobvMov(:,2),1);
lobvRew = [cat(1,matchData(cat(1,matchData.ilobv)).Rrew14) cat(1,matchData(cat(1,matchData.ilobv)).v4RewRatio)]; lobvRew = lobvRew(~isnan(lobvRew(:,1)),:);
oLobvRew = regress_hd(lobvRew(:,1),lobvRew(:,2),1);
sim2Mov = [cat(1,matchData(cat(1,matchData.isim2)).Rmov14) cat(1,matchData(cat(1,matchData.isim2)).v4MovRatio)]; sim2Mov = sim2Mov(~isnan(sim2Mov(:,1)),:);
oSim2Mov = regress_hd(sim2Mov(:,1),sim2Mov(:,2),1);
sim2Rew = [cat(1,matchData(cat(1,matchData.isim2)).Rrew14) cat(1,matchData(cat(1,matchData.isim2)).v4RewRatio)]; sim2Rew = sim2Rew(~isnan(sim2Rew(:,1)),:);
oSim2Rew = regress_hd(sim2Rew(:,1),sim2Rew(:,2),1);

figure;
subplot(2,2,1);
boundedline(oLobvMov.x_fit,oLobvMov.y_fit,oLobvMov.delta,'cmap',[0 0 0],'alpha'); hold on;
scatter(lobvMov(:,1),lobvMov(:,2));
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Lobule V: movement PCs');
xlabel('v1-v4 kernel correlation');
ylabel('Practiced movement ratio (Rhigh/Rlow)');
axis square;
axis([-1 1 0 3]);

subplot(2,2,2);
boundedline(oLobvRew.x_fit,oLobvRew.y_fit,oLobvRew.delta,'cmap',[0 0 0],'alpha'); hold on;
scatter(lobvRew(:,1),lobvRew(:,2));
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Lobule V: reward PCs');
xlabel('v1-v4 kernel correlation');
ylabel('Reward history ratio (O/C)');
axis square;
axis([-1 1 0 3]);

subplot(2,2,3);
boundedline(oSim2Mov.x_fit,oSim2Mov.y_fit,oSim2Mov.delta,'cmap',[0 0 0],'alpha'); hold on;
scatter(sim2Mov(:,1),sim2Mov(:,2));
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Simplex: movement PCs');
xlabel('v1-v4 kernel correlation');
ylabel('Practiced movement ratio');
axis square;
axis([-1 1 0 3]);

subplot(2,2,4);
boundedline(oSim2Rew.x_fit,oSim2Rew.y_fit,oSim2Rew.delta,'cmap',[0 0 0],'alpha'); hold on;
scatter(sim2Rew(:,1),sim2Rew(:,2));
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',14,'TickDir','out');
title('Simplex: reward PCs');
xlabel('v1-v4 kernel correlation');
ylabel('Reward history ratio (O/C)');
axis square;
axis([-1 1 0 3]);