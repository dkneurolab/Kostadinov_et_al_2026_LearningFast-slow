function fig3_sigSummary(sessions,dataFold,paperFold,glmParams,savebool)
%%
glmFold = sprintf('GLMfull_%s_%s_alpha%i',glmParams.dataType,glmParams.distr,glmParams.alpha*100);
figFold0 = fullfile(paperFold,'Figures','Fig3','Sub3_SummaryBars');
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
    
    glmSummary(iSesh).nPC = numel(glmLocal);
    glmSummary(iSesh).nSig = sum(cat(1,glmLocal.sigBool));
    glmSummary(iSesh).pPCs_sig = glmSummary(iSesh).nSig/glmSummary(iSesh).nPC*100;
    glmSummary(iSesh).nSigVis = sum(cat(1,glmLocal.sigBool_shufVis_proj));
    glmSummary(iSesh).pPCs_vis = glmSummary(iSesh).nSigVis/glmSummary(iSesh).nSig*100;
    glmSummary(iSesh).nSigMov = sum(cat(1,glmLocal.sigBool_shufMov_proj));
    glmSummary(iSesh).pPCs_mov = glmSummary(iSesh).nSigMov/glmSummary(iSesh).nSig*100;
    glmSummary(iSesh).nSigRew = sum(cat(1,glmLocal.sigBool_shufRew_proj));
    glmSummary(iSesh).pPCs_rew = glmSummary(iSesh).nSigRew/glmSummary(iSesh).nSig*100;
    glmSummary(iSesh).nSigLick = sum(cat(1,glmLocal.sigBool_shufLick_proj));
    glmSummary(iSesh).pPCs_lick = glmSummary(iSesh).nSigLick/glmSummary(iSesh).nSig*100;
    glmSummary(iSesh).devExp = cat(1,glmLocal.devExp);
    glmSummary(iSesh).devExp_full = cat(1,glmLocal.devExp_full);
    glmSummary(iSesh).devExp_shuf = cat(1,glmLocal.devExp_shuf);
    
    glmSummary(iSesh).pVisMov = sum(cat(1,glmLocal.sigBool_shufVis_proj) & cat(1,glmLocal.sigBool_shufMov_proj))/glmSummary(iSesh).nSig*100;
    glmSummary(iSesh).pVisMov_chance = sum(cat(1,glmLocal.sigBool_shufVis_proj))/glmSummary(iSesh).nSig * sum(cat(1,glmLocal.sigBool_shufMov_proj))/glmSummary(iSesh).nSig*100;
    
    glmSummary(iSesh).pVisRew = sum(cat(1,glmLocal.sigBool_shufVis_proj) & cat(1,glmLocal.sigBool_shufRew_proj))/glmSummary(iSesh).nSig*100;
    glmSummary(iSesh).pVisRew_chance = sum(cat(1,glmLocal.sigBool_shufVis_proj))/glmSummary(iSesh).nSig * sum(cat(1,glmLocal.sigBool_shufRew_proj))/glmSummary(iSesh).nSig*100;
    
    glmSummary(iSesh).pVisLick = sum(cat(1,glmLocal.sigBool_shufVis_proj) & cat(1,glmLocal.sigBool_shufLick_proj))/glmSummary(iSesh).nSig*100;
    glmSummary(iSesh).pVisLick_chance = sum(cat(1,glmLocal.sigBool_shufVis_proj))/glmSummary(iSesh).nSig * sum(cat(1,glmLocal.sigBool_shufLick_proj))/glmSummary(iSesh).nSig*100;
    
    glmSummary(iSesh).pMovRew = sum(cat(1,glmLocal.sigBool_shufMov_proj) & cat(1,glmLocal.sigBool_shufRew_proj))/glmSummary(iSesh).nSig*100;
    glmSummary(iSesh).pMovRew_chance = sum(cat(1,glmLocal.sigBool_shufMov_proj))/glmSummary(iSesh).nSig * sum(cat(1,glmLocal.sigBool_shufRew_proj))/glmSummary(iSesh).nSig*100;
    
    glmSummary(iSesh).pMovLick = sum(cat(1,glmLocal.sigBool_shufMov_proj) & cat(1,glmLocal.sigBool_shufLick_proj))/glmSummary(iSesh).nSig*100;
    glmSummary(iSesh).pMovLick_chance = sum(cat(1,glmLocal.sigBool_shufMov_proj))/glmSummary(iSesh).nSig * sum(cat(1,glmLocal.sigBool_shufLick_proj))/glmSummary(iSesh).nSig*100;
    
    glmSummary(iSesh).pRewLick = sum(cat(1,glmLocal.sigBool_shufRew_proj) & cat(1,glmLocal.sigBool_shufLick_proj))/glmSummary(iSesh).nSig*100;
    glmSummary(iSesh).pRewLick_chance = sum(cat(1,glmLocal.sigBool_shufRew_proj))/glmSummary(iSesh).nSig * sum(cat(1,glmLocal.sigBool_shufLick_proj))/glmSummary(iSesh).nSig*100;
        
    glmOld = load(fullfile(dataPath,'zz_FullDataOutputs_old','GLMsig.mat'));
    glmSummary(iSesh).nSigOld = sum(cat(1,glmOld.GLMoutput.sigBool));
    glmSummary(iSesh).pPCs_sigOld = glmSummary(iSesh).nSigOld/glmSummary(iSesh).nPC*100;
    
end
ilobv = logical(cat(1,glmSummary.iLobV));
isim2 = logical(cat(1,glmSummary.iSim2));

lobV = glmSummary(ilobv);
sim2 = glmSummary(isim2);

[~,bVis] = ttest(cat(1,lobV.pPCs_vis),cat(1,sim2.pPCs_vis));
[~,bMov] = ttest(cat(1,lobV.pPCs_mov),cat(1,sim2.pPCs_mov));
[~,bRew] = ttest(cat(1,lobV.pPCs_rew),cat(1,sim2.pPCs_rew));
[~,bLic] = ttest(cat(1,lobV.pPCs_lick),cat(1,sim2.pPCs_lick));

[aVis,~] = ranksum(cat(1,lobV.pPCs_vis),cat(1,sim2.pPCs_vis));
[aMov,~] = ranksum(cat(1,lobV.pPCs_mov),cat(1,sim2.pPCs_mov));
[aRew,~] = ranksum(cat(1,lobV.pPCs_rew),cat(1,sim2.pPCs_rew));
[aLic,~] = ranksum(cat(1,lobV.pPCs_lick),cat(1,sim2.pPCs_lick));
statSummary = v2struct(aVis, aMov, aRew, aLic);

if savebool
    save(fullfile(figFold,'sigSummary_v4n.mat'),'glmSummary','statSummary');
end

%% Joint coding probability
percVisMov = (100*cat(1,glmSummary.pVisMov)./cat(1,glmSummary.pVisMov_chance))-100;
percVisRew = (100*cat(1,glmSummary.pVisRew)./cat(1,glmSummary.pVisRew_chance))-100;
percVisLick = (100*cat(1,glmSummary.pVisLick)./cat(1,glmSummary.pVisLick_chance))-100;
percMovRew = (100*cat(1,glmSummary.pMovRew)./cat(1,glmSummary.pMovRew_chance))-100;
percMovLick = (100*cat(1,glmSummary.pMovLick)./cat(1,glmSummary.pMovLick_chance))-100;
percRewLick = (100*cat(1,glmSummary.pRewLick)./cat(1,glmSummary.pRewLick_chance))-100;

lobV_joint = zeros(4)*100; sim2_joint = lobV_joint;

lobV_joint(1,2) = mean(percVisMov(ilobv));
sim2_joint(1,2) = mean(percVisMov(isim2));
lobV_joint(1,3) = mean(percVisRew(ilobv));
sim2_joint(1,3) = mean(percVisRew(isim2));
lobV_joint(1,4) = mean(percVisLick(ilobv));
sim2_joint(1,4) = mean(percVisLick(isim2));
lobV_joint(2,3) = mean(percMovRew(ilobv));
sim2_joint(2,3) = mean(percMovRew(isim2));
lobV_joint(2,4) = mean(percMovLick(ilobv));
sim2_joint(2,4) = mean(percMovLick(isim2));
lobV_joint(3,4) = mean(percRewLick(ilobv));
sim2_joint(3,4) = mean(percRewLick(isim2));

% iUpper = logical(triu(ones(4),1));

lobV_joint = lobV_joint+lobV_joint';
sim2_joint = sim2_joint + sim2_joint';

% Now calculate who is 2 SDs away from mean
lobV_jointSD = zeros(4)*100; sim2_jointSD = lobV_jointSD;

lobV_jointSD(1,2) = std(percVisMov(ilobv));
sim2_jointSD(1,2) = std(percVisMov(isim2));
lobV_jointSD(1,3) = std(percVisRew(ilobv));
sim2_jointSD(1,3) = std(percVisRew(isim2));
lobV_jointSD(1,4) = std(percVisLick(ilobv));
sim2_jointSD(1,4) = std(percVisLick(isim2));
lobV_jointSD(2,3) = std(percMovRew(ilobv));
sim2_jointSD(2,3) = std(percMovRew(isim2));
lobV_jointSD(2,4) = std(percMovLick(ilobv));
sim2_jointSD(2,4) = std(percMovLick(isim2));
lobV_jointSD(3,4) = std(percRewLick(ilobv));
sim2_jointSD(3,4) = std(percRewLick(isim2));

lobV_jointSD = lobV_jointSD+lobV_jointSD';
sim2_jointSD = sim2_jointSD + sim2_jointSD';
lobV_jointSD(1:5:end) = 1;
sim2_jointSD(1:5:end) = 1;

lobV_jointSig = abs(lobV_joint)-2*lobV_jointSD;
sim2_jointSig = abs(sim2_joint)-2*sim2_jointSD;

lobV_jointSig_zscore = abs(lobV_joint)./lobV_jointSD;
sim2_jointSig_zscore = abs(sim2_joint)./sim2_jointSD;

lobV_jointSig_pVal = 2*(1-normcdf(lobV_jointSig_zscore));
sim2_jointSig_pVal = 2*(1-normcdf(sim2_jointSig_zscore));

lobV_jointSig(lobV_jointSig(:) <= 0) = 0; lobV_jointSig = logical(lobV_jointSig);
sim2_jointSig(sim2_jointSig(:) <= 0) = 0; sim2_jointSig = logical(sim2_jointSig);

fJoint = figure;
subplot(1,2,1)
imagesc(lobV_joint,[-50 50]);
% imagesc(lobV_joint./lobV_jointSD,[-2.5 2.5]);
axis square; colormap redblue;
title('LobV')
box off
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
xticks([.5:1:4.5]); yticks([.5:1:4.5]); xticklabels({}); yticklabels({});
colorbar
subplot(1,2,2)
imagesc(sim2_joint,[-50 50]);
% imagesc(sim2_joint./sim2_jointSD,[-2.5 2.5]);
axis square; colormap redblue;
title('sim2')
box off
xticks([.5:1:4.5]); yticks([.5:1:4.5]); xticklabels({}); yticklabels({});
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
colorbar

if savebool
    savename = sprintf('jointProbs');
    savefig(fJoint,fullfile(figFold,[savename,'.fig']));
    saveas(fJoint,fullfile(figFold,[savename,'.png']));
    print(fJoint,fullfile(figFold,[savename,'.eps']), '-depsc', '-vector');
end


%% Dev exp comparison
devExp = cat(1,glmSummary.devExp);
devExp_full = cat(1,glmSummary.devExp_full);
devExp_lobV = cat(1,lobV.devExp);
devExp_full_lobV = cat(1,lobV.devExp_full);
devExp_sim2 = cat(1,sim2.devExp);
devExp_full_sim2 = cat(1,sim2.devExp_full);
devMax = ceil(max([devExp; devExp_full]*10))/10;
devMin = floor(min([devExp; devExp_full]*10))/10;

xLnTotal = repmat([0.1;.9],1,numel(glmSummary));
pPCs_sig = cat(1,glmSummary.pPCs_sig);
pPCs_sigOld = cat(1,glmSummary.pPCs_sigOld);

fDevExp = figure;
subplot(1,3,1);
scatter(devExp,devExp_full,'MarkerFaceColor','b','MarkerEdgeColor','none',...
    'MarkerFaceAlpha',.1); axis square; hold on;
axis([devMin devMax devMin devMax]);
hRef = refline(1,0); hRef.Color = 'k';
box off
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
xlabel('Deviance explained: reduced model');
ylabel('Deviance explained: full model');

subplot(1,3,2);
scatter(devExp_lobV,devExp_full_lobV,'MarkerFaceColor','b','MarkerEdgeColor','none',...
    'MarkerFaceAlpha',.1); axis square; hold on;
scatter(devExp_sim2,devExp_full_sim2,'MarkerFaceColor','r','MarkerEdgeColor','none',...
    'MarkerFaceAlpha',.1); axis square;
axis([devMin devMax devMin devMax]);
hRef = refline(1,0); hRef.Color = 'k';
legend('Lobule V','Simplex');
box off
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
xlabel('Deviance explained: reduced model');
ylabel('Deviance explained: full model');

subplot(1,3,3);
hold on;
bar(1,mean(pPCs_sigOld),.3,'EdgeColor','none','FaceColor',[.5 .5 .5]);
bar(2,mean(pPCs_sig),.3,'EdgeColor','none','FaceColor',[.5 .5 .5]);
plot(1+xLnTotal,[pPCs_sigOld'; pPCs_sig'],'-k','LineWidth',2);
e1 = errorbar(1,mean(pPCs_sigOld),std(pPCs_sigOld)/sqrt(numel(pPCs_sigOld)));
e1.CapSize = 0; e1.LineWidth = 3.5; e1.Marker = 'none'; e1.Color = 'k';
e1 = errorbar(2,mean(pPCs_sig),std(pPCs_sig)/sqrt(numel(pPCs_sig)));
e1.CapSize = 0; e1.LineWidth = 3.5; e1.Marker = 'none'; e1.Color = 'k';
axis([0 3 0 100]); box off; axis square;
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
xticks([1 2]); xticklabels({'80% PCA var.','Best model'});
ylabel('Percent of PCs with significant fit');


if savebool
    savename = sprintf('devExpScatters');
    savefig(fDevExp,fullfile(figFold,[savename,'.fig']));
    saveas(fDevExp,fullfile(figFold,[savename,'.png']));
    print(fDevExp,fullfile(figFold,[savename,'.eps']), '-depsc', '-vector');
end

%%  Percent with sig fits

flobV = figure;
hold on;
axis([.5 5.5 0 100]);
yticks(0:20:100);
ylabel({'Percent of lobV PC dendrites','with significant fit'});
xticks(1:5);
xticklabels({'Full model','Trial ON','Movement','Reward','Licking'});
bar(1,cat(1,lobV.pPCs_sig),'EdgeColor','none','FaceColor',[.75 .75 .75]);
bar(2,cat(1,lobV.pPCs_vis),'EdgeColor','none','FaceColor',[.5 0 0]);
bar(3,cat(1,lobV.pPCs_mov),'EdgeColor','none','FaceColor',[.0 .3 .9]);
bar(4,cat(1,lobV.pPCs_rew),'EdgeColor','none','FaceColor',[0 .5 .125]);
bar(5,cat(1,lobV.pPCs_lick),'EdgeColor','none','FaceColor',[0.3010 0.7450 0.9330]);
box off
set(gca,'LineWidth',3,'XColor','k','YColor','k','Layer', 'Top','FontSize',40,'TickDir','out');

if savebool
    savename = sprintf('LobV_summaryBar');
    savefig(flobV,fullfile(figFold,[savename,'.fig']));
    saveas(flobV,fullfile(figFold,[savename,'.png']));
    print(flobV,fullfile(figFold,[savename,'.eps']), '-depsc', '-vector');
end

fsim2 = figure;
hold on;
axis([.5 5.5 0 100]);
yticks(0:20:100);
ylabel({'Percent of simplex PC dendrites','with significant fit'});
xticks(1:5);
xticklabels({'Full model','Trial ON','Movement','Reward','Licking'});
bar(1,cat(1,sim2.pPCs_sig),'EdgeColor','none','FaceColor',[.75 .75 .75]);
bar(2,cat(1,sim2.pPCs_vis),'EdgeColor','none','FaceColor',[.5 0 0]);
bar(3,cat(1,sim2.pPCs_mov),'EdgeColor','none','FaceColor',[.0 .3 .9]);
bar(4,cat(1,sim2.pPCs_rew),'EdgeColor','none','FaceColor',[0 .5 .125]);
bar(5,cat(1,sim2.pPCs_lick),'EdgeColor','none','FaceColor',[0.3010 0.7450 0.9330]);
box off
set(gca,'LineWidth',3,'XColor','k','YColor','k','Layer', 'Top','FontSize',40,'TickDir','out');

if savebool
    savename = sprintf('sim2_summaryBar');
    savefig(fsim2,fullfile(figFold,[savename,'.fig']));
    saveas(fsim2,fullfile(figFold,[savename,'.png']));
    print(fsim2,fullfile(figFold,[savename,'.eps']), '-depsc', '-vector');
end

%% Combined figure with single plots
xLn = repmat([0.1;.9],1,numel(lobV));
% xPts = linspace(-1/4,1/4,numel(lobV));

fboth = figure;
hold on;

% Visual
lobVplot = cat(1,lobV.pPCs_vis);
sim2plot = cat(1,sim2.pPCs_vis);
bar(1,mean(lobVplot),'EdgeColor','none','FaceColor',[.75 .0 .5]);
bar(2,mean(sim2plot),'EdgeColor','none','FaceColor',[0 .5 .125]);
plot(1+xLn,[lobVplot'; sim2plot'],'-k','LineWidth',2);
e1 = errorbar(1,mean(lobVplot),std(lobVplot)/sqrt(numel(lobVplot)));
e1.CapSize = 0; e1.LineWidth = 3.5; e1.Marker = 'none'; e1.Color = 'k';
e1 = errorbar(2,mean(sim2plot),std(sim2plot)/sqrt(numel(sim2plot)));
e1.CapSize = 0; e1.LineWidth = 3.5; e1.Marker = 'none'; e1.Color = 'k';
text(1.5,100,sprintf('p = %.2f',aVis),'FontSize',12,'HorizontalAlignment','center')
% scatter(1+xPts,lobVplot,'filled','MarkerFaceColor',[.75 .0 .5]/2);
% scatter(2+xPts,sim2plot,'filled','MarkerFaceColor',[0 .5 .125]/2);

% Movement
lobVplot = cat(1,lobV.pPCs_mov);
sim2plot = cat(1,sim2.pPCs_mov);
bar(4,mean(cat(1,lobV.pPCs_mov)),'EdgeColor','none','FaceColor',[.75 .0 .5]);
bar(5,mean(cat(1,sim2.pPCs_mov)),'EdgeColor','none','FaceColor',[0 .5 .125]);
plot(4+xLn,[lobVplot'; sim2plot'],'-k','LineWidth',2);
e1 = errorbar(4,mean(cat(1,lobV.pPCs_mov)),std(cat(1,lobV.pPCs_mov))/sqrt(numel(lobV)));
e1.CapSize = 0; e1.LineWidth = 3.5; e1.Marker = 'none'; e1.Color = 'k';
e1 = errorbar(5,mean(cat(1,sim2.pPCs_mov)),std(cat(1,sim2.pPCs_mov))/sqrt(numel(sim2)));
e1.CapSize = 0; e1.LineWidth = 3.5; e1.Marker = 'none'; e1.Color = 'k';
text(4.5,100,sprintf('p = %.2f',aMov),'FontSize',12,'HorizontalAlignment','center')
% scatter(4+xPts,cat(1,lobV.pPCs_mov),'filled','MarkerFaceColor',[0 .5 .125]/2);
% scatter(5+xPts,cat(1,sim2.pPCs_mov),'filled','MarkerFaceColor',[0 .5 .125]/2);

% Reward
lobVplot = cat(1,lobV.pPCs_rew);
sim2plot = cat(1,sim2.pPCs_rew);
bar(7,mean(cat(1,lobV.pPCs_rew)),'EdgeColor','none','FaceColor',[.75 .0 .5]);
bar(8,mean(cat(1,sim2.pPCs_rew)),'EdgeColor','none','FaceColor',[0 .5 .125]);
plot(7+xLn,[lobVplot'; sim2plot'],'-k','LineWidth',2);
e1 = errorbar(7,mean(cat(1,lobV.pPCs_rew)),std(cat(1,lobV.pPCs_rew))/sqrt(numel(lobV)));
e1.CapSize = 0; e1.LineWidth = 3.5; e1.Marker = 'none'; e1.Color = 'k';
e1 = errorbar(8,mean(cat(1,sim2.pPCs_rew)),std(cat(1,sim2.pPCs_rew))/sqrt(numel(sim2)));
e1.CapSize = 0; e1.LineWidth = 3.5; e1.Marker = 'none'; e1.Color = 'k';
text(7.5,100,sprintf('p = %.2f',aRew),'FontSize',12,'HorizontalAlignment','center')
% scatter(7+xPts,cat(1,lobV.pPCs_rew),'filled','MarkerFaceColor',[0 .5 .125]/2);
% scatter(8+xPts,cat(1,sim2.pPCs_rew),'filled','MarkerFaceColor',[0 .5 .125]/2);

% Lick
lobVplot = cat(1,lobV.pPCs_lick);
sim2plot = cat(1,sim2.pPCs_lick);
bar(10,mean(cat(1,lobV.pPCs_lick)),'EdgeColor','none','FaceColor',[.75 .0 .5]);
bar(11,mean(cat(1,sim2.pPCs_lick)),'EdgeColor','none','FaceColor',[0 .5 .125]);
plot(10+xLn,[lobVplot'; sim2plot'],'-k','LineWidth',2);
e1 = errorbar(10,mean(cat(1,lobV.pPCs_lick)),std(cat(1,lobV.pPCs_lick))/sqrt(numel(lobV)));
e1.CapSize = 0; e1.LineWidth = 3.5; e1.Marker = 'none'; e1.Color = 'k';
e1 = errorbar(11,mean(cat(1,sim2.pPCs_lick)),std(cat(1,sim2.pPCs_lick))/sqrt(numel(sim2)));
e1.CapSize = 0; e1.LineWidth = 3.5; e1.Marker = 'none'; e1.Color = 'k';
text(10.5,100,sprintf('p = %.2f',aLic),'FontSize',12,'HorizontalAlignment','center')
% scatter(10+xPts,cat(1,lobV.pPCs_lick),'filled','MarkerFaceColor',[0 .5 .125]/2);
% scatter(11+xPts,cat(1,sim2.pPCs_lick),'filled','MarkerFaceColor',[0 .5 .125]/2);

axis([.0 12 0 100]);
yticks(0:20:100);
ylabel({'Percent of PC dendrites','with significant fit'});
xticks(1.5:3:10.5);
xticklabels({'Trial ON','Movement','Reward','Licking'});
box off
set(gca,'LineWidth',3,'XColor','k','YColor','k','Layer', 'Top','FontSize',40,'TickDir','out');

if savebool
    savename = sprintf('LobVvsSim2_summaryBar');
    savefig(fboth,fullfile(figFold,[savename,'.fig']));
    saveas(fboth,fullfile(figFold,[savename,'.png']));
    print(fboth,fullfile(figFold,[savename,'.eps']), '-depsc', '-vector');
end

end
