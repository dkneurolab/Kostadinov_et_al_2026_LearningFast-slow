%% fig2_GLMsummaryStats

% Paths to lobV FOVs:
lobV(1).path =  '/Users/dimitar/Desktop/CF_learning_paper/01_Behav+imaging/Version_4/DK103/2018-02-22_lobv/GLM/glmRfiles';
lobV(2).path = '/Users/dimitar/Desktop/CF_learning_paper/01_Behav+imaging/Version_4/DK105/2018-03-03_lobv/GLM/glmRfiles';
lobV(3).path = '/Users/dimitar/Desktop/CF_learning_paper/01_Behav+imaging/Version_4/DK169/2019-12-06_lobv/GLM/glmRfiles';
lobV(4).path = '/Users/dimitar/Desktop/CF_learning_paper/01_Behav+imaging/Version_4/DK171/2019-12-02_lobv/GLM/glmRfiles';
lobV(5).path = '/Users/dimitar/Desktop/CF_learning_paper/01_Behav+imaging/Version_4/DK194/2020-10-02_lobv/GLM/glmRfiles';
lobV(6).path = '/Users/dimitar/Desktop/CF_learning_paper/01_Behav+imaging/Version_4/DK199/2020-09-07_lobv/GLM/glmRfiles';


% Paths to simplex FOVs:
sim2(1).path = '/Users/dimitar/Desktop/CF_learning_paper/01_Behav+imaging/Version_4/DK103/2018-03-19_sim2/GLM/glmRfiles';
sim2(2).path = '/Users/dimitar/Desktop/CF_learning_paper/01_Behav+imaging/Version_4/DK105/2018-02-23_sim2/GLM/glmRfiles';
sim2(3).path = '/Users/dimitar/Desktop/CF_learning_paper/01_Behav+imaging/Version_4/DK169/2019-12-02_sim2/GLM/glmRfiles';
sim2(4).path = '/Users/dimitar/Desktop/CF_learning_paper/01_Behav+imaging/Version_4/DK171/2019-12-12_sim2/GLM/glmRfiles';
sim2(5).path = '/Users/dimitar/Desktop/CF_learning_paper/01_Behav+imaging/Version_4/DK194/2020-09-15_sim2/GLM/glmRfiles';
sim2(6).path = '/Users/dimitar/Desktop/CF_learning_paper/01_Behav+imaging/Version_4/DK199/2020-09-22_sim2/GLM/glmRfiles';

saveDir = '/Users/dimitar/Dropbox/Workstuff/Presentations/2022-03 lab meeting';

%% Extract data for lobV bar graph

for ilobv = 1:numel(lobV)
    GLMout_path = fullfile(lobV(ilobv).path,'GLMsig+VisMovRewLick.mat');
    
    if exist(GLMout_path,'file')
        load(GLMout_path);
        lobV(ilobv).nPCs = numel(GLMoutput);
        lobV(ilobv).nPCs_sig = sum(cat(1,GLMoutput.sigBool));
        lobV(ilobv).nPCs_vis = sum(cat(1,GLMoutput.sigBool_shufVis));
        lobV(ilobv).nPCs_mov = sum(cat(1,GLMoutput.sigBool_shufMov));
        lobV(ilobv).nPCs_rew = sum(cat(1,GLMoutput.sigBool_shufRew));
        lobV(ilobv).nPCs_lick = sum(cat(1,GLMoutput.sigBool_shufLick));
        lobV(ilobv).pPCs_sig = 100*sum(cat(1,GLMoutput.sigBool))/lobV(ilobv).nPCs;
        lobV(ilobv).pPCs_vis = 100*sum(cat(1,GLMoutput.sigBool_shufVis))/lobV(ilobv).nPCs_sig;
        lobV(ilobv).pPCs_mov = 100*sum(cat(1,GLMoutput.sigBool_shufMov))/lobV(ilobv).nPCs_sig;
        lobV(ilobv).pPCs_rew = 100*sum(cat(1,GLMoutput.sigBool_shufRew))/lobV(ilobv).nPCs_sig;
        lobV(ilobv).pPCs_lick = 100*sum(cat(1,GLMoutput.sigBool_shufLick))/lobV(ilobv).nPCs_sig;
        
    else
        lobV(ilobv).nPCs = nan;
        lobV(ilobv).nPCs_sig = nan;
        lobV(ilobv).nPCs_vis = nan;
        lobV(ilobv).nPCs_mov = nan;
        lobV(ilobv).nPCs_rew = nan;
        lobV(ilobv).nPCs_lick = nan;
        lobV(ilobv).pPCs_sig = nan;
        lobV(ilobv).pPCs_vis = nan;
        lobV(ilobv).pPCs_mov = nan;
        lobV(ilobv).pPCs_rew = nan;
        lobV(ilobv).pPCs_lick = nan;
        
    end
end

%% Extract data for simplex bar graph

for isim2 = 1:numel(sim2)
    GLMout_path = fullfile(sim2(isim2).path,'GLMsig+VisMovRewLick.mat');
    
    if exist(GLMout_path,'file')
        load(GLMout_path);
        sim2(isim2).nPCs = numel(GLMoutput);
        sim2(isim2).nPCs_sig = sum(cat(1,GLMoutput.sigBool));
        sim2(isim2).nPCs_vis = sum(cat(1,GLMoutput.sigBool_shufVis));
        sim2(isim2).nPCs_mov = sum(cat(1,GLMoutput.sigBool_shufMov));
        sim2(isim2).nPCs_rew = sum(cat(1,GLMoutput.sigBool_shufRew));
        sim2(isim2).nPCs_lick = sum(cat(1,GLMoutput.sigBool_shufLick));
        sim2(isim2).pPCs_sig = 100*sum(cat(1,GLMoutput.sigBool))/sim2(isim2).nPCs;
        sim2(isim2).pPCs_vis = 100*sum(cat(1,GLMoutput.sigBool_shufVis))/sim2(isim2).nPCs_sig;
        sim2(isim2).pPCs_mov = 100*sum(cat(1,GLMoutput.sigBool_shufMov))/sim2(isim2).nPCs_sig;
        sim2(isim2).pPCs_rew = 100*sum(cat(1,GLMoutput.sigBool_shufRew))/sim2(isim2).nPCs_sig;
        sim2(isim2).pPCs_lick = 100*sum(cat(1,GLMoutput.sigBool_shufLick))/sim2(isim2).nPCs_sig;
    else
        sim2(isim2).nPCs = nan;
        sim2(isim2).nPCs_sig = nan;
        sim2(isim2).nPCs_vis = nan;
        sim2(isim2).nPCs_mov = nan;
        sim2(isim2).nPCs_rew = nan;
        sim2(isim2).nPCs_lick = nan;
        sim2(isim2).pPCs_sig = nan;
        sim2(isim2).pPCs_vis = nan;
        sim2(isim2).pPCs_mov = nan;
        sim2(isim2).pPCs_rew = nan;
        sim2(isim2).pPCs_lick = nan;
        
    end
end

%%
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

savename = sprintf('LobV_summaryBar');
savefig(flobV,fullfile(saveDir,[savename,'.fig']));
saveas(flobV,fullfile(saveDir,[savename,'.png']));
print(flobV,fullfile(saveDir,[savename,'.eps']), '-depsc', '-painters');


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

savename = sprintf('sim2_summaryBar');
savefig(fsim2,fullfile(saveDir,[savename,'.fig']));
saveas(fsim2,fullfile(saveDir,[savename,'.png']));
print(fsim2,fullfile(saveDir,[savename,'.eps']), '-depsc', '-painters');

%% Comparison

fboth = figure;
hold on;
bar(1,mean(cat(1,lobV.pPCs_vis)),'EdgeColor','none','FaceColor',[.75 .0 .5]);
e1 = errorbar(1,mean(cat(1,lobV.pPCs_vis)),std(cat(1,lobV.pPCs_vis))/sqrt(numel(lobV)));
e1.CapSize = 0; e1.LineWidth = 2.5; e1.Marker = 'none'; e1.Color = 'k';
bar(2,mean(cat(1,sim2.pPCs_vis)),'EdgeColor','none','FaceColor',[0 .5 .125]);
e1 = errorbar(2,mean(cat(1,sim2.pPCs_vis)),std(cat(1,sim2.pPCs_vis))/sqrt(numel(sim2)));
e1.CapSize = 0; e1.LineWidth = 2.5; e1.Marker = 'none'; e1.Color = 'k';

bar(4,mean(cat(1,lobV.pPCs_mov)),'EdgeColor','none','FaceColor',[.75 .0 .5]);
e1 = errorbar(4,mean(cat(1,lobV.pPCs_mov)),std(cat(1,lobV.pPCs_mov))/sqrt(numel(lobV)));
e1.CapSize = 0; e1.LineWidth = 2.5; e1.Marker = 'none'; e1.Color = 'k';
bar(5,mean(cat(1,sim2.pPCs_mov)),'EdgeColor','none','FaceColor',[0 .5 .125]);
e1 = errorbar(5,mean(cat(1,sim2.pPCs_mov)),std(cat(1,sim2.pPCs_mov))/sqrt(numel(sim2)));
e1.CapSize = 0; e1.LineWidth = 2.5; e1.Marker = 'none'; e1.Color = 'k';

bar(7,mean(cat(1,lobV.pPCs_rew)),'EdgeColor','none','FaceColor',[.75 .0 .5]);
e1 = errorbar(7,mean(cat(1,lobV.pPCs_rew)),std(cat(1,lobV.pPCs_rew))/sqrt(numel(lobV)));
e1.CapSize = 0; e1.LineWidth = 2.5; e1.Marker = 'none'; e1.Color = 'k';
bar(8,mean(cat(1,sim2.pPCs_rew)),'EdgeColor','none','FaceColor',[0 .5 .125]);
e1 = errorbar(8,mean(cat(1,sim2.pPCs_rew)),std(cat(1,sim2.pPCs_rew))/sqrt(numel(sim2)));
e1.CapSize = 0; e1.LineWidth = 2.5; e1.Marker = 'none'; e1.Color = 'k';

axis([.0 9 0 100]);
yticks(0:20:100);
ylabel({'Percent of PC dendrites','with significant fit'});
xticks(1.5:3:7.5);
xticklabels({'Trial ON','Movement','Reward'});
box off
set(gca,'LineWidth',3,'XColor','k','YColor','k','Layer', 'Top','FontSize',40,'TickDir','out');

savename = sprintf('LobVvsSim2_summaryBar');
savefig(fboth,fullfile(saveDir,[savename,'.fig']));
saveas(fboth,fullfile(saveDir,[savename,'.png']));
print(fboth,fullfile(saveDir,[savename,'.eps']), '-depsc', '-painters');

