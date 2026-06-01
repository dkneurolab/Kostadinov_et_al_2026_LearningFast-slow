function fig3_matchGLMbars_v44(sessionsv4_1,sessionsv4_n,glmParams,dataFold,paperFold,savebool)

%% Test v4_n GLM on v4_1 data
aa = 0;
GLMout = struct;
for iSesh = 1:numel(sessionsv4_n)    
    for jSesh = 1:numel(sessionsv4_1)
        if strcmpi(sessionsv4_n(iSesh).name,sessionsv4_1(jSesh).name) && ...
                strcmpi(sessionsv4_n(iSesh).fov,sessionsv4_1(jSesh).fov)
            aa = aa + 1;
            GLMoutTemp = fig3_matchGLM_v44(sessionsv4_n(iSesh),sessionsv4_1(jSesh),glmParams,dataFold,paperFold,savebool);            
            GLMout(aa).name = sessionsv4_1(jSesh).name;
            if strcmpi(sessionsv4_1(jSesh).fov,'lobv')
                GLMout(aa).ilobv = true;
                GLMout(aa).isim2 = false;
            end
            if strcmpi(sessionsv4_1(jSesh).fov,'sim2')
                GLMout(aa).ilobv = false;
                GLMout(aa).isim2 = true;
            end
            GLMout(aa).GLMstruct = GLMoutTemp;
%             GLMout{aa,1} = sessionsv1(jSesh).name;
%             GLMout{aa,2} = sessionsv1(jSesh).fov;
%             GLMout{aa,3} = GLMoutTemp;
        end
    end
end

% Save maybe
if savebool
    figFold = fullfile(paperFold,'Figures','Fig3');
    if ~exist(figFold,'dir'); mkdir(figFold); end
    save(fullfile(figFold,'GLM_v4_v4.mat'),'GLMout');
end

%% GLM_v1 summary
GLMout_lobv = GLMout(logical(cat(1,GLMout.ilobv)));
GLMout_sim2 = GLMout(logical(cat(1,GLMout.isim2)));

for ilobv = 1:numel(GLMout_lobv)
    glmStr = GLMout_lobv(ilobv).GLMstruct;
    glmStr = glmStr(logical(cat(1,glmStr.sigV4)));
    GLMout_lobv(ilobv).pSigV1 = sum(cat(1,glmStr.sigBool))/numel(glmStr)*100;
    boolVar = cat(1,glmStr.sigBool_shufVis);
    GLMout_lobv(ilobv).pSigV1vis = sum(boolVar,'omitnan')/sum(~isnan(boolVar))*100;
    boolVar = cat(1,glmStr.sigBool_shufMov);
    GLMout_lobv(ilobv).pSigV1mov = sum(boolVar,'omitnan')/sum(~isnan(boolVar))*100;
    boolVar = cat(1,glmStr.sigBool_shufRew);
    GLMout_lobv(ilobv).pSigV1rew = sum(boolVar,'omitnan')/sum(~isnan(boolVar))*100;
    boolVar = cat(1,glmStr.sigBool_shufLick);
    GLMout_lobv(ilobv).pSigV1lick = sum(boolVar,'omitnan')/sum(~isnan(boolVar))*100;    
end

for isim2 = 1:numel(GLMout_sim2)
    glmStr = GLMout_sim2(isim2).GLMstruct;
    glmStr = glmStr(logical(cat(1,glmStr.sigV4)));
    GLMout_sim2(isim2).pSigV1 = sum(cat(1,glmStr.sigBool))/numel(glmStr)*100;
    boolVar = cat(1,glmStr.sigBool_shufVis);
    GLMout_sim2(isim2).pSigV1vis = sum(boolVar,'omitnan')/sum(~isnan(boolVar))*100;
    boolVar = cat(1,glmStr.sigBool_shufMov);
    GLMout_sim2(isim2).pSigV1mov = sum(boolVar,'omitnan')/sum(~isnan(boolVar))*100;
    boolVar = cat(1,glmStr.sigBool_shufRew);
    GLMout_sim2(isim2).pSigV1rew = sum(boolVar,'omitnan')/sum(~isnan(boolVar))*100;
    boolVar = cat(1,glmStr.sigBool_shufLick);
    GLMout_sim2(isim2).pSigV1lick = sum(boolVar,'omitnan')/sum(~isnan(boolVar))*100;    
end

%% Significance
[~,bAll] = ttest(cat(1,GLMout_lobv.pSigV1),cat(1,GLMout_sim2.pSigV1));
[~,bVis] = ttest(cat(1,GLMout_lobv.pSigV1vis),cat(1,GLMout_sim2.pSigV1vis));
[~,bMov] = ttest(cat(1,GLMout_lobv.pSigV1mov),cat(1,GLMout_sim2.pSigV1mov));
[~,bRew] = ttest(cat(1,GLMout_lobv.pSigV1rew),cat(1,GLMout_sim2.pSigV1rew));
[~,bLic] = ttest(cat(1,GLMout_lobv.pSigV1lick),cat(1,GLMout_sim2.pSigV1lick));


%% Make figure
figBarV1 = figure;
subplot(1,2,1);
hold on;
axis([.5 5.5 0 80]);
yticks(0:10:80);
ylabel({'Percent of Lobule V PC dendrites','with significant fit'});
xticks(1:5);
xticklabels({'Full model','Trial ON','Movement','Reward','Licking'});
bar(1,mean(cat(1,GLMout_lobv.pSigV1)),'EdgeColor','none','FaceColor',[.75 .75 .75]);
e1 = errorbar(1,mean(cat(1,GLMout_lobv.pSigV1)),std(cat(1,GLMout_lobv.pSigV1))/sqrt(numel(GLMout_lobv)));
e1.CapSize = 0; e1.LineWidth = 3.5; e1.Marker = 'none'; e1.Color = 'k';
bar(2,mean(cat(1,GLMout_lobv.pSigV1vis)),'EdgeColor','none','FaceColor',[.5 0 0]);
e1 = errorbar(2,mean(cat(1,GLMout_lobv.pSigV1vis)),std(cat(1,GLMout_lobv.pSigV1vis))/sqrt(numel(GLMout_lobv)));
e1.CapSize = 0; e1.LineWidth = 3.5; e1.Marker = 'none'; e1.Color = 'k';
bar(3,mean(cat(1,GLMout_lobv.pSigV1mov)),'EdgeColor','none','FaceColor',[.0 .3 .9]);
e1 = errorbar(3,mean(cat(1,GLMout_lobv.pSigV1mov)),std(cat(1,GLMout_lobv.pSigV1mov))/sqrt(numel(GLMout_lobv)));
e1.CapSize = 0; e1.LineWidth = 3.5; e1.Marker = 'none'; e1.Color = 'k';
bar(4,mean(cat(1,GLMout_lobv.pSigV1rew)),'EdgeColor','none','FaceColor',[0 .5 .125]);
e1 = errorbar(4,mean(cat(1,GLMout_lobv.pSigV1rew)),std(cat(1,GLMout_lobv.pSigV1rew))/sqrt(numel(GLMout_lobv)));
e1.CapSize = 0; e1.LineWidth = 3.5; e1.Marker = 'none'; e1.Color = 'k';
bar(5,mean(cat(1,GLMout_lobv.pSigV1lick)),'EdgeColor','none','FaceColor',[0.3010 0.7450 0.9330]);
e1 = errorbar(5,mean(cat(1,GLMout_lobv.pSigV1lick)),std(cat(1,GLMout_lobv.pSigV1lick))/sqrt(numel(GLMout_lobv)));
e1.CapSize = 0; e1.LineWidth = 3.5; e1.Marker = 'none'; e1.Color = 'k';
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
axis square;

subplot(1,2,2);
hold on;
axis([.5 5.5 0 80]);
yticks(0:10:80);
ylabel({'Percent of Lobule simplex PC dendrites','with significant fit'});
xticks(1:5);
xticklabels({'Full model','Trial ON','Movement','Reward','Licking'});
bar(1,mean(cat(1,GLMout_sim2.pSigV1)),'EdgeColor','none','FaceColor',[.75 .75 .75]);
e1 = errorbar(1,mean(cat(1,GLMout_sim2.pSigV1)),std(cat(1,GLMout_sim2.pSigV1))/sqrt(numel(GLMout_sim2)));
e1.CapSize = 0; e1.LineWidth = 3.5; e1.Marker = 'none'; e1.Color = 'k';
bar(2,mean(cat(1,GLMout_sim2.pSigV1vis)),'EdgeColor','none','FaceColor',[.5 0 0]);
e1 = errorbar(2,mean(cat(1,GLMout_sim2.pSigV1vis)),std(cat(1,GLMout_sim2.pSigV1vis))/sqrt(numel(GLMout_sim2)));
e1.CapSize = 0; e1.LineWidth = 3.5; e1.Marker = 'none'; e1.Color = 'k';
bar(3,mean(cat(1,GLMout_sim2.pSigV1mov)),'EdgeColor','none','FaceColor',[.0 .3 .9]);
e1 = errorbar(3,mean(cat(1,GLMout_sim2.pSigV1mov)),std(cat(1,GLMout_sim2.pSigV1mov))/sqrt(numel(GLMout_sim2)));
e1.CapSize = 0; e1.LineWidth = 3.5; e1.Marker = 'none'; e1.Color = 'k';
bar(4,mean(cat(1,GLMout_sim2.pSigV1rew)),'EdgeColor','none','FaceColor',[0 .5 .125]);
e1 = errorbar(4,mean(cat(1,GLMout_sim2.pSigV1rew)),std(cat(1,GLMout_sim2.pSigV1rew))/sqrt(numel(GLMout_sim2)));
e1.CapSize = 0; e1.LineWidth = 3.5; e1.Marker = 'none'; e1.Color = 'k';
bar(5,mean(cat(1,GLMout_sim2.pSigV1lick)),'EdgeColor','none','FaceColor',[0.3010 0.7450 0.9330]);
e1 = errorbar(5,mean(cat(1,GLMout_sim2.pSigV1lick)),std(cat(1,GLMout_sim2.pSigV1lick))/sqrt(numel(GLMout_sim2)));
e1.CapSize = 0; e1.LineWidth = 3.5; e1.Marker = 'none'; e1.Color = 'k';
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
axis square;

if savebool
    savename = sprintf('LobVvsSim2_summaryBar_v4_v4');
    savefig(figBarV1,fullfile(figFold,[savename,'.fig']));
    saveas(figBarV1,fullfile(figFold,[savename,'.png']));
    print(figBarV1,fullfile(figFold,[savename,'.eps']), '-depsc', '-painters');    
end