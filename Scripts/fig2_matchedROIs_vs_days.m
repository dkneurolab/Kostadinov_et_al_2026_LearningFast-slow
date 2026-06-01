function matchData = fig1_matchedROIs_vs_days(dataFold,paperFold,matchColors,savebool)
% Calculate percentage of matching cells across days for all recordings


%% Set paths and get directories
% dataFold = '/Users/dimitar/Desktop/CF_learning_paper';
imFold = fullfile(dataFold,'02_Imaging');
imStructs = dir(fullfile(dataFold,'*good2.mat'));
% matchColors = distinguishable_colors(100,[1 1 1; 0 0 0]);
% matchColors(1:6,:) = matchColors([6,5,2,3,1,4],:);


%% Extract info from analysis files
a = 0;
for i = 1:numel(imStructs)
    mouseName = imStructs(i).name(1:5);
    mouseStruct = load(fullfile(dataFold,imStructs(i).name));
    mouseStruct = mouseStruct.(imStructs(i).name(1:10));
    alignFolds = dir(fullfile(imFold,mouseName,'REG*'));
    allDates = cat(1,mouseStruct.date)';
    for j = 1:numel(alignFolds)
        pairfold = alignFolds(j).name;
        datestr1 = [pairfold(5:6),'-',pairfold(7:8)];
        datestr2 = [pairfold(10:11),'-',pairfold(12:13)];
        iDate1 = ceil(strfind(allDates(:)',datestr1)/10);
        iDate2 = ceil(strfind(allDates(:)',datestr2)/10);
        aDays = abs(days(datetime(allDates(:,iDate1)')-datetime(allDates(:,iDate2)')));
        alignData = dir(fullfile(alignFolds(j).folder,alignFolds(j).name,'*_analysis.mat'));
        alignData = load(fullfile(alignData.folder,alignData.name));
        alignData = alignData.regi;
        
        a = a+1;
        matchData(a).name = mouseName;
        matchData(a).fov = mouseStruct(iDate1).fov;
        switch matchData(a).name
            case 'dk103'
                matchData(a).dateInj = '2017-12-16';
            case 'dk105'
                matchData(a).dateInj = '2017-12-16';
            case 'dk169'
                matchData(a).dateInj = '2019-10-17';
            case 'dk171'
                matchData(a).dateInj = '2019-10-17';
            case 'dk194'
                matchData(a).dateInj = '2020-06-29';
            case 'dk199'
                matchData(a).dateInj = '2020-06-29';
        end
        matchData(a).date1 = mouseStruct(iDate1).date;
        matchData(a).date2 = mouseStruct(iDate2).date;
        matchData(a).nRoi1 = numel(alignData.dat.files(1).iscell);
        matchData(a).nRoi2 = numel(alignData.dat.files(2).iscell);
%         matchData(a).nRoi1 = alignData.dat.files(1).num_rois;
%         matchData(a).nRoi2 = alignData.dat.files(2).num_rois;
        matchData(a).delDays = aDays;
        matchData(a).nMatched = size(alignData.rois.idcs,1);
        matchData(a).nTotal = min(cat(1,alignData.dat.files.num_rois));
        matchData(a).percMatched = matchData(a).nMatched/matchData(a).nTotal;
        if strcmpi(matchData(a).fov,'sim2')
            matchData(a).matchColor = matchColors(i,:);
        else
            matchData(a).matchColor = matchColors(i,:)*0.8;
        end
        
    end
    
    
end

%%
delDays = cat(1,matchData.delDays);
percMatch = cat(1,matchData.percMatched);
matchColor = cat(1,matchData.matchColor);
iGrps = kmeans(matchColor,size(unique(matchColor,'rows'),1));

%%
uGrps = unique(iGrps);

roisVtime = struct;
for j = 1:numel(uGrps)
    matchGrp = matchData(iGrps == j);
    dateInj = matchGrp(1).dateInj;
    dateExp = [cat(1,matchGrp.date1); cat(1,matchGrp.date2)];
    daysPostInj = datenum(dateExp)-datenum(dateInj);    
    nRoisLocal = [cat(1,matchGrp.nRoi1); cat(1,matchGrp.nRoi2)];
    nRoisNorm = nRoisLocal/max(nRoisLocal);
    
    % Take only unique FOVs and sort by days post injection
    [~,iFovs] = unique(daysPostInj);
    nRoisNorm = nRoisNorm(iFovs); daysPostInj = daysPostInj(iFovs);
    [~,iSort] = sort(daysPostInj);
    
    roisVtime(j).name = matchGrp(1).name;
    roisVtime(j).fov = matchGrp(1).fov;
    roisVtime(j).daysPostInj = daysPostInj(iSort);
    roisVtime(j).nRoisNorm = nRoisNorm(iSort);
    roisVtime(j).matchColor = matchGrp(1).matchColor;
end

%%
oMatch = regress_hd(delDays,percMatch,1);
oRoi = regress_hd(cat(1,roisVtime.daysPostInj),cat(1,roisVtime.nRoisNorm),1);

figout = figure;
subplot(2,2,1); hold on;
s1 = scatter(delDays,percMatch,108,matchColor,'filled');
s1.MarkerEdgeColor = 'none'; %s1.MarkerFaceAlpha = 0.5;
h1 = boundedline(oMatch.x_fit,oMatch.y_fit,oMatch.delta,'k','alpha');
set(h1,'linewidth',2);
axis square;
xlim([-1.5 75]); ylim([-.02 1]);
xticks([0:15:75]); yticks([0:0.2:1]);
ylabel('Proportion matched');
xlabel('Days between sessions');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
text(45,0.8,{sprintf('rsq = %0.3f',oMatch.rsq); sprintf('p = %0.3f',oMatch.corr_p); sprintf('50percDay = %i',round(oMatch.x_fit(find(oMatch.y_fit < .50,1))))});

subplot(2,2,2);
hold on;
for k = 1:numel(uGrps)
    p(k) = plot(delDays(iGrps==k),percMatch(iGrps==k),'o-');
    p(k).Color = [matchColor(find(iGrps==k,1),:) 0.3];
    p(k).MarkerFaceColor = matchColor(find(iGrps==k,1),:);
    p(k).MarkerEdgeColor = 'none';
    p(k).LineWidth = 2;
    p(k).MarkerSize = 9;
end
axis square;
xlim([-1.5 75]); ylim([-.02 1]);
xticks([0:15:75]); yticks([0:0.2:1]);
ylabel('Proportion matched');
xlabel('Days between sessions');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');

subplot(2,2,3);
hold on;
for k = 1:numel(uGrps)
    s1 = scatter(roisVtime(k).daysPostInj,roisVtime(k).nRoisNorm,108,roisVtime(k).matchColor,'filled');
    s1.MarkerEdgeColor = 'none'; %s1.MarkerFaceAlpha = 0.5;
end
h1 = boundedline(oRoi.x_fit,oRoi.y_fit,oRoi.delta,'k','alpha');
set(h1,'linewidth',2);
axis square;
xlim([-2 100]); ylim([-.02 1]);
xticks([0:20:100]); yticks([0:0.2:1]); xtickangle(0);
ylabel('Proportion of max ROIs');
xlabel('Days post-injection');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');
text(60,0.2,{sprintf('rsq = %0.3f',oRoi.rsq); sprintf('p = %0.3f',oRoi.corr_p)});

subplot(2,2,4);
hold on;
for k = 1:numel(uGrps)
    p(k) = plot(roisVtime(k).daysPostInj,roisVtime(k).nRoisNorm,'o-');
    p(k).Color = [roisVtime(k).matchColor 0.3];
    p(k).MarkerFaceColor = roisVtime(k).matchColor;
    p(k).MarkerEdgeColor = 'none';
    p(k).LineWidth = 2;
    p(k).MarkerSize = 9;
end  
axis square;
xlim([-2 100]); ylim([-.02 1]);
xticks([0:20:100]); yticks([0:0.2:1]); xtickangle(0);
ylabel('Proportion of max ROIs');
xlabel('Days post-injection');
set(gca,'LineWidth',2,'XColor','k','YColor','k','Layer', 'Top','FontSize',20,'TickDir','out');


if savebool
    figFold = fullfile(dataFold,'paperFigs','Fig2','Sub1_matchROIs');
    if ~exist(figFold,'dir'); mkdir(figFold); end
    figname = 'matchROIvsTime';
    savefig(figout, fullfile(figFold,[figname,'.fig']));
    saveas(figout,fullfile(figFold,[figname,'.png']));
    print(figout,fullfile(figFold,[figname,'.eps']), '-depsc', '-painters');
    matchROIs = matchData;
    save(fullfile(figFold,'matchROIvsTime.mat'),'matchROIs','oMatch','oRoi');
    
    figFold2 = fullfile(paperFold,'Figures','Fig2','Sub1_matchROIs');
    if ~exist(figFold2,'dir'); mkdir(figFold2); end
    savefig(figout, fullfile(figFold2,[figname,'.fig']));
    saveas(figout,fullfile(figFold2,[figname,'.png']));
    print(figout,fullfile(figFold2,[figname,'.eps']), '-depsc', '-painters');
    save(fullfile(figFold2,'matchROIvsTime.mat'),'matchROIs','oMatch','oRoi');
end

end