function [matchData, matchFOVs] = fig2_roiMatch(datasets,matchROIs,dataFold,paperFold,matchColors, savebool)
%% Concatenate sessions and
mouseStructs = fieldnames(datasets);
sessions = [];
for i = 1:numel(mouseStructs)
    sessions = [sessions datasets.(mouseStructs{i})]; %#ok<AGROW>
end
% sessionsv1 = sessions(cat(1,sessions.version) == 1);
% sessionsv4 = sessions(cat(1,sessions.version) == 4);
% sessionsv5 = sessions(cat(1,sessions.version) == 5);

alignFold = fullfile(dataFold,'01_Behav+imaging');
imFold = fullfile(dataFold,'02_Imaging');
behavFold = fullfile(dataFold,'04_Behav-processed');
figFold = fullfile(paperFold,'Figures','Fig4','Sub3_behavFigs');
if ~exist(figFold,'dir'); mkdir(figFold); end
figFold2 = fullfile(dataFold,'paperFigs','Fig4','Sub3_behavFigs');
if ~exist(figFold2,'dir'); mkdir(figFold2); end
% matchData = rmfield(matchData0,{'delDays','nMatched','nTotal','percMatched'});

%% Load datasets
if ~exist(fullfile(figFold2,'matchDataBehav.mat'),'file')
    
    matchFOVs = struct;
    ifov = false(size(sessions))';
    idate1 = ifov; idate2 = ifov; imouse = ifov;
    
    a = 0;
    for i = 1:numel(matchROIs)
        imouse(:) = false; idate1(:) = false; idate2(:) = false; ifov(:) = false; seshinds = [];
        for j = 1:numel(sessions)
            if strcmpi(sessions(j).name,matchROIs(i).name); imouse(j) = true; end
            if strcmpi(sessions(j).date,matchROIs(i).date1); idate1(j) = true; end
            if strcmpi(sessions(j).date,matchROIs(i).date2); idate2(j) = true; end
            if strcmpi(sessions(j).fov,matchROIs(i).fov); ifov(j) = true; end
        end
        
        try % If there is a session pair that meets criteria of ROI matching, populate
            sessionpair(1) = sessions(imouse & idate1 & ifov);
            sessionpair(2) = sessions(imouse & idate2 & ifov);
            seshinds = [find(imouse & idate1 & ifov); find(imouse & idate2 & ifov)];
        catch
            continue
        end
        if datenum(matchROIs(i).date1) > datenum(matchROIs(i).date2)
            sessionpair = sessionpair([2,1]);
            seshinds = seshinds([2,1]);
            tempdate = matchROIs(i).date1;
            matchROIs(i).date1 = matchROIs(i).date2;
            matchROIs(i).date2 = tempdate; clear tempdate
        end
        
        % Populate struct that has 1 entry for each FOV
        if i == 1
            a = a + 1;
            matchFOVs(a).name = matchROIs(i).name;
            matchFOVs(a).fov = matchROIs(i).fov;
            matchFOVs(a).v1date = [];
            matchFOVs(a).iv1sesh = [];
            matchFOVs(a).v4_1date = [];
            matchFOVs(a).iv4_1sesh = [];
            matchFOVs(a).v4_Ndate = [];
            matchFOVs(a).iv4_Nsesh = [];
            matchFOVs(a).v5date = [];
            matchFOVs(a).iv5sesh = [];
        elseif ~strcmpi(matchROIs(i).name,matchROIs(i-1).name) || ~strcmpi(matchROIs(i).fov,matchROIs(i-1).fov)
            a = a + 1;
            matchFOVs(a).name = matchROIs(i).name;
            matchFOVs(a).fov = matchROIs(i).fov;
        end
        
        %
        if strcmpi(matchROIs(i).fov,'lobv')
            matchFOVs(a).ilobv = true;
            matchFOVs(a).isim2 = false;
        elseif strcmpi(matchROIs(i).fov,'sim2')
            matchFOVs(a).ilobv = false;
            matchFOVs(a).isim2 = true;
        end
        
        % Populate matchFOV struct with dates, these should be in temporal order
        matchROIs(i).v1 = sessionpair(1).version; matchROIs(i).v2 = sessionpair(2).version;
        if sessionpair(1).version == 1 && sessionpair(2).version == 4 % v1+v4 pairs
            matchFOVs(a).v1date = sessionpair(1).date;
            matchFOVs(a).iv1sesh = seshinds(1);
            matchFOVs(a).v4_Ndate = sessionpair(2).date;
            matchFOVs(a).iv4_Nsesh = seshinds(2);
        elseif sessionpair(1).version == 4 && sessionpair(2).version == 4
            matchFOVs(a).v4_1date = sessionpair(1).date;
            matchFOVs(a).iv4_1sesh = seshinds(1);
            matchFOVs(a).v4_Ndate = sessionpair(2).date;
            matchFOVs(a).iv4_Nsesh = seshinds(2);
        elseif sessionpair(1).version == 4 && sessionpair(2).version == 5
            matchFOVs(a).v4_Ndate = sessionpair(1).date;
            matchFOVs(a).iv4_Nsesh = seshinds(1);
            matchFOVs(a).v5date = sessionpair(2).date;
            matchFOVs(a).iv5sesh = seshinds(2);
        end
        % Load data into to make figs
        matchData(i) = fig1_roiMatch_load(matchROIs(i),alignFold,imFold,behavFold,sessionpair); %#ok<AGROW>
        
    end
    matchData = rmfield(matchData,{'delDays','nMatched','nTotal','percMatched'});
    clear a i idate1 idate2 ifov imouse j seshinds sessionpair;
    
    matchFOVs = fig1_roiMatch_cleanup(matchFOVs);
    
    if savebool
        save(fullfile(figFold2,'matchDataBehav.mat'),'matchData','-v7.3');
        save(fullfile(figFold2,'matchDataFOVs.mat'),'matchFOVs');
        save(fullfile(figFold,'matchDataFOVs.mat'),'matchFOVs');
    end
    
else
    matchData = load(fullfile(figFold2,'matchDataBehav.mat')); matchData = matchData.matchData;
    matchFOVs = load(fullfile(figFold2,'matchDataFOVs.mat')); matchFOVs = matchFOVs.matchFOVs;
end

%% For each session pair, plot FOV and mean responses

for i = 1:numel(matchData)
    if ~isempty(matchData(i).name)
        figs = fig1_roiMatch_plot(matchData(i));
        if savebool
            figFoldfov1 = fullfile(figFold,sprintf('%s_%s_v%i_v%i',matchData(i).name,matchData(i).fov,matchData(i).v1,matchData(i).v2));
            if ~exist(figFoldfov1,'dir'); mkdir(figFoldfov1); end
            figFoldfov2 = fullfile(figFold2,sprintf('%s_%s_v%i_v%i',matchData(i).name,matchData(i).fov,matchData(i).v1,matchData(i).v2));
            if ~exist(figFoldfov2,'dir'); mkdir(figFoldfov2); end
            figname1 = '01_FOVs_raw';
            savefig(figs(1), fullfile(figFoldfov1,[figname1,'.fig']));
            saveas(figs(1),fullfile(figFoldfov1,[figname1,'.png']));
            print(figs(1),fullfile(figFoldfov1,[figname1,'.eps']), '-depsc', '-painters');
            savefig(figs(1), fullfile(figFoldfov2,[figname1,'.fig']));
            saveas(figs(1),fullfile(figFoldfov2,[figname1,'.png']));
            print(figs(1),fullfile(figFoldfov2,[figname1,'.eps']), '-depsc', '-painters');
            
            figname2 = '02_FOVs_matched';
            savefig(figs(2), fullfile(figFoldfov1,[figname2,'.fig']));
            saveas(figs(2),fullfile(figFoldfov1,[figname2,'.png']));
            print(figs(2),fullfile(figFoldfov1,[figname2,'.eps']), '-depsc', '-painters');
            savefig(figs(2), fullfile(figFoldfov2,[figname2,'.fig']));
            saveas(figs(2),fullfile(figFoldfov2,[figname2,'.png']));
            print(figs(2),fullfile(figFoldfov2,[figname2,'.eps']), '-depsc', '-painters');
            
            figname3 = '03_OnMovRewTmap';
            savefig(figs(3), fullfile(figFoldfov1,[figname3,'.fig']));
            saveas(figs(3),fullfile(figFoldfov1,[figname3,'.png']));
            print(figs(3),fullfile(figFoldfov1,[figname3,'.eps']), '-depsc', '-painters');
            savefig(figs(3), fullfile(figFoldfov2,[figname3,'.fig']));
            saveas(figs(3),fullfile(figFoldfov2,[figname3,'.png']));
            print(figs(3),fullfile(figFoldfov2,[figname3,'.eps']), '-depsc', '-painters');
            
            figname4 = '04_OnMovRewTmean';
            savefig(figs(4), fullfile(figFoldfov1,[figname4,'.fig']));
            saveas(figs(4),fullfile(figFoldfov1,[figname4,'.png']));
            print(figs(4),fullfile(figFoldfov1,[figname4,'.eps']), '-depsc', '-painters');
            savefig(figs(4), fullfile(figFoldfov2,[figname4,'.fig']));
            saveas(figs(4),fullfile(figFoldfov2,[figname4,'.png']));
            print(figs(4),fullfile(figFoldfov2,[figname4,'.eps']), '-depsc', '-painters');
            
            figname5 = '05_rewRrewTrewCmap';
            savefig(figs(5), fullfile(figFoldfov1,[figname5,'.fig']));
            saveas(figs(5),fullfile(figFoldfov1,[figname5,'.png']));
            print(figs(5),fullfile(figFoldfov1,[figname5,'.eps']), '-depsc', '-painters');
            savefig(figs(5), fullfile(figFoldfov2,[figname5,'.fig']));
            saveas(figs(5),fullfile(figFoldfov2,[figname5,'.png']));
            print(figs(5),fullfile(figFoldfov2,[figname5,'.eps']), '-depsc', '-painters');
            
            figname6 = '06_rewRrewTrewCmean';
            savefig(figs(6), fullfile(figFoldfov1,[figname6,'.fig']));
            saveas(figs(6),fullfile(figFoldfov1,[figname6,'.png']));
            print(figs(6),fullfile(figFoldfov1,[figname6,'.eps']), '-depsc', '-painters');
            savefig(figs(6), fullfile(figFoldfov2,[figname6,'.fig']));
            saveas(figs(6),fullfile(figFoldfov2,[figname6,'.png']));
            print(figs(6),fullfile(figFoldfov2,[figname6,'.eps']), '-depsc', '-painters');
            
            figname7 = '07_exTraces';
            savefig(figs(7), fullfile(figFoldfov1,[figname7,'.fig']));
            saveas(figs(7),fullfile(figFoldfov1,[figname7,'.png']));
            print(figs(7),fullfile(figFoldfov1,[figname7,'.eps']), '-depsc', '-painters');
            savefig(figs(7), fullfile(figFoldfov2,[figname7,'.fig']));
            saveas(figs(7),fullfile(figFoldfov2,[figname7,'.png']));
            print(figs(7),fullfile(figFoldfov2,[figname7,'.eps']), '-depsc', '-painters');
        end
    end
end

end