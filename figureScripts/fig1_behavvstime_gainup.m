function statout = fig1_behavvstime_gainup(savebool)

%% Load in all datasets:

% savebool = 0;
numgrp = 3;
maxnum = 180;
blocksize = 20;

parentfolder = '/Users/dimitar/Desktop/CF_learning_paper/04_Behav-processed';
cd(parentfolder);

targetfolder = '/Users/dimitar/Desktop/CF_learning_paper/paperFigs/Fig1';
sub1 = fullfile(targetfolder,'Sub1_percMovCorr');
if ~exist(sub1,'dir'); mkdir(sub1); end

v5sessions0 = load('v5_analysis_sessions_gainup.mat'); v5sessions0 = v5sessions0.v5all_stab;

% v5sessions0 = v5sessions0([2:end]);

a = 0; 
for i = 1:numel(v5sessions0)
    if sum([~strcmp(v5sessions0(i).name,'DK105'), ~strcmp(v5sessions0(i).date(1:10),'2018-03-07')]) > 0
        a = a + 1;
        v5s(a) = v5sessions0(i);        
    elseif strcmp(v5sessions0(i).name,'DK105') && strcmp(v5sessions0(i).date,'2018-03-07a')
    	a = a + 1;
        v5s(a) = fixdk1050307(v5sessions0(i:i+1));        
    end
end

gainmean = [1.4*1+1.1*3+1.2*2+1.3*2+1*1+0.9*2+2*1+2*1.1+2*1]/17;
gaindiv = 37.5/360*4.5*pi/gainmean;
% gaindiv = 1;

%% Reindex trials for each session so gain change happens at the same trial number!!
gaintrial = zeros(size(v5s))';
gaintrial(1) = 100; % DK052
gaintrial(2) = 100; % DK063
gaintrial(3) = 100; % DK070
gaintrial(4) = 80; % DK103 0223
gaintrial(5) = 80; % DK103 0303
gaintrial(6) = 80; % DK103 0306
gaintrial(7) = 80; % DK103 0320
gaintrial(8) = 119; % DK105 0307 - 109 trials from first recording and 10 baselines from 2nd
gaintrial(9) = 80; % DK105 0307 - 109 trials from first recording and 10 baselines from 2nd
gaintrial(10) = 60; %
gaintrial(11) = 60; %
gaintrial(12) = 60; %
gaintrial(13) = 60; %
gaintrial(14) = 60; %
gaintrial(15) = 60; %
gaintrial(16) = 60; %
gaintrial(17) = 60; %

changetrial = min(gaintrial);

for i = 1:numel(v5s)
    v5s2(i).name = v5s(i).name;
    v5s2(i).date = v5s(i).date;
    
    newtrialinds = v5s(i).params.trialinds-gaintrial(i)+changetrial > 0;
    v5s2(i).params.trialnum = sum(newtrialinds);
    trialindsold = v5s(i).params.trialinds(newtrialinds);
    v5s2(i).params.trialinds = v5s(i).params.trialinds(newtrialinds)-gaintrial(i)+changetrial;
    
    v5s2(i).params.trialsU = v5s(i).params.trialsU(ismember(v5s(i).params.trialsU,trialindsold))-gaintrial(i)+changetrial;
    v5s2(i).params.trialsC = v5s(i).params.trialsC(ismember(v5s(i).params.trialsC,trialindsold))-gaintrial(i)+changetrial;
    v5s2(i).params.trialsO = v5s(i).params.trialsO(ismember(v5s(i).params.trialsO,trialindsold))-gaintrial(i)+changetrial;
    
    v5s2(i).params.rhoU = v5s(i).params.rhoU(ismember(v5s(i).params.trialsU,trialindsold));
    v5s2(i).params.rhoC = v5s(i).params.rhoC(ismember(v5s(i).params.trialsC,trialindsold));
    v5s2(i).params.rhoO = v5s(i).params.rhoO(ismember(v5s(i).params.trialsO,trialindsold));
    
    
    v5s2(i).offxall = v5s(i).offxall(:,newtrialinds);
    v5s2(i).movvall = v5s(i).movvall(:,newtrialinds);
    
end
numsession = numel(v5s2);

%% Make matrix of all sessions and trials

numsess = zeros(size(v5s2))';

for i = 1:numel(numsess); numsess(i) = max(v5s2(i).params.trialinds); end

maxsess = max(numsess);

offxlate = zeros(400,maxsess,numel(numsess))*nan;
movvlate = offxlate;
movxlate = offxlate;
Rlate =  zeros(maxsess,numel(numsess))*nan;
Clate = Rlate;

% Collect statout
statout.sessions0 = v5s;
statout.sessions = v5s2;

%% Populate big matrices with data from each session

for i = 1:numel(numsess)
    offxtemp = v5s2(i).offxall;
    movvtemp = -v5s2(i).movvall*gaindiv;
    movxtemp = cumsum(v5s2(i).movvall)/100+1;
    trialindtemp = v5s2(i).params.trialinds;
    
    offxlate(:,trialindtemp,i) = offxtemp;
    movvlate(:,trialindtemp,i) = movvtemp;
    movxlate(:,trialindtemp,i) = movxtemp;
    
    Rlate(v5s2(i).params.trialsU,i) = v5s2(i).params.rhoU;
    Rlate(v5s2(i).params.trialsC,i) = v5s2(i).params.rhoC;
    Rlate(v5s2(i).params.trialsO,i) = v5s2(i).params.rhoO;
    
    Clate(v5s2(i).params.trialsU,i) = 0;
    Clate(v5s2(i).params.trialsC,i) = 1;
    Clate(v5s2(i).params.trialsO,i) = 0;
end

% Collect statout
statout.movxlateall = movxlate;
statout.movvlateall = movvlate;
statout.offxlateall = offxlate;

%% Now make averages

offxlatemu =  nanmean(offxlate,3);
offxlatesd =  nanstd(offxlate,[],3);
movvlatemu =  nanmean(movvlate,3);
movvlatesd =  nanstd(movvlate,[],3);
movxlatemu =  nanmean(movxlate,3);
movxlatesd =  nanstd(movxlate,[],3);

% Rlatemu = nanmean(Rlate,2);
% Rlatesd = nanstd(Rlate,[],2)/sqrt(size(Rlate,2));
Clatemu = nanmean(Clate,2);
Clatesd = nanstd(Clate,[],2)/sqrt(size(Clate,2));

numblocks = ceil(size(offxlatemu,2)/blocksize);
offxblockmu = zeros(size(offxlatemu,1)/2,numblocks);
offxblocksd = offxblockmu; 
movvblockmu = offxblockmu; movvblocksd = offxblockmu;
movxblockmu = offxblockmu; movxblocksd = offxblockmu;

for i = 1:numblocks
    blockinds = ((i-1)*blocksize+1:i*blocksize)';
    if blockinds(end) > size(offxlatemu,2); blockinds = ((i-1)*blocksize+1:size(offxlatemu,2))'; end
    offxblockmu(:,i) = mean(offxlatemu(1:200,blockinds),2);
    offxblocksd(:,i) = mean(offxlatesd(1:200,blockinds),2);
    movvblockmu(:,i) = mean(movvlatemu(101:300,blockinds),2);
    movvblocksd(:,i) = mean(movvlatesd(101:300,blockinds),2);
    movxblockmu(:,i) = mean(movxlatemu(101:300,blockinds),2);
    movxblocksd(:,i) = mean(movxlatesd(101:300,blockinds),2);
end

offxblockmuplot = offxblockmu(:);
offxblocksdplot = offxblocksd(:)/sqrt(numsession);
offxblockmuplot_gain = offxblockmuplot;
offxblockmuplot_gain(200*changetrial/blocksize+1:end) = (offxblockmuplot_gain(200*changetrial/blocksize+1:end)-1)*1.6+1;
offxblocksdplot_gain = offxblocksdplot;
offxblocksdplot_gain(200*changetrial/blocksize+1:end) = offxblocksdplot_gain(200*changetrial/blocksize+1:end)*1.6;

movxblockmuplot = movxblockmu(:);
movxblocksdplot = movxblocksd(:)/sqrt(numsession);

offxendmu = nanmean(squeeze(nanmean(offxlate(end/2-9:end/2,:,:),1)),2);
offxendsd = nanstd(squeeze(nanmean(offxlate(end/2-9:end/2,:,:),1)),[],2)/sqrt(numsession);
offxendmu_gain = offxendmu;
offxendmu_gain(changetrial+1:end) = (offxendmu_gain(changetrial+1:end)-1)*1.6+1;
offxendsd_gain = offxendsd;
offxendsd_gain(changetrial+1:end) = offxendsd_gain(changetrial+1:end)*1.6;

movvblockmu2 = movvblockmu;
movvbase = nanmean(nanmean(movvlate(101:300,1:changetrial,:),2),3);
for i = 1:numblocks
    movvblockmu2(:,i) = movvblockmu2(:,i) - movvbase;
end
movvblockmuplot = movvblockmu(:);
movvblockmuplot2 = movvblockmu2(:);
movvblocksdplot = movvblocksd(:)/sqrt(numsession);

% Plot endpoints (single trials)
offxendmu2 = zeros(maxnum/numgrp,1); offxendsd2 = offxendmu2; offxendmu_gain2 = offxendmu2; offxendsd_gain2 = offxendmu2; 
Clatemu2 = offxendmu2; Clatesd2 = offxendmu2;
for i = 1:maxnum/numgrp
    offxendmu2(i) = mean(offxendmu((i-1)*numgrp+1:i*numgrp));
    offxendsd2(i) = mean(offxendsd((i-1)*numgrp+1:i*numgrp));
    offxendmu_gain2(i) = mean(offxendmu_gain((i-1)*numgrp+1:i*numgrp));
    offxendsd_gain2(i) = mean(offxendsd_gain((i-1)*numgrp+1:i*numgrp));
    
    Clatemu2(i) = mean(Clatemu((i-1)*numgrp+1:i*numgrp));
    Clatesd2(i) = mean(Clatesd((i-1)*numgrp+1:i*numgrp));
end
toff2 = (numgrp:numgrp:maxnum)' - numgrp/2;

% Collect statout
statout.Clate = Clate;
statout.Clatemu = Clatemu;
statout.Clatesd = Clatesd;
statout.Clatemu2 = Clatemu2;
statout.Clatesd2 = Clatesd2;

statout.offxblockmu = offxblockmu;
statout.offxblocksd = offxblocksd;
statout.offxblockmuplot = offxblockmuplot;
statout.offxblocksemplot = offxblocksdplot;
statout.offxblockmuplot_gain = offxblockmuplot_gain;
statout.offxblocksemplot_gain = offxblocksdplot_gain;

statout.offxendmu = offxendmu;
statout.offxendsem = offxendsd;
statout.offxendmu2 = offxendmu2;
statout.offxendsem2 = offxendsd2;
statout.offxendmu_gain = offxendmu_gain;
statout.offxendsem_gain = offxendsd_gain;
statout.offxendmu_gain2 = offxendmu_gain2;
statout.offxendsem_gain2 = offxendsd_gain2;

statout.movxblockmu = movxblockmu;
statout.movxblocksd = movxblocksd;
statout.movxblockmuplot = movxblockmuplot;
statout.movxblocksemplot = movxblocksdplot;

statout.movvblockmu = movvblockmu;
statout.movvblocksd = movvblocksd;
statout.movvblockmuplot = movvblockmuplot;
statout.movvblocksemplot = movvblocksdplot;

statout.movvblockmudelta = movvblockmu2;
statout.movvblocksddelta = movvblocksd;
statout.movvblockmudeltaplot = movvblockmuplot2;
statout.movvblocksemdeltaplot = movvblocksdplot;

%% stats:
subset = 1:17;
abcd = squeeze(nanmean(offxlate(end/2-9:end/2,:,:),1));
% abcd(changetrial+1:end,:) = (abcd(changetrial+1:end,:)-1)*1.6+1;
abcd2 = [nanmean(abcd(1:60,subset))' nanmean(abcd(61:80,subset))' nanmean(abcd(161:180,subset))'];
[~,~,stats1] = kruskalwallis(abcd2); %,'displayopt','off');
% [~,~,stats1] = anova1(abcd2); %,'displayopt','off');
cex1 = multcompare(stats1,'Ctype','bonferroni');

statout.cex = cex1;

%% Plotting
maxt = maxnum/blocksize*2;

f1 = figure;
% Plot velocity-movement (blocked)
subplot(3,2,1);
plot([0 maxt],[0 0],'--k');
b1 = boundedline((1:numel(movvblockmuplot))/100,movvblockmuplot,movvblocksdplot,'k','alpha');
set(b1,'linewidth',1.5); box off;
set(gca,'LineWidth',1,'XColor','k','YColor','k','Layer', 'Top','FontSize',12,'XGrid','on');
axis([0 maxt -2.5 12.5]);
yticks([0:2.5:10]);
xticks([0:2:maxt]);
xticklabels([1,blocksize:blocksize:maxnum]);
ylabel('Wheel velocity (cm/s)');
xlabel('Trial number (blocked)');

% Plot positions-movement (blocked)
subplot(3,2,2);
rectangle('Position',[0 -1/3 maxt 2/3],'FaceColor',[0 0 .5 .1],'EdgeColor','none')
axis([0 maxt -1.25 1.25]); hold on;
plot([0 maxt],[0 0],'--k');
b2(1) = boundedline((1:numel(offxblocksdplot))/100,offxblockmuplot,offxblocksdplot,'r','alpha');
set(b2(1),'linewidth',1.5); box off;
b2(2) = boundedline((1:numel(offxblocksdplot))/100,offxblockmuplot_gain,offxblocksdplot_gain,'k','alpha');
set(b2(2),'linewidth',1.5); box off;
set(gca,'LineWidth',1,'XColor','k','YColor','k','Layer', 'Top','FontSize',12,'XGrid','on');
yticks([-1:.5:1]);
xticks([0:2:maxt]);
xticklabels([1,blocksize:blocksize:maxnum]);
ylabel('Object position - movement offset');
xlabel('Trial number (blocked)');

% Plot endpoints (single trials)
subplot(3,2,4);
rectangle('Position',[0 -1/3 maxnum 2/3],'FaceColor',[0 0 .5 .1],'EdgeColor','none')
axis([0 maxnum -1.25 1.25]); hold on;
plot([0 maxnum],[0 0],'--k');
b4(1) = boundedline((1:numel(offxendmu)),offxendmu,offxendsd,'r','alpha');
set(b4(1),'linewidth',1.5); box off;
b4(2) = boundedline((1:numel(offxendmu)),offxendmu_gain,offxendsd_gain,'k','alpha');
set(b4(2),'linewidth',1.5); box off;
set(gca,'LineWidth',1,'XColor','k','YColor','k','Layer', 'Top','FontSize',12,'XGrid','on');
yticks([-1:.5:1]);
xticks([0:blocksize:maxnum]);
xticklabels([1,blocksize:blocksize:maxnum]);
ylabel('Object endpoint');
xlabel('Trial number');

% Plot velocities - delta
subplot(3,2,3);
b3 = boundedline([1:numel(movvblockmuplot)]'/100,movvblockmuplot2,movvblocksdplot,'k');
set(b3,'linewidth',1.5); box off;
set(gca,'LineWidth',1,'XColor','k','YColor','k','Layer', 'Top','FontSize',12,'XGrid','on');
axis([0 maxt -2.5 2.5]);
yticks([-2:2]);
xticks([0:2:maxt]);
xticklabels([1,blocksize:blocksize:maxnum]);
ylabel({'Wheel velocity';'delta (cm/s)'});
xlabel('Trial number (blocked)');

% Plot percent correct
subplot(3,2,5);
b5 = boundedline(1:numel(Clatemu),Clatemu,Clatesd,'cmap',[0 .5 0],'alpha');
set(b5,'linewidth',1.5); box off;
set(gca,'LineWidth',1,'XColor','k','YColor','k','Layer', 'Top','FontSize',12,'XGrid','on');
ylim([0 1]);
yticks([0:.2:1]);
xlim([0 maxnum]);
xticks([0:blocksize:maxnum]);
xticklabels([1,blocksize:blocksize:maxnum]);
ylabel('Percent correct');
xlabel('Trial number');

mousenames0 = cat(1,v5s2.name);
mousenames = unique(mousenames0,'rows');
suptitle(['Gain change behaviour: n = ',num2str(numel(v5s2)),' sessions from ',num2str(size(mousenames,1)),' mice.']);

if savebool
savename = ['02_Gain_change_behaviour'];
savefig(f1,fullfile(sub1,[savename,'.fig']),'compact')
saveas(f1,fullfile(sub1,[savename,'.png']));
print(fullfile(sub1,[savename,'.eps']), '-depsc', '-painters');
end

%% Replot with trials binned for endpoint and percent correct:

f2 = figure;
% Plot positions-movement (blocked)
subplot(3,2,1);
rectangle('Position',[0 -1/3 maxt 2/3],'FaceColor',[0 0 .5 .1],'EdgeColor','none')
axis([0 maxt -1.25 1.25]); hold on;
plot([0 maxt],[0 0],'--k');
bb1 = boundedline((1:numel(movxblockmuplot))/100,movxblockmuplot,movxblocksdplot,'k','alpha');
set(bb1,'linewidth',1.5); box off;
set(gca,'LineWidth',1,'XColor','k','YColor','k','Layer', 'Top','FontSize',12,'XGrid','on');
yticks([-1:.5:1]);
xticks([0:2:maxt]);
xticklabels([1,blocksize:blocksize:maxnum]);
ylabel('Object position - movement onset');
xlabel('Trial number (blocked)');

subplot(3,2,2);
rectangle('Position',[0 -1/3 maxt 2/3],'FaceColor',[0 0 .5 .1],'EdgeColor','none')
axis([0 maxt -1.25 1.25]); hold on;
plot([0 maxt],[0 0],'--k');
bb2(1) = boundedline((1:numel(offxblocksdplot))/100,offxblockmuplot,offxblocksdplot,'r','alpha');
set(bb2(1),'linewidth',1.5); box off;
bb2(2) = boundedline((1:numel(offxblocksdplot))/100,offxblockmuplot_gain,offxblocksdplot_gain,'k','alpha');
set(bb2(2),'linewidth',1.5); box off;
set(gca,'LineWidth',1,'XColor','k','YColor','k','Layer', 'Top','FontSize',12,'XGrid','on');
yticks([-1:.5:1]);
xticks([0:2:maxt]);
xticklabels([1,blocksize:blocksize:maxnum]);
ylabel('Object position - movement offset');
xlabel('Trial number (blocked)');

subplot(3,2,4);
rectangle('Position',[0 -1/3 maxnum 2/3],'FaceColor',[0 0 .5 .1],'EdgeColor','none')
axis([0 maxnum -1.25 1.25]); hold on;
plot([0 maxnum],[0 0],'--k');
bb4(1) = boundedline(toff2,offxendmu2,offxendsd2,'r','alpha');
set(bb4(1),'linewidth',1.5); box off;
bb4(2) = boundedline(toff2,offxendmu_gain2,offxendsd_gain2,'k','alpha');
set(bb4(2),'linewidth',1.5); box off;
set(gca,'LineWidth',1,'XColor','k','YColor','k','Layer', 'Top','FontSize',12,'XGrid','on');
yticks([-1:.5:1]);
xticks([0:blocksize:maxnum]);
xticklabels([1,blocksize:blocksize:maxnum]);
ylabel('Object endpoint');
xlabel('Trial number');

% Plot velocities
subplot(3,2,3);
bb3 = boundedline([1:numel(movvblockmuplot2)]'/100,movvblockmuplot2,movvblocksdplot,'k');
set(bb3,'linewidth',1.5); box off;
set(gca,'LineWidth',1,'XColor','k','YColor','k','Layer', 'Top','FontSize',12,'XGrid','on');
axis([0 maxt -2.5 2.5]);
yticks([-2:2]);
xticks([0:2:maxt]);
xticklabels([1,blocksize:blocksize:maxnum]);
ylabel({'Wheel velocity';'delta (cm/s)'});
xlabel('Trial number (blocked)');

% Plot percent correct
subplot(3,2,5);
bb5 = boundedline(toff2,Clatemu2,Clatesd2,'cmap',[0 .5 0],'alpha');
set(bb5,'linewidth',1.5); box off;
set(gca,'LineWidth',1,'XColor','k','YColor','k','Layer', 'Top','FontSize',12,'XGrid','on');
ylim([0 1]);
yticks([0:.2:1]);
xlim([0 maxnum]);
xticks([0:blocksize:maxnum]);
xticklabels([1,blocksize:blocksize:maxnum]);
ylabel('Percent correct');
xlabel('Trial number');

suptitle(['Gain change behaviour: n = ',num2str(numel(v5s2)),' sessions from ',num2str(size(mousenames,1)),' mice.']);

cd(targetfolder);
if savebool
    savename = ['02_Gain_change_behaviour_binned'];
    savefig(f2,fullfile(sub1,[savename,'.fig']),'compact')
    saveas(f2,fullfile(sub1,[savename,'.png']));
    print(fullfile(sub1,[savename,'.eps']), '-depsc', '-painters');
    
    save(fullfile(targetfolder,'Gaindata.mat'),'statout');
end

cd(sub1);
end
