function outstruct = Rot2p_behavstability(trials0,inits,sets,savebool)
% Rot2p_behavstability(paqdata.behav.trials,inits,sets,1);
%% Set initial condtions for things

if isfield(sets,'rotgain')
    gaindiv = 37.5/360*4.5*pi/sets.rotgain;
else
    gaindiv = 37.5/360*4.5*pi; % Set gain to 1 if it's not there, may want to fix eventually
end

paqfs = inits.paqfs;
mintime = inits.mintime;

realtrials = logical(cat(1,trials0.realtrial));
trials = trials0(realtrials);

outcomes = cat(1,trials.outcome);

LUstruct = trials(logical(outcomes(:,1)));
LCstruct = trials(logical(outcomes(:,2)));
LOstruct = trials(logical(outcomes(:,3)));

%% In earlier version, was offsetting trialoffset by reward delay, but it is better not to do this, offx just ends at end of movement!
rewoffset = 0;

% if isfield(sets,'rewarddelay')
%     rewoffset = round(sets.rewarddelay/60*inits.paqfs);
% else
%     rewdelay = zeros(size(LCstruct))';
%     for i = 1:numel(rewdelay)        
%         % 1 second before start and end of trial should be enough
%         % do not have to index specifically, because it is relative timing
%         % but indices need tobe the same for trialbool and rewbool!
%         trialbool = LCstruct(i).data(inits.paqfs*(inits.mintime-1):end-inits.paqfs*(inits.mintime-1),3); 
%         rewbool = LCstruct(i).data(inits.paqfs*(inits.mintime-1):end-inits.paqfs*(inits.mintime-1),4);
%         trialoffi = threshdetect_dk(trialbool,.5,-1);
%         rewbooli = threshdetect_dk(rewbool,.5,11);
%         rewbooli = rewbooli((rewbooli > trialoffi));
%         rewdelay(i) = rewbooli(1)-trialoffi;        
%     end
%     % Note that we need to add a frame here because of they way that virmen
%     % updates things - it sends the trialoff trigger 1 frame later.
%     rewoffset = round((mean(rewdelay)/inits.paqfs+1/60)*inits.paqfs,-2);
% end
    
%% Build average behavioural responses - position, velocity, correlation for trial types

% Position data
LCoffx = double(cat(2,LCstruct.offx));
Lofft = (1:size(LCoffx,1))'/inits.paqfs; Lofft = Lofft(1:end/2)-2; % halfway here is the reward time, so take everything before this
LCoffxmu = mean(LCoffx,2); LCoffxmu = LCoffxmu(rewoffset+(1:end/2));
LCoffxsd = std(LCoffx,[],2); LCoffxsd = LCoffxsd(rewoffset+(1:end/2));

LUoffx = double(cat(2,LUstruct.offx));
LUoffxmu = mean(LUoffx,2); LUoffxmu = LUoffxmu(rewoffset+(1:end/2));
LUoffxsd = std(LUoffx,[],2); LUoffxsd = LUoffxsd(rewoffset+(1:end/2));

LOoffx = double(cat(2,LOstruct.offx));
LOoffxmu = mean(LOoffx,2); LOoffxmu = LOoffxmu(rewoffset+(1:end/2));
LOoffxsd = std(LOoffx,[],2); LOoffxsd = LOoffxsd(rewoffset+(1:end/2));

% Velocity data
LCmovv = double(cat(2,LCstruct.movv)*-gaindiv);
Lmovt = (1:size(LCmovv,1))'/inits.paqfs; Lmovt = Lmovt(end/4+1:end/4*3)-2; % halfway here is the movement onset time, so take 1 s around it
LCmovvmu = mean(LCmovv,2); LCmovvmu = LCmovvmu(end/4+1:end/4*3);
LCmovvsd = std(LCmovv,[],2); LCmovvsd = LCmovvsd(end/4+1:end/4*3);

LUmovv = double(cat(2,LUstruct.movv)*-gaindiv);
LUmovvmu = mean(LUmovv,2); LUmovvmu = LUmovvmu(end/4+1:end/4*3);
LUmovvsd = std(LUmovv,[],2); LUmovvsd = LUmovvsd(end/4+1:end/4*3);

LOmovv = double(cat(2,LOstruct.movv)*-gaindiv);
LOmovvmu = mean(LOmovv,2); LOmovvmu = LOmovvmu(end/4+1:end/4*3);
LOmovvsd = std(LOmovv,[],2); LOmovvsd = LOmovvsd(end/4+1:end/4*3);

% Correlation data
corrrange = 8e3:15e3;
rhoC = zeros(size(LCmovv,2),1);
for i = 1:size(LCmovv,2)
    rhodummy = corrcoef(LCmovvmu(corrrange-5e3),LCmovv(corrrange,i));
    rhoC(i) = rhodummy(2,1);
end
rhoCmu = mean(rhoC); rhoCsd = std(rhoC);
rhoU = zeros(size(LUmovv,2),1);
for i = 1:size(LUmovv,2)
    rhodummy = corrcoef(LCmovvmu(corrrange-5e3),LUmovv(corrrange,i));
    rhoU(i) = rhodummy(2,1);    
end
rhoUmu = mean(rhoU); rhoUsd = std(rhoU);
rhoO = zeros(size(LOmovv,2),1);
for i = 1:size(LOmovv,2)
    rhodummy = corrcoef(LCmovvmu(corrrange-5e3),LOmovv(corrrange,i));
    rhoO(i) = rhodummy(2,1);    
end
rhoOmu = mean(rhoO); rhoOsd = std(rhoO);


%% Makre figure left panels (summaries)
numcols = 5;
f1 = figure;

% Object position
subplot(3,numcols,1);
rectangle('Position',[-2 -1/3 2 2/3],'FaceColor',[0 0 .5 .1],'EdgeColor','none')
axis([-2 0 -1.25 1.25]); hold on;
plot([-2 0],[0 0],'--k');
hx(1) = boundedline(Lofft,LUoffxmu,LUoffxsd,'cmap',[1 0 0]*.75,'alpha');
set(hx(1),'linewidth',1.5);
hx(2) = boundedline(Lofft,LOoffxmu,LOoffxsd,'cmap',[0 1 1]*.75,'alpha');
set(hx(2),'linewidth',1.5);
hx(3) = boundedline(Lofft,LCoffxmu,LCoffxsd,'cmap','k','alpha');
set(hx(3),'linewidth',1.5);
ylabel('Object position');
xlabel('Time from wheel stop (s)');
set(gca,'LineWidth',1,'XColor','k','YColor','k','Layer', 'Top','FontSize',12);

% Wheel velocity
subplot(3,numcols,numcols+1);
hv(1) = boundedline(Lmovt,LUmovvmu,LUmovvsd,'cmap',[1 0 0]*.75,'alpha');
set(hv(1),'linewidth',1.5);
hv(2) = boundedline(Lmovt,LOmovvmu,LOmovvsd,'cmap',[0 1 1]*.75,'alpha');
set(hv(2),'linewidth',1.5);
hv(3) = boundedline(Lmovt,LCmovvmu,LCmovvsd,'cmap','k','alpha');
set(hv(3),'linewidth',1.5);
axis([-1 1 -2.5 22.5]);
ylabel({'Wheel velocity';'(cm/s)'});
xlabel('Time from wheel movement (s)');
set(gca,'LineWidth',1,'XColor','k','YColor','k','Layer', 'Top','FontSize',12);

% Correlations
subplot(3,numcols,2*numcols+1); 
yyaxis left; hold on;
errorbar(1,rhoUmu,rhoUsd,'o','CapSize',0,'LineWidth',1,'Color',[1 0 0]*.75)
errorbar(2,rhoCmu,rhoCsd,'o','CapSize',0,'LineWidth',1,'Color',[0 0 0]*.75)
errorbar(3,rhoOmu,rhoOsd,'o','CapSize',0,'LineWidth',1,'Color',[0 1 1]*.75)
xticks([1:3]);
ylabel({'Correlation with mean';'correct movement'});
xlabel('Trial outcome');
xticklabels({'Under.','Corr.','Over.'});
axis([.5 3.5 0 1]);
set(gca,'LineWidth',1,'XColor','k','YColor','k','Layer', 'Top','FontSize',12);

outUCO = sum(outcomes);
yyaxis right;
bar(1,outUCO(1),.5,'FaceColor',[.75 0 0],'EdgeColor','none','FaceAlpha',.3);
bar(2,outUCO(2),.5,'FaceColor',[0 0 0],'EdgeColor','none','FaceAlpha',.3);
bar(3,outUCO(3),.5,'FaceColor',[0 .75 .75],'EdgeColor','none','FaceAlpha',.3);
ylim([0 sum(outUCO)]);
% ylabel(['Number of trials']);
set(gca,'LineWidth',1,'XColor','k','YColor','k','Layer', 'Top','FontSize',12);

%% Concatenate trialstructs and reorder

% With new input 'trials', everything is in the correct order
trialinds = cat(1,trials.trialnum);
% trialsU = cat(1,LUstruct.trialind);
trialsC = cat(1,LCstruct.trialnum);
% trialsO = cat(1,LOstruct.trialind);

blocksize = 20;

numC = zeros(ceil(max(trialinds)/blocksize),2); blockis = numC;
offxmu = []; offxsd = []; movvmu = []; movvsd = []; rhomovmu = []; rhomovsd = []; xend = [];

numblocks = 0;
for i = 1:blocksize:max(trialinds)
    rollrange = [i,i+blocksize-1];
    numblocks = numblocks + 1;
    blockis(numblocks,:) = rollrange;
    substruct = trials(trialinds >= rollrange(1) & trialinds <= rollrange(2));
    numC(numblocks,:) = [sum(trialsC >= rollrange(1) & trialsC <= rollrange(2)) numel(substruct)];
    offx = double(cat(2,substruct.offx)); offx = offx(rewoffset+(1:end/2),:);
    offxmu = [offxmu; decimate(mean(offx,2),50)]; %#ok<*AGROW>
    offxsd = [offxsd; decimate(std(offx,[],2),50)]; 
    xend = [xend; mean(offx(end-399:end,:))'];
    
    movv = double(cat(2,substruct.movv)); movv = movv(end/4+1:end/4*3,:)*-gaindiv;
    movvmu = [movvmu; decimate(mean(movv,2),50)]; %#ok<*AGROW>
    movvsd = [movvsd; decimate(std(movv,[],2),50)];
    
    rhomov = zeros(size(movv,2),1);
    for j = 1:size(movv,2)
        rhodummy = corrcoef(LCmovvmu(corrrange-paqfs),movv(corrrange-paqfs,j));
        rhomov(j) = rhodummy(2,1);    
    end
    rhomovmu = [rhomovmu; mean(rhomov)]; rhomovsd = [rhomovsd; std(rhomov)];
end

maxt = ceil(numel(offxmu)/100);

% Plot positions
subplot(3,numcols,2:numcols);
rectangle('Position',[0 -1/3 maxt 2/3],'FaceColor',[0 0 .5 .1],'EdgeColor','none')
axis([0 maxt -1.25 1.25]); hold on;
plot([0 maxt],[0 0],'--k');
b1 = boundedline([1:numel(offxmu)]'/100,offxmu,offxsd,'k','alpha');
set(b1,'linewidth',1.5); box off;
set(gca,'LineWidth',1,'XColor','k','YColor','k','Layer', 'Top','FontSize',12,'XGrid','on');
yticks([-1:.5:1]);
xticks([0:2:maxt]);
xticklabels([]);

% Plot velocities
subplot(3,numcols,numcols+(2:numcols));
b2 = boundedline([1:numel(movvmu)]'/100,movvmu,movvsd,'k');
set(b2,'linewidth',1.5); box off;
set(gca,'LineWidth',1,'XColor','k','YColor','k','Layer', 'Top','FontSize',12,'XGrid','on');
axis([0 maxt -2.5 22.5]);
yticks([0:5:20]);
xticks([0:2:maxt]);
xticklabels([]);

% Plot correlations
subplot(3,numcols,numcols*2+(2:numcols));
yyaxis left;
e1 = errorbar(2*(1:numel(rhomovmu))-1,rhomovmu,rhomovsd,'o','CapSize',0,'LineWidth',2,'Color',[0 0 0],'LineStyle','none');
e1.MarkerFaceColor = 'k';
ylim([0 1]);
yticks([0:.2:1]);
set(gca,'LineWidth',1,'XColor','k','YColor','k','Layer', 'Top','FontSize',12);
yyaxis right;
plot(1:2:2*numblocks,numC(:,1)./numC(:,2),'o','Color',[0 .5 0],'MarkerFaceColor',[0 .5 0])
ylim([0 1]);
yticks([0:.2:1]);
ylabel(['Percent correct per block']);
set(gca,'LineWidth',1,'XColor','k','YColor',[0 .5 0],'Layer', 'Top','FontSize',12,'XGrid','on','YGrid','on');
xlim([0 maxt]);
% axis([0 maxt .4 1]);
xticks(0:2:maxt);
xticklabels([]);
box off
xlabel(['Trial block size: ',num2str(blocksize)]);

suptitle([inits.mousename,' on ',inits.date]);

%% Save fig maybe
if savebool
    savename = [inits.mousename,'_',inits.date,'_behavstability'];
    savefig(f1,fullfile(inits.targetfolder,[savename,'.fig']),'compact')
    saveas(f1,fullfile(inits.targetfolder,[savename,'.png']));
    print(fullfile(inits.targetfolder,[savename,'.eps']), '-depsc', '-painters');
end

%% Package output structure
outstruct.name = inits.mousename;
outstruct.date = inits.date;
outstruct.params.trialnum = numel(trials);
outstruct.params.trialinds = trialinds;
outstruct.params.trialsU = cat(1,LUstruct.trialnum);
outstruct.params.trialsC = cat(1,LCstruct.trialnum);
outstruct.params.trialsO = cat(1,LOstruct.trialnum);
outstruct.params.blocksize = blocksize;
outstruct.params.rhoU = rhoU;
outstruct.params.rhoC = rhoC;
outstruct.params.rhoO = rhoO;
outstruct.params.rhoUmu = rhoUmu;
outstruct.params.rhoUsd = rhoUsd;
outstruct.params.rhoCmu = rhoCmu;
outstruct.params.rhoCsd = rhoCsd;
outstruct.params.rhoOmu = rhoOmu;
outstruct.params.rhoOsd = rhoOsd;
outstruct.params.numC = numC;
outstruct.params.xendall = unique(xend,'stable');
outstruct.params.gaindiv = gaindiv;


offxall0 = double(cat(2,trials.offx));
movvall0 = double(cat(2,trials.movv));

offxall = zeros(size(offxall0,1)/50,size(offxall0,2)); movvall = offxall;
for k = 1:size(offxall0,2)
    offxall(:,k) = decimate(offxall0(:,k),50);
    movvall(:,k) = decimate(movvall0(:,k),50);
end

outstruct.offxall = offxall;
outstruct.offxmu = reshape(offxmu,[],numblocks);
outstruct.offxsd = reshape(offxsd,[],numblocks);
outstruct.movvall = movvall;
outstruct.movvmu = reshape(movvmu,[],numblocks);
outstruct.movvsd = reshape(movvsd,[],numblocks);
outstruct.percC = numC(:,1)./numC(:,2);
outstruct.xend = mean(outstruct.offxmu(151:160,:))';
outstruct.rhomovmu = rhomovmu;
outstruct.rhomovsd = rhomovsd;


end
