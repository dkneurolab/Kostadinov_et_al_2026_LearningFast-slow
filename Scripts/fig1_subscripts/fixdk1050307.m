function v5out = fixdk1050307(v5in)

%% Put in basic info that will not require much fudging

v5out.name = v5in(1).name;
v5out.date = v5in(1).date(1:10);

v5out.params.trialnum = v5in(1).params.trialnum + v5in(2).params.trialnum;
maxtrials1 = max(v5in(1).params.trialinds);
v5out.params.trialinds = [v5in(1).params.trialinds; v5in(2).params.trialinds + maxtrials1];
v5out.params.trialsU = [v5in(1).params.trialsU; v5in(2).params.trialsU + maxtrials1];
v5out.params.trialsC = [v5in(1).params.trialsC; v5in(2).params.trialsC + maxtrials1];
v5out.params.trialsO = [v5in(1).params.trialsO; v5in(2).params.trialsO + maxtrials1];
v5out.params.blocksize = v5in(1).params.blocksize;
v5out.params.rhoU = [v5in(1).params.rhoU; v5in(2).params.rhoU];
v5out.params.rhoC = [v5in(1).params.rhoC; v5in(2).params.rhoC];
v5out.params.rhoO = [v5in(1).params.rhoO; v5in(2).params.rhoO];
v5out.params.rhoUmu = mean(v5out.params.rhoU);
v5out.params.rhoUsd = std(v5out.params.rhoU);
v5out.params.rhoCmu = mean(v5out.params.rhoC);
v5out.params.rhoCsd = std(v5out.params.rhoC);
v5out.params.rhoOmu = mean(v5out.params.rhoO);
v5out.params.rhoOsd = std(v5out.params.rhoO);

numblocks = 0;
alltrials = v5out.params.trialinds;
trialsC = v5out.params.trialsC;
for i = 1:v5out.params.blocksize:max(v5out.params.trialinds)
    numblocks = numblocks + 1;
    rollrange = [i,i+v5out.params.blocksize-1];    
    numC(numblocks,:) = [sum(trialsC >= rollrange(1) & trialsC <= rollrange(2)) sum(alltrials >= rollrange(1) & alltrials <= rollrange(2))];
end

v5out.params.numC = numC;
v5out.params.xendall = [v5in(1).params.xendall; v5in(2).params.xendall];

v5out.offxall = cat(2,v5in.offxall);
v5out.offxmu = cat(2,v5in.offxmu);
v5out.offxsd = cat(2,v5in.offxsd);
v5out.movvall = cat(2,v5in.movvall);
v5out.movvmu = cat(2,v5in.movvmu);
v5out.movvsd = cat(2,v5in.movvsd);
v5out.percC = cat(1,v5in.percC);
v5out.xend = cat(1,v5in.xend);
v5out.rhomovmu = cat(1,v5in.rhomovmu);
v5out.rhomovsd = cat(1,v5in.rhomovsd);

a = 4;

end