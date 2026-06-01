function sesh = fig1_v4learningMovRsq(sesh)

movvC = mean(cat(2,sesh(end-1:end).movvC),2);
offxC = mean(cat(2,sesh(end-1:end).offxC),2);

for iSesh = 1:numel(sesh)
    Rmov0 = corrcoef([sesh(iSesh).movv movvC]);
    sesh(iSesh).RmovC = Rmov0(1:end-1,end);
    sesh(iSesh).RmovCmu = mean(sesh(iSesh).RmovC);
    Roffx0 = corrcoef([sesh(iSesh).offx offxC]);
    sesh(iSesh).RoffxC = Roffx0(1:end-1,end);
    sesh(iSesh).RoffxCmu = mean(sesh(iSesh).RoffxC);
end

end