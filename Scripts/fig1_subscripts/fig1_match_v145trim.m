function dataout = fig1_match_v145trim(datain,iSort,fs,offset)
%%
spk1 = cat(2,datain(1:3:end).spkMap);
spk1 = spk1(offset*fs.imfs+[end/2-fs.imfs+1:end/2+fs.imfs],iSort);
spk4 = cat(2,datain(2:3:end).spkMap);
spk4 = spk4(offset*fs.imfs+[end/2-fs.imfs+1:end/2+fs.imfs],iSort);
spk5 = cat(2,datain(3:3:end).spkMap);
spk5 = spk5(offset*fs.imfs+[end/2-fs.imfs+1:end/2+fs.imfs],iSort);

dataout.spkMap = fs.imfs*cat(1,spk1,spk4,spk5);

%%
spkMu0 = cat(2,datain.spkMu);
spkMu0 = spkMu0(offset*fs.imfs+[end/2-fs.imfs+1:end/2+fs.imfs],:);
spkMu1 = spkMu0(:,1:3:end);
spkMu4 = spkMu0(:,2:3:end);
spkMu5 = spkMu0(:,3:3:end);

dataout.spkMu = fs.imfs*cat(1,spkMu1, spkMu4, spkMu5);

%%
lickMu0 = cat(2,datain.lickMu);
lickMu0 = lickMu0(offset*fs.behavfs+[end/2-fs.behavfs+1:end/2+fs.behavfs],:);
lickMu1 = lickMu0(:,1:3:end);
lickMu4 = lickMu0(:,2:3:end);
lickMu5 = lickMu0(:,3:3:end);

dataout.lickMu = cat(1,lickMu1, lickMu4, lickMu5);

%%
movMu0 = cat(2,datain.movMu);
movMu0 = movMu0(offset*fs.behavfs+[end/2-fs.behavfs+1:end/2+fs.behavfs],:);
movMu1 = movMu0(:,1:3:end);
movMu4 = movMu0(:,2:3:end);
movMu5 = movMu0(:,3:3:end);

dataout.movMu = cat(1,movMu1, movMu4, movMu5);

end