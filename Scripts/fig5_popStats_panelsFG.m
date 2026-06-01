load('C:\Users\Dimitar\Dropbox\DK papers\2022_Kostadinov_FastSlowLearning\Figures\Fig5\Sub2_behavSummary\v5_adaptSummary.mat')
imfs = 30;

% Lob V mov
lvMov_baseFOV = arrayfun(@(s) mean(s.spkMovBaseFR), lobvdata)*imfs;
lvMov_earlyFOV = arrayfun(@(s) mean(s.spkMovEarlyFR), lobvdata)*imfs;
lvMov_midFOV = arrayfun(@(s) mean(s.spkMovMidFR), lobvdata)*imfs;
lvMov_lateFOV = arrayfun(@(s) mean(s.spkMovLateFR), lobvdata)*imfs;
lvMov_washFOV = arrayfun(@(s) mean(s.spkMovWashFR), lobvdata)*imfs;
lvMov_FOV = cat(1,lvMov_baseFOV,lvMov_earlyFOV,lvMov_midFOV,lvMov_lateFOV,lvMov_washFOV)';
[p,tbl,stats] = friedman(lvMov_FOV,1);
multcompare(stats)

% Lob sim mov
lsMov_baseFOV = arrayfun(@(s) mean(s.spkMovBaseFR), sim2data)*imfs;
lsMov_earlyFOV = arrayfun(@(s) mean(s.spkMovEarlyFR), sim2data)*imfs;
lsMov_midFOV = arrayfun(@(s) mean(s.spkMovMidFR), sim2data)*imfs;
lsMov_lateFOV = arrayfun(@(s) mean(s.spkMovLateFR), sim2data)*imfs;
lsMov_washFOV = arrayfun(@(s) mean(s.spkMovWashFR), sim2data)*imfs;
lsMov_FOV = cat(1,lsMov_baseFOV,lsMov_earlyFOV,lsMov_midFOV,lsMov_lateFOV,lsMov_washFOV)';
[p,tbl,stats] = friedman(lsMov_FOV,1);
multcompare(stats)

% Lob V rew
lvRew_baseFOV = arrayfun(@(s) mean(s.spkRewBaseFR), lobvdata)*imfs;
lvRew_earlyFOV = arrayfun(@(s) mean(s.spkRewEarlyFR), lobvdata)*imfs;
lvRew_midFOV = arrayfun(@(s) mean(s.spkRewMidFR), lobvdata)*imfs;
lvRew_lateFOV = arrayfun(@(s) mean(s.spkRewLateFR), lobvdata)*imfs;
lvRew_washFOV = arrayfun(@(s) mean(s.spkRewWashFR), lobvdata)*imfs;
lvRew_FOV = cat(1,lvRew_baseFOV,lvRew_earlyFOV,lvRew_midFOV,lvRew_lateFOV,lvRew_washFOV)';
[p,tbl,stats] = friedman(lvRew_FOV,1);
multcompare(stats)

% Lob sim rew
lsRew_baseFOV = arrayfun(@(s) mean(s.spkRewBaseFR), sim2data)*imfs;
lsRew_earlyFOV = arrayfun(@(s) mean(s.spkRewEarlyFR), sim2data)*imfs;
lsRew_midFOV = arrayfun(@(s) mean(s.spkRewMidFR), sim2data)*imfs;
lsRew_lateFOV = arrayfun(@(s) mean(s.spkRewLateFR), sim2data)*imfs;
lsRew_washFOV = arrayfun(@(s) mean(s.spkRewWashFR), sim2data)*imfs;
lsRew_FOV = cat(1,lsRew_baseFOV,lsRew_earlyFOV,lsRew_midFOV,lsRew_lateFOV,lsRew_washFOV)';
[p,tbl,stats] = friedman(lsRew_FOV,1);
multcompare(stats)




