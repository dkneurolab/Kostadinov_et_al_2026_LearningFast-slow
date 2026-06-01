perC_allGood2 = perC_allGood2';

[movP, movTable, movStats] = kruskalwallis(perC_allGood2)
Rmov_stats = multcompare(movStats, 'Ctype','bonferroni')



percC_grp = cat(2,vAdaptationCurves.sessionsv5.percC)';

[percC_P, percC_Table, percC_Stats] = kruskalwallis(percC_grp)
percC_Stats_post = multcompare(percC_Stats, 'Ctype','bonferroni')

for iSesh = 1:numel(glmSummary)
    glmSummary(iSesh).pPCs_vis = glmSummary(iSesh).nSigVis/glmSummary(iSesh).nSig*100;
    glmSummary(iSesh).pPCs_mov = glmSummary(iSesh).nSigMov/glmSummary(iSesh).nSig*100;
    glmSummary(iSesh).pPCs_rew = glmSummary(iSesh).nSigRew/glmSummary(iSesh).nSig*100;
    glmSummary(iSesh).pPCs_lick = glmSummary(iSesh).nSigLick/glmSummary(iSesh).nSig*100;
end
lobV = glmSummary(ilobv);
sim2 = glmSummary(isim2);

[aVis,~] = signrank(cat(1,lobV.nSigVis)./cat(1,lobV.nSig),cat(1,sim2.nSigVis)./cat(1,sim2.nSig));
[aMov,~] = signrank(cat(1,lobV.nSigMov)./cat(1,lobV.nSig),cat(1,sim2.nSigMov)./cat(1,sim2.nSig));
[aRew,~] = signrank(cat(1,lobV.nSigRew)./cat(1,lobV.nSig),cat(1,sim2.nSigRew)./cat(1,sim2.nSig));
[aLic,~] = signrank(cat(1,lobV.nSigLick)./cat(1,lobV.nSig),cat(1,sim2.nSigLick)./cat(1,sim2.nSig));