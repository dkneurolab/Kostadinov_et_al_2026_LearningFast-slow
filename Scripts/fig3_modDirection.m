function fig3_modDirection(matchData,glmParams,dataFold,paperFold,savebool)
% Calculate percentage of neurons that are movement and reward activated
% and suppressed
imfs = 30;

%% Load in sig summary based on info:
figFold0 = fullfile(paperFold,'Figures','Fig3','Sub3_SummaryBars');
figFold = fullfile(figFold0,sprintf('%s_%s_alpha%i',glmParams.dataType,glmParams.distr,glmParams.alpha*100));

sigSummary_v4n = load(fullfile(figFold,'sigSummary_v4n.mat'));
glmData = sigSummary_v4n.glmSummary;

%% Loop through and pull out data
glmNames = fieldnames(glmData);
modDirs = rmfield(glmData, glmNames(11:end));

for iFov = 1:numel(glmData)
    % Look for matching session and extract mean traces per cell
    for jMatch = 1:numel(matchData)        
        if strcmpi(modDirs(iFov).name, matchData(jMatch).name) & strcmpi(modDirs(iFov).fov, matchData(jMatch).fov) ...
                && strcmpi(modDirs(iFov).date, matchData(jMatch).date1)
            modDirs(iFov).movData = mean(matchData(jMatch).mov1.spk(:,modDirs(iFov).sigBoolMov,:),3);
            modDirs(iFov).rewData = mean(matchData(jMatch).rewT1.spk(:,modDirs(iFov).sigBoolRew,:),3);
        elseif strcmpi(modDirs(iFov).name, matchData(jMatch).name) & strcmpi(modDirs(iFov).fov, matchData(jMatch).fov) ...
                && strcmpi(modDirs(iFov).date, matchData(jMatch).date2)            
            modDirs(iFov).movData = mean(matchData(jMatch).mov2.spk(:,modDirs(iFov).sigBoolMov,:),3);
            modDirs(iFov).rewData = mean(matchData(jMatch).rewT2.spk(:,modDirs(iFov).sigBoolRew,:),3);
        end
    end

    % Now go through cell-wise and compute if movement activated or suppressed:
    movBase = mean(modDirs(iFov).movData(1:imfs,:),1)'*imfs;
    movMod = mean(modDirs(iFov).movData(imfs*2-imfs*0.3+1:imfs*2,:),1)'*imfs;
    modDirs(iFov).movUp = movMod > movBase;
    modDirs(iFov).movDown = movMod < movBase;

    % Now go through cell-wise and compute if reward activated or suppressed:
    rewBase = mean(modDirs(iFov).rewData(imfs*2-imfs*0.5+1:imfs*2,:),1)'*imfs;
    rewMod = mean(modDirs(iFov).rewData(imfs*2+1:imfs*2+0.2*imfs,:),1)'*imfs;
    modDirs(iFov).rewUp = rewMod > rewBase;
    modDirs(iFov).rewDown = rewMod < rewBase;

    fprintf('nMov up: %i of %i\n',sum(modDirs(iFov).movUp),numel(modDirs(iFov).movUp))
    fprintf('nRew up: %i of %i\n',sum(modDirs(iFov).rewUp),numel(modDirs(iFov).rewUp))
end

modDirs_lobv = modDirs(cat(1,modDirs.iLobV));
modDirs_sim2 = modDirs(cat(1,modDirs.iSim2));

movUp_lobv = cat(1,modDirs_lobv.movUp);
movDown_lobv = cat(1,modDirs_lobv.movDown);
fprintf('LV: nMov up: %i of %i\n',sum(movUp_lobv),numel(movUp_lobv));
fprintf('LV: nMov down: %i of %i\n',sum(movDown_lobv),numel(movDown_lobv));
movUp_sim2 = cat(1,modDirs_sim2.movUp);
movDown_sim2 = cat(1,modDirs_sim2.movDown);
fprintf('LS: nMov up: %i of %i\n',sum(movUp_sim2),numel(movUp_sim2));
fprintf('LS: nMov down: %i of %i\n',sum(movDown_sim2),numel(movDown_sim2));
rewUp_lobv = cat(1,modDirs_lobv.rewUp);
rewDown_lobv = cat(1,modDirs_lobv.rewDown);
fprintf('LV: nRew up: %i of %i\n',sum(rewUp_lobv),numel(rewUp_lobv));
fprintf('LV: nRew down: %i of %i\n',sum(rewDown_lobv),numel(rewDown_lobv));
rewUp_sim2 = cat(1,modDirs_sim2.rewUp);
rewDown_sim2 = cat(1,modDirs_sim2.rewDown);
fprintf('LS: nRew up: %i of %i\n',sum(rewUp_sim2),numel(rewUp_sim2));
fprintf('LS: nRew down: %i of %i\n',sum(rewDown_sim2),numel(rewDown_sim2));



%% Save!
if savebool
    a = 4

end

end
