function [GLMdata, GLMholdout,DMdataBoxSmall,trialbasis] = fig2_make_GLMdata_dk(DMboxdata, trialbasis, mzones, xvalfold,trialout)
%% Get rid of predictors we don't want

badFields = {'TrialOFFunder', 'TrialOFFover', 'RewTunder','RewTover','CueBeep','CueReward','RandReward'}';

iBadFields = false(numel(DMboxdata.fields),1);
for i = 1:numel(iBadFields)
    for j = 1:numel(badFields)
        if strcmpi(DMboxdata.fields{i},badFields{j})
            iBadFields(i) = true;
        end
    end
end

DMboxdata.fields = DMboxdata.fields(~iBadFields);
DMboxdata.behavfields = DMboxdata.behavfields(~iBadFields);
DMboxdata.nBases = DMboxdata.nBases(~iBadFields);

% trialbasis.bs = rmfield(trialbasis.bs,badFields);
trialbasis.bsbox = rmfield(trialbasis.bsbox,badFields);

%% Get indices of spk variables 

% trialsLU0 = trialout.trialsLU0;
trialsLC0 = trialout.trialsLC0;
% trialsLO0 = trialout.trialsLO0;

% Get indices of spk variables - assume that the spk vars are all at the end!!!!
% Note that we probably won't use these
DMboxspki = [0; cumsum(DMboxdata.nBases)];
DMboxspki = DMboxspki([~DMboxdata.behavfields; true]);

trialsALL0 = cat(1,trialbasis.expt.trialnum);

trialsGOOD = ismember(trialsALL0, trialsLC0);


trialbasis.expt = trialbasis.expt(trialsGOOD);
trialbasis.bs = trialbasis.bs(trialsGOOD);
trialbasis.bsbox = trialbasis.bsbox(trialsGOOD);

% if ~trialsGOOD(1)
%     trialbasis.bs(1).CueBeep = zeros(size(trialbasis.bs(1).TrialON));
%     trialbasis.bs(1).CueReward = zeros(size(trialbasis.bs(1).TrialON));
%     trialbasis.bs(1).CueReward = trialbasis.bs(1).CueReward(:,1:end/2+1);
%     trialbasis.bs(1).RandReward = zeros(size(trialbasis.bs(1).TrialON));
%     trialbasis.bsbox(1).CueBeep = zeros(size(trialbasis.bsbox(1).TrialON));
%     trialbasis.bsbox(1).CueReward = zeros(size(trialbasis.bsbox(1).TrialON));
%     trialbasis.bsbox(1).CueReward = trialbasis.bsbox(1).CueReward(:,1:end/2);
%     trialbasis.bsbox(1).RandReward = zeros(size(trialbasis.bsbox(1).TrialON));
% end

trialsALL = cat(1,trialbasis.expt.trialnum);
% trialsLUall = find(ismember(trialsALL,trialsLU0));
trialsLCall = find(ismember(trialsALL,trialsLC0));
% trialsLOall = find(ismember(trialsALL,trialsLO0));

starti = 2;

% trialsLU_holdout = trialsLUall(starti:3:end);
trialsLC_holdout = trialsLCall(starti:3:end);
% trialsLO_holdout = trialsLOall(starti:3:end);

% trialsLU = trialsLUall(~ismember(trialsLUall,trialsLU_holdout));
trialsLC = trialsLCall(~ismember(trialsLCall,trialsLC_holdout));
% trialsLO = trialsLOall(~ismember(trialsLOall,trialsLO_holdout));

iTrials0 = (1:numel(trialbasis.expt))';
% iHoldTrials = sort([trialsLU_holdout; trialsLC_holdout; trialsLO_holdout]);
iHoldTrials = sort(trialsLC_holdout,'ascend');
iTrials = iTrials0(~ismember(iTrials0,iHoldTrials));

% iRews0 = (1:numel(rewbasis.expt))';
% iHoldRews = iRews0(starti:3:end);
% iRews = iRews0(~ismember(iRews0,iHoldRews));

clear trialsALL0 trialsLU0 trialsLC0 trialsLO0 trialsLUall trialsLCall trialsLOall iTrials0 iRews0

%% Loop through and make x-val matrices

GLMdata = struct;
for i = 1:xvalfold
%     testU = trialsLU(i:xvalfold:end);
    testC = trialsLC(i:xvalfold:end);
%     testO = trialsLO(i:xvalfold:end);
%     testRew = iRews(i:xvalfold:end);
    
%     testtrials = sort([testU; testC; testO]);
    testtrials = sort(testC,'ascend');
    traintrials = iTrials(~ismember(iTrials,testtrials))';
%     trainRew = iRews(~ismember(iRews,testRew));
    
%     DMtest0trial = makeDesignMatrix_dk(trialbasis, testtrials, 'bs');
%     DMtrain0trial = makeDesignMatrix_dk(trialbasis, traintrials, 'bs');
    DMtestbox0trial = makeDesignMatrix_dk(trialbasis, testtrials, 'bsbox');
    DMtrainbox0trial = makeDesignMatrix_dk(trialbasis, traintrials, 'bsbox');
    
%     DMtest0rew = makeDesignMatrix_dk(rewbasis, testRew, 'bs');
%     DMtrain0rew = makeDesignMatrix_dk(rewbasis, trainRew, 'bs');
%     DMtestbox0rew = makeDesignMatrix_dk(rewbasis, testRew, 'bsbox');
%     DMtrainbox0rew = makeDesignMatrix_dk(rewbasis, trainRew, 'bsbox');
    
%     DMtest0 = [DMtest0trial; DMtest0rew];
%     DMtrain0 = [DMtrain0trial; DMtrain0rew];
%     DMtestbox0 = [DMtestbox0trial; DMtestbox0rew];
%     DMtrainbox0 = [DMtrainbox0trial; DMtrainbox0rew];
    DMtestbox0 = [DMtestbox0trial];
    DMtrainbox0 = [DMtrainbox0trial];
    
%     DMtest = DMtest0(:,1:sum(DMdata.nBases(DMdata.behavfields)));
%     DMtrain = DMtrain0(:,1:sum(DMdata.nBases(DMdata.behavfields)));
    DMtestbox = DMtestbox0(:,1:sum(DMboxdata.nBases(DMboxdata.behavfields)));
    DMtrainbox = DMtrainbox0(:,1:sum(DMboxdata.nBases(DMboxdata.behavfields)));
    
    GLMdata(i).testtrials = testtrials;
    GLMdata(i).traintrials = traintrials;
%     GLMdata(i).testRew = testRew;
%     GLMdata(i).trainRew = trainRew;
    

    for j = 1:numel(mzones)        

        GLMdata(i).YtestTrial(:,j) = getSpikes_dk(trialbasis,testtrials,mzones{j});
        GLMdata(i).YtrainTrial(:,j) = getSpikes_dk(trialbasis,traintrials,mzones{j});
%         GLMdata(i).YtestRew(:,j) = getSpikes_dk(rewbasis,testRew,mzones{j});
%         GLMdata(i).YtrainRew(:,j) = getSpikes_dk(rewbasis,trainRew,mzones{j});
        
%         GLMdata(i).Ytest(:,j) = [GLMdata(i).YtestTrial(:,j); GLMdata(i).YtestRew(:,j)];
%         GLMdata(i).Ytrain(:,j) = [GLMdata(i).YtrainTrial(:,j); GLMdata(i).YtrainRew(:,j)];
        GLMdata(i).Ytest(:,j) = GLMdata(i).YtestTrial(:,j);
        GLMdata(i).Ytrain(:,j) = GLMdata(i).YtrainTrial(:,j);
        
        if numel(DMboxspki) > 1 % this is pprobably obsolete
%             GLMdata(i).DMtest(:,:,j)        = [DMtest DMtest0(:,DMspki(j)+1:DMspki(j+1))];
%             GLMdata(i).DMtrain(:,:,j)       = [DMtrain DMtrain0(:,DMspki(j)+1:DMspki(j+1))];
            GLMdata(i).DMtestbox(:,:,j)     = [DMtestbox DMtestbox0(:,DMboxspki(j)+1:DMboxspki(j+1))];
            GLMdata(i).DMtrainbox(:,:,j)    = [DMtrainbox DMtrainbox0(:,DMboxspki(j)+1:DMboxspki(j+1))];
        else
%             GLMdata(i).DMtest       = DMtest;
%             GLMdata(i).DMtrain      = DMtrain;
            GLMdata(i).DMtestbox    = DMtestbox;
            GLMdata(i).DMtrainbox   = DMtrainbox;
        end
    end
end

clear DMtest DMtest0 DMtest0rew DMtest0trial DMtestbox DMtestbox0 DMtestbox0rew DMtestbox0trial DMtrain DMtrain0 DMtrain0rew ...
    DMtrain0trial DMtrainbox DMtrainbox0 DMtrainbox0rew DMtrainbox0trial testC testO testU testRew testtrials trainRew traintrials

%% Make DM for withheld data

% GLMholdout.iTrials = v2struct(trialsALL, trialsGOOD,trialsLU, trialsLU_holdout, trialsLC, trialsLC_holdout, trialsLO, trialsLO_holdout,...
%     iTrials, iHoldTrials, iRews, iHoldRews);
GLMholdout.iTrials = v2struct(trialsALL, trialsGOOD, trialsLC, trialsLC_holdout, iTrials, iHoldTrials);


% [DMholdTrial, DMholdData] = makeDesignMatrix_dk(trialbasis, iHoldTrials, 'bs');
[DMholdBoxTrial, DMholdBoxData] = makeDesignMatrix_dk(trialbasis, iHoldTrials, 'bsbox');
% [DMholdRew] = makeDesignMatrix_dk(rewbasis, iHoldRews, 'bs');
% [DMholdBoxRew] = makeDesignMatrix_dk(rewbasis, iHoldRews, 'bsbox');

% GLMholdout.DMtrial_full = DMholdTrial;
GLMholdout.DMtrialBox_full = DMholdBoxTrial;
% GLMholdout.DMrew_full = DMholdRew;
% GLMholdout.DMrewBox_full = DMholdBoxRew;
% GLMholdout.DMdata_full = DMholdData;
GLMholdout.DMtrialBoxData_full = DMholdBoxData;

for j = 1:numel(mzones)
    GLMholdout.Ytrial(:,j) = getSpikes_dk(trialbasis,iHoldTrials,mzones{j});
%     GLMholdout.Yrew(:,j) = getSpikes_dk(rewbasis,iHoldRews,mzones{j});
end

[DMdataBoxSmall,iBoxSmall] = reducePredictorsBox(DMholdBoxData);

GLMholdout.DMtrialBox = GLMholdout.DMtrialBox_full(:,iBoxSmall);
% GLMholdout.DMrewBox = GLMholdout.DMrewBox_full(:,iBoxSmall);

for i = 1:xvalfold
    GLMdata(i).DMtestbox_full = GLMdata(i).DMtestbox;
    GLMdata(i).DMtestbox = GLMdata(i).DMtestbox(:,iBoxSmall);    
    GLMdata(i).DMtrainbox_full = GLMdata(i).DMtrainbox;
    GLMdata(i).DMtrainbox = GLMdata(i).DMtrainbox(:,iBoxSmall);    
end

end