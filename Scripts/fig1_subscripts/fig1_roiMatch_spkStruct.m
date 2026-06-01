function outStruct = fig1_roiMatch_spkStruct(matchData,fname)
%% Make basic matrices of all data
fname1 = [fname,'1'];
fname2 = [fname,'2'];

outStruct(1).spk = 30*matchData.(fname1).spk(:,matchData.iROIs1,:);
outStruct(1).lick = matchData.(fname1).lick;
outStruct(1).mov = matchData.(fname1).mov;

outStruct(2).spk = 30*matchData.(fname2).spk(:,matchData.iROIs2,:);
outStruct(2).lick = matchData.(fname2).lick;
outStruct(2).mov = matchData.(fname2).mov; 

%% Decide who is reference dataset and sort ROIs based on activity in 1/2 of its trials

if matchData.tForm.tform_direction(1) == 1
    % If ref sesh is 1, use odd trials to sort even trials
    spk0 = mean(outStruct(1).spk(:,:,1:2:end),3,'omitnan');
    coeff = pca(spk0(end/2-14:end/2+15,:));
    [~, iSort] = sort(coeff(:,1),'descend');
    
    outStruct(1).spkMapSort = mean(outStruct(1).spk(:,iSort,1:2:end),3);
    outStruct(1).spkMap = mean(outStruct(1).spk(:,iSort,2:2:end),3);
    spkfov1 = mean(outStruct(1).spk(:,iSort,2:2:end),2);
    outStruct(1).spkMu = mean(spkfov1,3);
    outStruct(1).spkSD = std(spkfov1,[],3);
    outStruct(1).lickMu = mean(outStruct(1).lick(:,2:2:end),2);
    outStruct(1).lickSD = std(outStruct(1).lick(:,2:2:end),[],2);
    outStruct(1).movMu = mean(outStruct(1).mov(:,2:2:end),2);
    outStruct(1).movSD = std(outStruct(1).mov(:,2:2:end),[],2);
    
    % If ref sesh is 1, take full data from other sesh, but resort
    outStruct(2).spkMap = mean(outStruct(2).spk(:,iSort,:),3);
    spkfov2 = mean(outStruct(2).spk(:,iSort,:),2);
    outStruct(2).spkMu = mean(spkfov2,3);
    outStruct(2).spkSD = std(spkfov2,[],3);
    outStruct(2).lickMu = mean(outStruct(2).lick,2);
    outStruct(2).lickSD = std(outStruct(2).lick,[],2);
    outStruct(2).movMu = mean(outStruct(2).mov,2);
    outStruct(2).movSD = std(outStruct(2).mov,[],2);
else
    % If ref sesh is 2, use odd trials to sort even trials
    spk0 = mean(outStruct(2).spk(:,:,1:2:end),3,'omitnan');
    coeff = pca(spk0(end/2-14:end/2+15,:));
    [~, iSort] = sort(coeff(:,1),'descend');
    
    outStruct(2).spkMapSort = mean(outStruct(2).spk(:,iSort,1:2:end),3);
    outStruct(2).spkMap = mean(outStruct(2).spk(:,iSort,2:2:end),3);
    spkfov2 = mean(outStruct(2).spk(:,iSort,2:2:end),2);
    outStruct(2).spkMu = mean(spkfov2,3);
    outStruct(2).spkSD = std(spkfov2,[],3);
    outStruct(2).lickMu = mean(outStruct(2).lick(:,2:2:end),2);
    outStruct(2).lickSD = std(outStruct(2).lick(:,2:2:end),[],2);
    outStruct(2).movMu = mean(outStruct(2).mov(:,2:2:end),2);
    outStruct(2).movSD = std(outStruct(2).mov(:,2:2:end),[],2);
    
    % If ref sesh is 2, take full data from other sesh, but resort
    outStruct(1).spkMap = mean(outStruct(1).spk(:,iSort,:),3);
    spkfov1 = mean(outStruct(1).spk(:,iSort,:),2);
    outStruct(1).spkMu = mean(spkfov1,3);
    outStruct(1).spkSD = std(spkfov1,[],3);
    outStruct(1).lickMu = mean(outStruct(1).lick,2);
    outStruct(1).lickSD = std(outStruct(1).lick,[],2);
    outStruct(1).movMu = mean(outStruct(1).mov,2);
    outStruct(1).movSD = std(outStruct(1).mov,[],2);
end

end