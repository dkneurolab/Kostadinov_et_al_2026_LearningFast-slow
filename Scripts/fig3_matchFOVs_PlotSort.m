function outStruct = fig3_matchFOVs_PlotSort(matchStruct,sortFOV,movBool)

if ~movBool || size(matchStruct(1).lick,2) < 2 || size(matchStruct(2).lick,2) < 2
    outStruct.spkMapSort = matchStruct(sortFOV).spkMapSort(end/2-14:end/2+15,:);
    outStruct.spkMap = cat(1,matchStruct(1).spkMap(end/2-29:end/2+30,:),matchStruct(2).spkMap(end/2-29:end/2+30,:));
    outStruct.spkMu = cat(1,matchStruct(1).spkMu(end/2-29:end/2+30),matchStruct(2).spkMu(end/2-29:end/2+30));
    outStruct.movMu = cat(1,matchStruct(1).movMu(end/2-99:end/2+100),matchStruct(2).movMu(end/2-99:end/2+100));
    outStruct.lickMu = cat(1,matchStruct(1).lickMu(end/2-99:end/2+100),matchStruct(2).lickMu(end/2-99:end/2+100));
else
    movMat = cat(2,matchStruct.mov);
    movMat = movMat(end/2-19:end/2+80,:);
    movCorr0 = corrcoef(movMat);
    movCorr0 = movCorr0(1:size(matchStruct(1).mov,2),size(matchStruct(1).mov,2)+1:end);
    movDist0 = squareform(pdist(movMat','cosine'));
    movDist0 = movDist0(1:size(matchStruct(1).mov,2),size(matchStruct(1).mov,2)+1:end);
    
    % Find movements that match the best and use them
    movCorr = movCorr0; movDist = movDist0;
    nMovs = min(size(movCorr));
    iMovs = zeros(nMovs,3); iMovs2 = iMovs;  
    for i = 1:nMovs
        corrMax = max(movCorr(:));
        iCorrMax = find(movCorr0(:) == corrMax,1);
        [iCorr1, iCorr2] = ind2sub(size(movCorr0),iCorrMax);
        iMovs(i,:) = [iCorr1 iCorr2 corrMax];
        movCorr(iCorr1,:) = 0; movCorr(:,iCorr2) = 0;
        
        distMin = min(movDist(:));
        iDistMin = find(movDist0(:) == distMin,1);
        [iDist1, iDist2] = ind2sub(size(movDist0),iDistMin);
        iMovs2(i,:) = [iDist1 iDist2 distMin];
        movDist(iDist1,:) = nan; movDist(:,iDist2) = nan;
        
    end    
    
    iMovs2 = iMovs2(iMovs2(:,3) < prctile(movDist0(:),20),:);
    [~,iMovs2sort] = sort(iMovs2(:,sortFOV),'ascend');
    iMovs2 = iMovs2(iMovs2sort,:);
    
    outStruct.spkMapSort = mean(matchStruct(sortFOV).spk(end/2-14:end/2+15,:,iMovs2(1:2:end,sortFOV)),3);
    outStruct.spkMap = cat(1,mean(matchStruct(1).spk(end/2-29:end/2+30,:,iMovs2(2:2:end,1)),3),...
        mean(matchStruct(2).spk(end/2-29:end/2+30,:,iMovs2(2:2:end,2)),3));
    outStruct.spkMu = mean(outStruct.spkMap,2);
    outStruct.movMu = cat(1,mean(matchStruct(1).mov(end/2-99:end/2+100,iMovs2(2:2:end,1)),2),...
        mean(matchStruct(2).mov(end/2-99:end/2+100,iMovs2(2:2:end,2)),2));
    outStruct.lickMu = cat(1,mean(matchStruct(1).lick(end/2-99:end/2+100,iMovs2(2:2:end,1)),2),...
        mean(matchStruct(2).lick(end/2-99:end/2+100,iMovs2(2:2:end,2)),2));

end

end