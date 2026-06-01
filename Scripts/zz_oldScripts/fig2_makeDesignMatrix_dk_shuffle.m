function [DMcells, DMdata] = fig2_makeDesignMatrix_dk_shuffle(basisfunct, trialIndices, filtfunc,nShuffles,iFieldShuf)
expt = basisfunct.expt;
bs = basisfunct.(filtfunc);
clear basisfunct

%% preallocate the matrix
totalT = sum([expt(trialIndices).duration]); % total duration of trials
endTrialIndices = [0, cumsum([expt(trialIndices).duration])];

fields = fieldnames(bs);
nBases = zeros(size(fields));
behavfields = false(size(fields));
for i = 1:numel(fields)
    j = 1;    
    if isempty(strfind(fields{i},'spk')); behavfields(i) = true; end
    fieldtemp = bs(j).(fields{i});    
    while isempty(fieldtemp)
        j = j+1;
        fieldtemp = bs(j).(fields{i});
    end
    nBases(i) = size(fieldtemp,2);    
end

totalnBases = sum(nBases);
DMcells = cell(nShuffles,1);
totaliBases = [0; cumsum(nBases)];

rng(69);
%% Place DM values in the right place in the big DM array
for iShuf = 1:nShuffles
    DM = zeros(totalT,totalnBases); % make it sparse?
    for j = 1:numel(trialIndices)
        for k = 1:numel(fields)
            fieldtemp = bs(trialIndices(j)).(fields{k});
            if ~isempty(fieldtemp) && iFieldShuf(k)
                fieldtemp = sum(fieldtemp,3); % Should not have anything more than 3d, so sum across the 3rd dimension
                ishift = randi(size(fieldtemp,1),1);
                fieldtemp_shift = circshift(fieldtemp,ishift,1);
                DMdim1 = endTrialIndices(j)+1:endTrialIndices(j+1);
                DMdim2 = totaliBases(k)+1:totaliBases(k+1);
                DM(DMdim1,DMdim2) = fieldtemp_shift;
            elseif ~isempty(fieldtemp)
                DMdim1 = endTrialIndices(j)+1:endTrialIndices(j+1);
                DMdim2 = totaliBases(k)+1:totaliBases(k+1);
                DM(DMdim1,DMdim2) = sum(fieldtemp,3);
            end
        end
    end
    DMcells{iShuf,1} = sparse(DM);
end

DMdata = v2struct(totalT,endTrialIndices,fields,behavfields,nBases,nShuffles,iFieldShuf);

end