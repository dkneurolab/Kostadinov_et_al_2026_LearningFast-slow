function fig3_Supp_eventRateSD(seshim, dataFold)
% Calculate SD of event amplitudes for each cell and FOV

imfs = 30;

spkRate = [];
spkSD = [];

for iSesh = 1:numel(seshim)
    spkmat = load(fullfile(dataFold,'01_Behav+imaging', sprintf('Version_%i',seshim(iSesh).version), seshim(iSesh).name, sprintf('%s_%s',seshim(iSesh).date, seshim(iSesh).fov),'spkmat.mat'));
    
    spkmat = spkmat.spkmat;

    nCells = size(spkmat,2);
    spkSD_local = zeros(nCells,1);
    spkRate_local = zeros(nCells,1);
    for jCell = 1:nCells
        spks0 = spkmat(spkmat(:,jCell) > 0,jCell);        
        spks = spks0/mean(spks0);
        spkSD_local(jCell) = std(spks);
        spkRate_local(jCell) = numel(spks)/size(spkmat,1)*imfs;
    end    
    spkRate = [spkRate; spkRate_local];
    spkSD = [spkRate; spkSD_local];

    if mod(iSesh,5) == 0
        fprintf(['Done ',num2str(iSesh),' of ',num2str(numel(seshim)), ' sessions.\n'])
    end
end

fprintf('Event detection stats across all sessions:\n')
fprintf('  Rate mean:   %.3f Hz\n', mean(spkRate))
fprintf('  Rate std:    %.3f Hz\n', std(spkRate))
fprintf('  Amp SD mean: %.3f Hz\n', mean(spkSD))
fprintf('  Amp SD std:  %.3f Hz\n', std(spkSD))


end