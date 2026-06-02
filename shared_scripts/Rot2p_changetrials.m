function sets = Rot2p_changetrials(sets)
%% Put a definitive number on the change back trial for v5 analyses
if ~isfield(sets,'changeback')
    sets.changeback = 0;
    sets.changebacktrial = inf;
elseif sets.changeback == 0 % If there is no change back - set the trial num to inf
    sets.changebacktrial = inf;
else % If there is a change back (i.e. == 1), compute
    sets.changebacktrial = sets.changetrial*sets.changebackmult;
end

end