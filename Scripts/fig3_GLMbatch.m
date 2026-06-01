function fig3_GLMbatch(dkdatabatch, dataFold, paperFold, comp, savebool)
%%
fprintf('\n---------- GLM analysis code ----------\n\n');
a = 0;
for i = 1:numel(dkdatabatch)
    
if dkdatabatch(i).good
    
a = a + 1;
%% These should come from data params input
inits.date              = dkdatabatch(i).date; % Experiment date
inits.name              = dkdatabatch(i).name; % Animal num
inits.paqnum            = dkdatabatch(i).paqnum; % Paq file num
inits.vrnum             = dkdatabatch(i).vrnum; % VR file num
if isfield(dkdatabatch(i),'fov')
    inits.fov           = dkdatabatch(i).fov;
else
    inits.fov           = [];
end
inits.dataFold          = dataFold;
inits.paperFold         = paperFold;
%% General settings - where to find data and whether to save
inits.computer          = comp; % 'mac' or 'pc' for 2 computers, 'npix' for probe  recordings
inits.savebool          = savebool; % Save data outputs - TRUE if wanting real data
inits.plotbool          = 0; % Plot example figs as analysis progresses - TRUE if testing, FALSE if running a lot of data
inits.nullbool          = 0;
inits.grpspks           = 1; % 0 for ungrouped spikes and 1 for grouped ones
inits.cosbool           = 0; % Make cosine glm plots
inits.boxbool           = 1; % Make boxcar glm plots

%% Specific settings - may change with data sets and tweaking of analysis [UNITS IN BRACKETS]
inits.paqfs             = []; % set inside 5e3; % [Hz] Sampling frequency for behaviour
inits.imfs              = []; % set inside 30; % [Hz] Sampling frequency for imaging        
inits.prebuffsec        = 1; % [s] Amount of time to take before trial onset
inits.postbuffsec       = 2; % [s] Amount of time to take before trial offset
inits.basisdur          = 1000; % [ms] Basis fxn duration for symmetrical filter sets i.e. 800 leads to 400 ms before and 400 ms after - keep this value even to get properly centered filters
inits.basiswid          = []; % set later - 1000/inits.imfs; % [ms] Basis fxn width/spacing - should divide evenly into basisdur - nBases = basisdur/basiswid/2
inits.lambda            = [.001 .01 .025 .05 .1 .25 .5 1 2.5 5 10 25 50 100]; % Range of lambda values for ridge regression
inits.xvalfold          = 10; % Cross-validation parameter
inits.cellgrping        = 'cells'; % 'zones' or 'cells';
inits.smoothspk         = false; % Smooth spikes with a causal gaussian? T/F
inits.smoothspksd       = 2; % Gaussian SD (in imaging frames);

inits.spkhistbool       = 0; % [BOOL] Whether to use spiking history as filter set
inits.basisdurSPK       = 800; % [ms] Duration of SPK history bases
inits.nBasesSPK         = 4; % [INT] Number of basis fxns to fit into basisdurSPK
inits.nloffsetSPK       = 2; % [frames] Number of frames to offset spk history (only for cosine bases)

%% Run GLM script
fig3_GLMmain(inits)

close all;

end

end