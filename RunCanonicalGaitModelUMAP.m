%% RunCanonicalGaitModelUMAP.m


%% Set parameters

% Select whether a new simulation should be run if a cache exists
% Note that this should be set to true if a new tripodFraction is set since
% the current code does not check for caches with different tripodFractions
runSim = false;

% Set the sample size
n = 1e5;

% Set the fraction of canonical tripod gait
tripodFraction = 1/3;

%%  Set plotting options

% Colormaps
gaitCmp = [0.3467,0.5360,0.6907;0.9153,0.2816,0.2878;0.4416,0.7490,0.4322;1.0000,0.5984,0.2000];
cmpFreq = copper(2^16);
cmpLimb = jet(2^16);
cmpDown = 0.8 * jet(7);

% Viewport for 3D plotting
viewport = [-37.5, 30];

%% Run simulation or load cached results

% Add necessary folders to path
addpath('utils','models');

% Search for a cache file
cacheFilePath = dir('cache/canonical_gait_data_*.mat');

% Check if a cache file exists
% If one does, load the simulation data from the cache.
% If not, re-run the simulation
if isempty(cacheFilePath) || runSim
    
    % Run the simulation
    [ dataTable, gaitInd, freq, L, D, parpoolTocBytes ] = GenerateCanonicalGaitModelData(n, tripodFraction);
    
    % Cache the results
    if ~isfolder('cache'), mkdir('cache'); end
    spath =  sprintf('canonical_gait_data_n%d_%.2fpct_tripod_%s.mat', n, tripodFraction, datestr(now, 'yyyy-mm-dd_HH-MM-SS'));
    save(fullfile('cache',spath), 'dataTable','gaitInd','freq','L','D', 'n', 'tripodFraction', 'parpoolTocBytes', '-v7.3');
    
else
    % Load data from cache
    cacheFilePath = fullfile(cacheFilePath.folder, cacheFilePath.name);
    fprintf('Loading cached simulation results from %s\n', cacheFilePath);
    load(cacheFilePath, 'L','D','freq', 'gaitInd', 'tripodFraction');
end

%% Preprocess the limb data

% Preprocess limb data
tic;
L = L - nanmean(L, 1);
X_L = reshape(permute(L, [3 1 2]), [], size(L,1)*size(L,2));
X_L = (X_L - nanmean(X_L,1))./nanstd(X_L,0,1);
fprintf('Preprocessed limb data in %f seconds.\n', toc);

%% Run UMAP

% Write data to file
tic;
datapath = fullfile('cache', 'umap_data.h5');
if isfile(datapath)
    delete(datapath);
end
h5create(datapath, '/x_l', size(X_L));
h5write(datapath, '/x_l', X_L);
h5create(datapath, '/limb_components', 1);
h5write(datapath, '/limb_components', 3);
fprintf('Wrote data to h5 file %s in %f seconds\n', datapath, toc);

%% Load UMAP results

tic;
resultpath = fullfile('cache', 'umap_result.h5');
limbUMAP = h5read(resultpath, '/embedding_l');
limbUMAP = limbUMAP';
fprintf('Read result from h5 file %s in %f seconds\n', resultpath, toc);

%% Plot embedding

% Define the labels
gaitList = {'tripod','left tetrapod','right tetrapod', 'wave'};

% Plot the stepping frequency distributions for each gait
xq = (0:0.1:14)';

% Compute the distributions
f = nan(length(xq)-1,4);
for ind = 1:4
    f(:,ind) = histcounts(freq(gaitInd==ind), xq, 'normalization','pdf');
end

% Plot the distributions
MakeFigure;
hold on;
set(gca, 'colororder', gaitCmp);
plot(xq(1:end-1), f, 'linewidth', 2);
hold off;
legend(gaitList);
axis('square');
xlabel('frequency (Hz)');
ylabel('pdf (s)');
ConfAxis('fontSize', 14);

% Scatter embedding by gait
MakeFigure;
scattern(limbUMAP, 10, gaitInd, 'filled');
axis('equal');
colormap(gaitCmp);
cbar = colorbar();
caxis([1 4]);
cbar.Ticks =  (1+(4-1)/(2*4)):((4-1)/4):4;
cbar.TickLabels = gaitList;
ylabel(cbar, 'gait type');
xlabel('u_1 (a.u.)');
ylabel('u_2 (a.u.)');
zlabel('u_3 (a.u.)');
ConfAxis('fontSize', 14);
view(viewport);

% Scatter embedding by frequency
MakeFigure;
scattern(limbUMAP, 10, freq, 'filled');
axis('equal');
caxis([0 15]);
colormap(cmpFreq);
cbar = colorbar();
ylabel(cbar, 'mean frequency (Hz)');
xlabel('u_1 (a.u.)');
ylabel('u_2 (a.u.)');
zlabel('u_3 (a.u.)');
ConfAxis('fontSize', 14);
view(viewport);

% Scatter embedding colored by limb position
limbList = {'L1','L2','L3','R1','R2','R3'};
for ind = 1:6
    MakeFigure;
    scattern(limbUMAP, 10, squeeze(L(floor(size(L,1)/2)+1, ind, :)), 'filled');
    axis('equal');
    cbar = colorbar();
    caxis([-1 1]);
    colormap(cmpLimb);
    ylabel(cbar, sprintf('%s(0) (a.u.)', limbList{ind}));
    cbar.Ticks = [-1 0 1];
    xlabel('u_1 (a.u.)');
    ylabel('u_2 (a.u.)');
    zlabel('u_3 (a.u.)');
    ConfAxis('fontSize', 14);
    view(viewport);
end

% Scatter embedding colored by the total number of feet down
ndown = squeeze(sum(D(floor(size(D,1)/2)+1,:,:),2));
MakeFigure;
scattern(limbUMAP, 10, ndown, 'filled');
axis('equal');
caxis([0 6]);
colormap(cmpDown);
cbar = colorbar();
ylabel(cbar, 'total feet down');
caxis([0 6]);
cbar.Ticks = (0+(6-1)/(2*6)):((7-1)/7):6;
cbar.TickLabels = 0:6;
xlabel('u_1 (a.u.)');
ylabel('u_2 (a.u.)');
zlabel('u_3 (a.u.)');
ConfAxis('fontSize', 14);
view(viewport);