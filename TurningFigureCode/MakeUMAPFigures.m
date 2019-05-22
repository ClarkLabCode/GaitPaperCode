%% MAKEUMAPFIGURES.m
% This script generates the UMAP figures from 
% "The manifold structure of limb coordination in walking Drosophila"

% Set data paths
SetTurningDataPaths;

% Set minimum forward speed for yaw extremum detection
vfMin = 0.5;

% Define cache file paths
umapdatapath = fullfile(basepath, 'umap_data.h5');
umapresultpath = fullfile(basepath, 'umap_result.h5');

%% Load data from file and prepare data

load(wtdatapath, 'newData');
[ newData ] = PrepareTurningData(newData, vfMin);

%% Set UMAP sampling + embedding parameters

% Set whether to use smoothed or raw centroid kinematic data
smoothing = true;

% Set the half-window length
winlen = 15;

% Sample size
ns = 1e5;

% Set the sampling method (see contents of DefineRandomSampleOfTrajectorySegments)
samplingMethod = 'all';

% If turns are explicitly included, select whether to fold left and right
foldTurns = true;

% Sampling parameters
vfMin = 0.5;
vfMax = 40;
vfBin = 1;

% Set the target dimensionalities
centroid_components = 2;
limb_components = 3;

% Set whether to mean-subtract trajectories
meanSubtractTrajectories = true;

%% Prepare UMAP input data

% This script prepares everything 
PrepareDataForUMAP;

%% Make PCA visualizations

MakeLimbDataPCAVisualizations;

%% Write cache paths to a text file so that the Python code can find them
% This is a terrible bodge, and assumes that this code is run with the base
% TurningPaperCode directory as the working directory.

tic;
cfgpath = fullfile('AnalysisUtilities','config.txt');
fid = fopen(cfgpath, 'w');
fprintf(fid, strrep([umapdatapath,',',umapresultpath], '\','\\'));
fclose(fid);
fprintf('Wrote UMAP cache file paths to configuration file in %f seconds.\n', toc);

%% Write data to file to run UMAP from Python
% This block creates an HDF5 cache file containing the data expected by the
% Python UMAP script

tic;
% If a temporary data file exists, delete it
if exist(umapdatapath, 'file') == 2
    delete(umapdatapath);
end

% Write data arrays to the file
h5create(umapdatapath, '/x_v', size(X_V));
h5write(umapdatapath, '/x_v', X_V);
h5create(umapdatapath, '/x_l', size(X_L));
h5write(umapdatapath, '/x_l', X_L);

% Write embedding parameters to the file
h5create(umapdatapath, '/centroid_components', 1);
h5write(umapdatapath, '/centroid_components', centroid_components);
h5create(umapdatapath, '/limb_components', 1);
h5write(umapdatapath, '/limb_components', limb_components);

fprintf('Wrote data to HDF5 file in %f seconds.\nPath: %s\n', toc, umapdatapath);

%% Run the UMAP script in Python
% One could call Python from inside MATLAB, but it's generally easier to
% just run the Python script manually.

%% Load UMAP results into MATLAB from file
% The Python UMAP script generates an HDF5 cache file containing the
% embeddings. 

tic;
centroidUMAP = h5read(umapresultpath, '/embedding_v');
limbUMAP = h5read(umapresultpath, '/embedding_l');

centroidUMAP = centroidUMAP';
limbUMAP = limbUMAP';

fprintf('Read result from HDF5 file in %f seconds.\nPath: %s\n', toc, umapresultpath); 

%% Make UMAP scatter plots
% It's probably advisible to run this script block-by-block, as the size of
% the plots in memory will be quite large.

MakeCentroidAndLimbUMAPScatterPlots;

%% Make cylindrical coordinate plots
MakeUMAPCylindricalCoordinatePlots;
