% GaitPaperFiguresScript.m: Script For organizing and generating figures 
% for "The manifold structure of limb coordination in walking Drosophila"

% Define the parent path for the location of all data files. Keep all 
% datasets in this directory.
parentPath = '';

%% Load the relevant datasets

% Load the wild type data
tic;
path = parentPath;
file = '20181025_20180530-20180614_IsoD1_Glass_MaskedModel_1000PCs_amplitude_phase_down_downcam_Steps(_down_cam).mat';
data_path = strcat(path,file);
load(data_path, 'newData');
fprintf('Loaded the desired dataset: %f seconds.\n',toc);

%% Set some parameters for the plots

% Phase calculation type (Hilbert or Haar)
phaseType = 'Hilbert';
if strcmp(phaseType, 'Haar')
    [ newData ] = ComputeHaarWaveletInstantaneousPhase( newData );
end

% Linewidths
figureLineWidth = 1;

% Smoothing
figureSmoothing = true;

% Colormap: Body Velocity Components
cmap_bodyvel = lines(6);

% Colormap: Six Limbs
cmap_sixlegs = linspecer(6);

% Colormap: Limb Relative Phase
colorMatrix = [1,0,0;0,0,1];
x1 = colorMatrix;
x2 = colorMatrix .* 0.75;
x3 = colorMatrix .* 0.5;
rbmap(1,:) = x1(1,:);
rbmap(2,:) = x2(1,:);
rbmap(3,:) = x3(1,:);
rbmap(4,:) = x1(2,:);
rbmap(5,:) = x2(2,:);
rbmap(6,:) = x3(2,:);
cmap_limbrel = rbmap;

% Colormap: Density Plots
[ ~, cmap_density, ~ ] = MakeDensityColormaps(); % White to Red colormap
% cmap_density = copper(256);

% Colormap: Number of feet down
cmap_numfeet = jet(7) .* 0.8;

%% Preprocess the dataset

% Calculate smoothed versions of the phase and frequency variables and the centroid velocity components
[ newData ] = filterData( newData, [], [], [], [] );

% Remove frames in which the fly stops
[ newData_filtered ] = filterFrames( newData, [], [] );

fprintf('Preprocessed the dataset: %f seconds.\n',toc);

%% Figure 1: Measurements of body and limb kinematics in freely-walking Drosophila.

% A: NO PLOT -- Diagram of rig

% B: NOT GENERATED IN THIS SCRIPT -- Frame from example video

% C: PDFs of body variables with bootstrapped errorbars
% Set the parameters
num_bootstraps = 1000;
speedBinEdges = linspace(0, 35, 70);
vrBinEdges = linspace(-300, 300, 100);
vfBinEdges = linspace(-5, 35, 100);
vpBinEdges = linspace(-10, 10, 100);
corder = cmap_bodyvel;
cmap = cmap_density;
linewidth = figureLineWidth;
smoothing = figureSmoothing;
PlotBootstrapVelocityDistributions( newData, num_bootstraps, speedBinEdges, vrBinEdges, vfBinEdges, vpBinEdges, corder, cmap, linewidth, smoothing);

% D: Joint densities of body velocity components
% Set the parameters
speedBinEdges = linspace(0, 50, 200);
vrBinEdges = linspace(-300, 300, 100);
vfBinEdges = linspace(-5, 35, 100);
vpBinEdges = linspace(-10, 10, 100);
cmap = cmap_density;
logscale = true;
linewidth = figureLineWidth;
smoothing = figureSmoothing;
PlotCentroidVelocityJointDistributions( newData, speedBinEdges, vfBinEdges, vpBinEdges, vrBinEdges, cmap, logscale, linewidth, smoothing);

% E: Limb Measurements
num_bootstraps = 1000;
vfedges = [5:1:35];
corder_all = cmap_sixlegs;
[ f, gof, output ] = calcPowerLaw( newData );
[ swing_counts, stance_counts ] = PlotStepStatisticDistributionsBootstraps( newData, num_bootstraps, vfedges, corder_all, f, gof );

%% Figure 2: Drosophila uses a two-cycle gait across all walking speeds.

% A,B: Example Gaits
PlotCanonicalGaits();

% C,D: Stance Configuration Fractions
analysisType = 'forward';
binEdges = 0:0.5:35;
corder = cmap_numfeet;
faceAlpha = 0.6;
minorLineWidth = figureLineWidth;
smoothing = false;
removeDiscontinuities = true;
limbVarType = 'Camera';
suppressExtraPlots = true;
PlotStanceConfigurationFractions( newData, analysisType, binEdges, corder, faceAlpha, minorLineWidth, smoothing, removeDiscontinuities, limbVarType, suppressExtraPlots );

% E: Stance Duration Plots
vfBinEdges = [10.1983, 19.0112]; % 33.3% and 66.6% Percentiles in newData.forwardSpeed_mmPerSec > 0.5 mm/s
corder = cmap_numfeet;
limbVarType = 'Camera';
PlotStanceDurationDwellTimesVsSpeed( newData, vfBinEdges, corder, limbVarType );

% F and Figure 2, figure supp. 2: Forward Speed vs. Phase number of feet down
nbins = 20;
vfBinEdges = [10.1983, 19.0112]; % 33.3% and 66.6% Percentiles in newData.forwardSpeed_mmPerSec > 0.5 mm/s
cmap = cmap_density;
corder = cmap_numfeet;
limbVarType = 'Camera';
PlotFeetDownVsLimbPhase( newData, nbins, vfBinEdges, cmap, corder, limbVarType);

%% Figure 3: Relative phase measurements reveal a continuum across all walking speeds with contralateral limbs in antiphase.

% A: NO PLOT -- Diagram of relative phases

% B: Polar plot of limb relative phases
nbins = 180;
corder_pairwise = cmap_limbrel;
corder_all = cmap_sixlegs;
rlims = [];
vfBinEdges = [0, 10.1983, 19.0112]; % 33.3% and 66.6% Percentiles in newData.forwardSpeed_mmPerSec > 0.5 mm/s
smoothing = figureSmoothing;
suppressExtraPlots = true;
PlotRelativeInstantaneousPhasePolarHistograms( newData_filtered, nbins, corder_pairwise, corder_all, rlims, vfBinEdges, smoothing, suppressExtraPlots );
PlotRelativeInstantaneousPhasePolarPDF( newData_filtered, nbins, corder_pairwise, corder_all, rlims, vfBinEdges, smoothing, suppressExtraPlots );

% C,D: Forward Speed vs. Circular Mean and Angular Deviation of Relative Phases
vfBinEdges = 5:1:35;
vrBinEdges = 0:10:350;
nbins = 100;
cmap = cmap_density;
logscale = true;
corder = cmap_limbrel;
linewidth = figureLineWidth;
smoothing = figureSmoothing;
showDistributions = false;
suppressExtraPlots = true;
num_bootstraps = 100;
PlotRelativeInstantaneousPhaseVsVelocityComponentsBootstrapped( newData, vfBinEdges, vrBinEdges, nbins, cmap, logscale, corder, linewidth, smoothing, showDistributions, suppressExtraPlots, num_bootstraps);

% E: Relative Phase Joint Distribution
nbins = 100;
cmap = cmap_density;
logscale = true;
clims = [];
vfBinEdges = [0, 10.1983, 19.0112]; % 33.3% and 66.6% Percentiles in newData.forwardSpeed_mmPerSec > 0.5 mm/s
smoothing = figureSmoothing;
pointColor = 'k';
textColor = 'k';
suppressExtraPlots = false;
PlotRelativeInstantaneousPhaseJointDistributions_nonSymmetrized( newData_filtered, nbins, cmap, logscale, clims, vfBinEdges, smoothing, pointColor, textColor, suppressExtraPlots );

% F: Example trajectory step plots
colors = cmap_sixlegs;
colors_all = cmap_limbrel;
PlotTripodExample
PlotTetrapodExample
PlotNonCanonicalExample

% G: Stance duration with example trajectories overlaid
logscale = false;
clims = [];
cmap = cmap_density;
plot_type = 'conditional';
PlotStanceDurationDensity( newData, logscale, clims, cmap, plot_type );

% H,I,J & Figure 3, figure supp. 2C
cmp = cmap_density;
corder = cmap_sixlegs(1:4,:);
PlotCanonicalGaitModeCoherences( newData, cmp, corder );

%% Supplemental Figure Components: Additional Measures of Limb Phase

% Figure 2, figure supp. 1: Analytic Signal Decomposition
path = parentPath;
file = 'Tripod-gait.mat';
data_path = strcat(path,file);
load(data_path, 'singleFly');
PlotHilbertTransformDecomposition();

% Figure 3, figure supp. 2A: Speed Partitions of Relative Phase Joint Probability Distributions
nbins = 100;
cmap = cmap_density;
logscale = true;
clims = [];
vfBinEdges = [0, 10.1983, 19.0112]; % 33.3% and 66.6% Percentiles in newData.forwardSpeed_mmPerSec > 0.5 mm/s
smoothing = figureSmoothing;
pointColor = 'k';
textColor = 'k';
suppressExtraPlots = false;
PlotRelativeInstantaneousPhaseJointDistributions_nonSymmetrized( newData_filtered, nbins, cmap, logscale, clims, vfBinEdges, smoothing, pointColor, textColor, suppressExtraPlots );

% Figure 3, figure supp. 2B: Variation in relative phase as function of limb phase
speed_edges = [5:5:35];
corder = flipud(cmap_numfeet);
smoothing = figureSmoothing;
PlotRelativePhaseVsPhaseMultipleSpeedsFrames( newData, speed_edges, corder, smoothing );

%% Figure 4: Dimensionality reduction reveals the manifold structure of limb coordination patterns.

% See UMAP code.

%% Figure 5: A single-parameter model with speed-independent coupling predicts continuum of inter-limb coordination patterns.

% Decide if you want to add noise to the model
addNoise = false;

% Clear the IsoD1 dataset
clear newData

% Load the model data
path = parentPath;
file = 'step_aspect_ratio_model_integrated_20181031-MatchExpDataFormat.mat';
data_path = strcat(path,file);
load(data_path, 'newData');
modelData = newData;
clear newData

% Add noise to limb and phase variables
if addNoise
    disp('Adding noise to model limb positions and phases.');
    limbVarList = {'L1','L2','L3','R1','R2','R3'};
    yPosList = strcat(limbVarList,'_yPlot');
    X = modelData{:,yPosList};
    SNR = 20; % Estimate of the white noise in our data
    X = awgn(X,SNR);
    modelData{:,yPosList} = X;
    [ modelData ] = computeAnalyticSignal( modelData );
else
    disp('Using noise-free model of limb positions and phases.');
end

% A,B: Model Diagram (One side first, Complete model second)
% No Code

% A,B: Example trajectories at different forward speeds
sample_length = 100;
PlotModelExampleTrajectories( modelData, sample_length );

% C: Stance Configuration Fractions
analysisType = 'forward';
binEdges = 5:.5:25;
corder = cmap_numfeet;
faceAlpha = 0.6;
minorLineWidth = figureLineWidth;
smoothing = false;
removeDiscontinuities = true;
limbVarType = 'Camera';
suppressExtraPlots = true;
PlotStanceConfigurationFractions( modelData, analysisType, binEdges, corder, faceAlpha, minorLineWidth, smoothing, removeDiscontinuities, limbVarType, suppressExtraPlots );

% D & Figure 5, Figure supp. 1: Two-cycle plots (Speed Bins and Density)
nbins = 20;
vfBinEdges = [13.8363, 20.1187]; % 33.3% and 66.6% Percentiles in modelData.forwardSpeed_mmPerSec
cmap = cmap_density;
corder = cmap_numfeet;
limbVarType = 'Camera';
PlotFeetDownVsLimbPhase( modelData, nbins, vfBinEdges, cmap, corder, limbVarType);

% E,F: Forward Speed vs. Circular Mean and Angular Deviation of Relative Phases
vfBinEdges = 5:2:25;
vrBinEdges = 0:10:350;
nbins = 100;
cmap = cmap_density;
logscale = true;
corder = cmap_limbrel;
linewidth = figureLineWidth;
smoothing = false;
showDistributions = false;
suppressExtraPlots = true;
PlotRelativeInstantaneousPhaseVersusCentroidVelocityComponents( modelData, vfBinEdges, vrBinEdges, nbins, cmap, logscale, corder, linewidth, smoothing, showDistributions, suppressExtraPlots);

% G: Model UMAP
% See UMAP code

% Additional Plots: Relative Phase Joint Distribution
nbins = 50;
cmap = cmap_density;
logscale = true;
clims = [];
vfBinEdges = [0, 13.8363, 20.1187]; % 33.3% and 66.6% Percentiles in modelData.forwardSpeed_mmPerSec
smoothing = false; % True in experimental data
pointColor = 'r';
textColor = 'w';
suppressExtraPlots = false;
PlotRelativeInstantaneousPhaseJointDistributions_nonSymmetrized( modelData, nbins, cmap, logscale, clims, vfBinEdges, smoothing, pointColor, textColor, suppressExtraPlots );

% Additional Plots: Relative Phase vs. Phase
% Model
speed_edges = [5:5:25];
corder = flipud(cmap_numfeet);
smoothing = false;
PlotRelativePhaseVsPhaseMultipleSpeedsFrames( modelData, speed_edges, corder, smoothing );

%% Figure 6: Experimental perturbations of walking speed modulate stance duration.

% Clear the Model dataset
clear modelData

% Load the moonwalker data
path = parentPath;
file = '20181025_20180501-20180511_Moonwalker-8ms_Glass_MaskedModel_1000PCs_amplitude_phase_down_downcam_Steps(_down_cam).mat';
data_path = strcat(path,file);
load(data_path, 'newData');
moonData = newData;
clear newData
% A,B: Moonwalker Plots
suppressPlots = true;
num_bootstraps = 1000;
[p_ctrl, h_ctrl, stats_ctrl, p_exp, h_exp, stats_exp, h_moon, p_moon, ks2stat_moon] = PlotMoonwalkerGaitFigures( moonData, figureLineWidth, suppressPlots, num_bootstraps );

% C,D: Select and plot moonwalker example
lineThickness = figureLineWidth;
cmap = cmap_sixlegs;
picked_ID = 42;
singleFly = PlotExampleTrajectory( moonData, 'Optogenetic', 'Moonwalker', lineThickness, cmap, picked_ID );

% Load the visual stimulus data
path = parentPath;
file = '20181025_VisualStimulus-isod1_50ms_random_dots_20180126-merged_down_downcam_Steps(_down_cam).mat';
data_path = strcat(path,file);
load(data_path, 'newData');
vizData = newData;
clear newData
% E,F: Visual Stimulus Plots
suppressPlots = true;
num_bootstraps = 1000;
[p_ctrl, h_ctrl, stats_ctrl, p_exp, h_exp, stats_exp, h_viz, p_viz, ks2stat_viz] = PlotVisualStimulusGaitFigures( vizData, figureLineWidth, suppressPlots, num_bootstraps );

% G,H: Select and plot visual example
lineThickness = figureLineWidth;
cmap = cmap_sixlegs;
picked_ID = 29;
singleFly = PlotExampleTrajectory( vizData, 'Visual', 'Visual', lineThickness, cmap, picked_ID );

%% Figure 1, figure supp. 1: Linear regression accurately identifies footfall positions.

% See feature extraction code

%% Figure 3, figure supp. 1: Synthetic canonical gaits differ in relative phase measurements from those of free-walking Drosophila.

% Load the data
path = parentPath;
file = '20181126_CanonicalGaits-DataTable.mat';
data_path = strcat(path,file);
load(data_path, 'newData');

% Partition the data into each of the individual gaits
tripodData = newData(newData.gaitType == 1,:);
leftTetraData = newData(newData.gaitType == 2,:);
rightTetraData = newData(newData.gaitType == 3,:);
waveData = newData(newData.gaitType == 4,:);

% A: Relative phase polar plots
nbins = 180;
corder_pairwise = cmap_limbrel;
corder_all = cmap_sixlegs;
rlims = [];
vfBinEdges = [0, 10.1983, 19.0112]; % Bounds for experiments. Dummy Variable here.
smoothing = false;
suppressExtraPlots = true;
% Tripod 
PlotRelativeInstantaneousPhasePolarHistograms( tripodData, nbins, corder_pairwise, corder_all, rlims, vfBinEdges, smoothing, suppressExtraPlots );
title('Tripod');
PlotRelativeInstantaneousPhasePolarPDF( tripodData, nbins, corder_pairwise, corder_all, rlims, vfBinEdges, smoothing, suppressExtraPlots );
title('Tripod');
% Left Tetrapod
PlotRelativeInstantaneousPhasePolarHistograms( leftTetraData, nbins, corder_pairwise, corder_all, rlims, vfBinEdges, smoothing, suppressExtraPlots );
title('Left Tetrapod');
PlotRelativeInstantaneousPhasePolarPDF( leftTetraData, nbins, corder_pairwise, corder_all, rlims, vfBinEdges, smoothing, suppressExtraPlots );
title('Left Tetrapod');
% Right Tetrapod
PlotRelativeInstantaneousPhasePolarHistograms( rightTetraData, nbins, corder_pairwise, corder_all, rlims, vfBinEdges, smoothing, suppressExtraPlots );
title('Right Tetrapod');
PlotRelativeInstantaneousPhasePolarPDF( rightTetraData, nbins, corder_pairwise, corder_all, rlims, vfBinEdges, smoothing, suppressExtraPlots );
title('Right Tetrapod');
% Wave
PlotRelativeInstantaneousPhasePolarHistograms( waveData, nbins, corder_pairwise, corder_all, rlims, vfBinEdges, smoothing, suppressExtraPlots );
title('Wave');
PlotRelativeInstantaneousPhasePolarPDF( waveData, nbins, corder_pairwise, corder_all, rlims, vfBinEdges, smoothing, suppressExtraPlots );
title('Wave');

% B: Relative phase joint distribution
nbins = 100;
cmap = cmap_density;
logscale = true;
clims = [];
vfBinEdges = [0, 10.1983, 19.0112]; % Bounds for experiments.
smoothing = false;
pointColor = 'k';
textColor = 'k';
suppressExtraPlots = true;
% Tripod 
PlotRelativeInstantaneousPhaseJointDistributions_nonSymmetrized( tripodData, nbins, cmap, logscale, clims, vfBinEdges, smoothing, pointColor, textColor, suppressExtraPlots );
% Left Tetrapod
PlotRelativeInstantaneousPhaseJointDistributions_nonSymmetrized( leftTetraData, nbins, cmap, logscale, clims, vfBinEdges, smoothing, pointColor, textColor, suppressExtraPlots );
% Right Tetrapod
PlotRelativeInstantaneousPhaseJointDistributions_nonSymmetrized( rightTetraData, nbins, cmap, logscale, clims, vfBinEdges, smoothing, pointColor, textColor, suppressExtraPlots );
% Wave
PlotRelativeInstantaneousPhaseJointDistributions_nonSymmetrized( waveData, nbins, cmap, logscale, clims, vfBinEdges, smoothing, pointColor, textColor, suppressExtraPlots );


%% Figure 4, figure supp. 1: Manifold structure of synthetic canonical gaits differ qualitatively from structure of free-walking Drosophila.

% A: UMAP Embedding - Equal Proportions
% See UMAP code

% B: UMAP Embedding - Low frequency of non-tripod
% See UMAP code

