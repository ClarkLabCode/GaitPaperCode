%% ReviewResponseFigureScript.m: This script generates all the reviewer response figures

%% Load the data

% Define the parent path for the location of all data files. Keep all 
% datasets in this directory.
% parentPath = '';
parentPath = '/Users/bdeangelis/Desktop/Datasets/GaitPaperDatasets/';

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

%% Figure 1: Swing-Stance Sensitivity

% Compute the updated swing-stance logicals
threshold = 10; % 10 pixels/frame --> 50% decrease from paper
[ newData_UpDownChange ] = addSmoothedCameraFrameUpDown_variableThreshold( newData, threshold );

% Plot the figure components

% % A,B: Example Gaits
% PlotCanonicalGaits();

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
PlotStanceConfigurationFractions( newData_UpDownChange, analysisType, binEdges, corder, faceAlpha, minorLineWidth, smoothing, removeDiscontinuities, limbVarType, suppressExtraPlots );

% E: Stance Duration Plots
vfBinEdges = [10.1983, 19.0112]; % 33.3% and 66.6% Percentiles in newData.forwardSpeed_mmPerSec > 0.5 mm/s
corder = cmap_numfeet;
limbVarType = 'Camera';
PlotStanceDurationDwellTimesVsSpeed( newData_UpDownChange, vfBinEdges, corder, limbVarType );

% F and Figure 2, figure supp. 2: Forward Speed vs. Phase number of feet down
nbins = 20;
vfBinEdges = [10.1983, 19.0112]; % 33.3% and 66.6% Percentiles in newData.forwardSpeed_mmPerSec > 0.5 mm/s
cmap = cmap_density;
corder = cmap_numfeet;
limbVarType = 'Camera';
PlotFeetDownVsLimbPhase( newData_UpDownChange, nbins, vfBinEdges, cmap, corder, limbVarType);

% Clear the threshold changed dataset
clear newData_UpDownChange

%% Figure 3: Example Trajectories with Mendes 2013 metrics

% Calculate and plot the example trajectories with the desired metrics added
calcGaitIndex( newData );