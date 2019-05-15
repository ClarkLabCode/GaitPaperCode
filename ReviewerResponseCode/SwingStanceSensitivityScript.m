% SwingStanceSensitivityScript.m: Sensitivity Analysis for Swing Stance Parameters
% This script performs a quick sensitivity analysis on the swing/stance
% determination and recreates the panels in Figure 2 to show that
% qualitatively results don't change.

% NOTE: Need to update these variables based on the calculation of swing and stance
% varList = {'L1_down_cam', 'L2_down_cam', 'L3_down_cam','R1_down_cam', 'R2_down_cam', 'R3_down_cam'};

%% Load the dataset
% Load the wild type data
tic;
path = '/Users/bdeangelis/Desktop/Datasets/_CURRENT_DATASETS/IsoD1-Masked-1000PCs/';
file = '20181025_20180530-20180614_IsoD1_Glass_MaskedModel_1000PCs_amplitude_phase_down_downcam_Steps(_down_cam).mat';
data_path = strcat(path,file);
load(data_path, 'newData');
fprintf('Loaded the desired dataset: %f seconds.\n',toc);


%% Set parameters of the current run
% NOTE: Default value is 20
% threshold = {10, 20, 30};
threshold = 30;

%% Compute the dataset for a given forward velocity threshold in camera frame
newData = appendMovingAverageUpDownLogical_SensitivityUtility( newData, threshold );

%% Define some defaults for the analysis

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
cmap_density = copper(256);

% Colormap: Number of feet down
cmap_numfeet = jet(7) .* 0.8;

%% Compute the corresponding figures

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

%% Save the figures

