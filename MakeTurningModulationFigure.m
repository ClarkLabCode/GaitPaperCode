%% MAKETURNINGMODULATIONFIGURE.m
% This script generates the turning modulation figures from 
% "The manifold structure of limb coordination in walking Drosophila"

% Set data paths
SetDataPaths;

% Set minimum forward speed for yaw extremum detection
vfMin = 0.5;

%% Load data from file and prepare data

load(wtdatapath, 'newData');
[ newData ] = PrepareTurningData(newData, vfMin);

%% Run the analyses

% Set range of forward velocities on which to condition
vfInterval = [15, 20];

% Swing & stance durations conditioned on centroid kinematics
PlotSwingStanceStatisticsConditionedOnCentroidKinematics( newData, vfInterval );

% Step length & direction conditioned on centroid kinematics
PlotAverageStepLengthAndDirection( newData, vfInterval );