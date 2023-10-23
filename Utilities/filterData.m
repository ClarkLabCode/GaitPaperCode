function [ newData ] = filterData( newData, so, sl, gl, gs, ms, use_moving_average )
%FILTERDATA This function applies all needed smoothing and filtering to
%various variables in the data table and appends the results to the table.

%% Check inputs

if ~exist('so','var') || isempty(so)
    so = 3;
end

if ~exist('sl','var') || isempty(sl)
    sl = 15;
end

if ~exist('gl','var') || isempty(gl)
    gl = 15;
end

if ~exist('gs','var') || isempty(gs)
    gs = 3;
end

if ~exist('ms','var') || isempty(ms)
    ms = 5;
end

% Select whether to use a Gaussian kernel or moving average filter to
% smooth centroid data
if ~exist('use_moving_average','var') || isempty(use_moving_average)
    use_moving_average = false;
end

%% Smooth phases and approximate frequencies

limbList = {'L1','L2','L3','R1','R2','R3'};
limbList = [strcat(limbList, 'x'), strcat(limbList, 'y')];
pvl = strcat('InstantaneousPhase_', limbList);
fvl = strcat('InstantaneousFrequency_', limbList);
avl = strcat('InstantaneousAmplitude_', limbList);


spvl = strcat('smooth_', pvl);
sfvl = strcat('smooth_', fvl);
savl = strcat('smooth_', avl);

% Approximate instantaneous frequencies
tic;
[ newData{:, sfvl} ] = savitzkyGolayFilter( newData, pvl, 1, so, sl );
fprintf('Approximated instantaneous frequencies in %f seconds.\n', toc);

% Smooth instantaneous phases
tic;
[ newData{:, spvl} ] = savitzkyGolayFilter( newData, pvl, 0, so, sl );
fprintf('Smoothed instantaneous phases in %f seconds.\n', toc);

% Smooth instantaneous amplitudes
tic;
[ newData{:, savl} ] = SmoothWithMovingAverageFilter(newData, avl, ms);
fprintf('Smoothed instantaneous amplitudes in %f seconds.\n', toc);

%% Smooth centroid dynamics

vl = {'xCOM','yCOM','Orient_Rad_FullRotation',...
    'angVel_radPerSec','linVel_mmPerSec',...
    'forwardSpeed_mmPerSec','translationalSpeed_mmPerSec', 'AngleBody2Displacement'};
svl = strcat('smooth_', vl);

% Smooth centroid traces
if use_moving_average
    tic;
    [ newData{:, svl} ] = SmoothWithMovingAverageFilter(newData, vl, ms);
    fprintf('Smoothed centroid velocities using a moving average filter in %f seconds.\n', toc);
else
    tic;
    [ newData{:, svl}, ~ ] = smoothWithGaussianKernel( newData, vl, gl, 'sigma', gs, 'method', 'filtfilt');
    fprintf('Smoothed centroid velocities using a Gaussian filter in %f seconds.\n', toc);
end


end

