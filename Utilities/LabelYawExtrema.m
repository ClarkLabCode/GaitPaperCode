function [ yawExtremum ] = LabelYawExtrema(newData, turn_method, vF_threshold, smoothing, min_prominence, minDist)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here

%% Parse arguments
if nargin < 3
    vF_threshold = 0;
end

if nargin < 4
    smoothing = true;
end

if nargin < 5
    min_prominence = 1;
end

if nargin < 6
    minDist = 0; % frames
end

%% Locate yaw extrema

tic;

% Extract the needed data
g = findgroups(newData.uniqueFlyTrajID);
if smoothing
    vf = newData.smooth_forwardSpeed_mmPerSec;
    vr = newData.smooth_angVel_radPerSec;
else
    vf = newData.forwardSpeed_mmPerSec;
    vr = newData.angVel_radPerSec;
end

% Locate maxima in the rotational velocity
switch turn_method
    case 'findpeaks'
        [ trigger_R ] = findLocalExtrema('maxima', vr, g, minDist, 0, min_prominence);
        [ trigger_L ] = findLocalExtrema('minima', vr, g, minDist, 0, min_prominence);
    case 'katsov'
        [ trigger_R ] = findLocalExtremaKatsov('maxima', vr, g, minDist, 0);
        [ trigger_L ] = findLocalExtremaKatsov('minima', vr, g, minDist, 0);
    case 'rdp'
        trigger_R = (newData.RDPangle > 0);
        trigger_L = (newData.RDPangle < 0);
    otherwise
        error('Invalid extremum detection method: %s', turn_method);
end

% Exclude turns with walking speed less than the desired walking speed
trigger_R = trigger_R & (vf > vF_threshold);
trigger_L = trigger_L & (vf > vF_threshold);

yawExtremum = trigger_R - trigger_L;

fprintf('Labeled yaw extrema in %f seconds.\n', toc);

end

