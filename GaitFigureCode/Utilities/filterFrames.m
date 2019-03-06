function [ data ] = filterFrames( data, lowThresh, highThresh )
% This function filters out frames in the dataset that fall outside the
% desired linear velocity bounds

%% Set thresholds (mm per second)

if ~exist('lowThresh','var' ) || isempty(lowThresh)
    lowThresh = 0.5;
end

if ~exist('highThresh','var') || isempty(highThresh)
    highThresh = 70; % Inf in a previous version
end

fprintf('\nFiltering data by linear velocity.\n');
fprintf('Low threshold set at %f millimeters per second.\n',lowThresh);
fprintf('High threshold set at %f millimeters per second.\n',highThresh);

%% Filter the data

% Filter out the low velocity frames
data(data.linVel_mmPerSec <= lowThresh, :) = [];

% Filter out the high velocity frames
data(data.linVel_mmPerSec >= highThresh,:) = [];

end

