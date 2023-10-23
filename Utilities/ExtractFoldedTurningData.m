function [ ID, V, L, D, CM, vL, Phi, dPhi ] = ExtractFoldedTurningData( newData, winlen, excludeIncompleteTraces, smoothing )
%UNTITLED16 Summary of this function goes here
%   Detailed explanation goes here

%% Set parameters

minDist = 10; % frames
vR_threshold = deg2rad(0); % deg/s
vF_threshold = 1; % mm/s

fps = 150;
mmPerPix = 0.043;

if nargin < 3 || isempty(excludeIncompleteTraces)
    excludeIncompleteTraces = true;
end

% Select whether to smooth extracted data. Turns are always detected based
% on smoothed data
if nargin < 4 || isempty(smoothing)
    smoothing = true;
end

extractPhases = (nargout > 5);
min_prominence = 1;

%% Define variables to extract

idVarList = {'uniqueFlyTrajID','Frame', 'videoID'};

limbList = {'L1','L2','L3','R1','R2','R3'};

phaseVarList = strcat('InstantaneousPhase_', limbList, 'y');
frequencyVarList = strcat('InstantaneousFrequency_', limbList, 'y');

centroidVarList = {'xCOM','yCOM'};
velVarList = {'angVel_radPerSec','forwardSpeed_mmPerSec','translationalSpeed_mmPerSec', 'Orient_Rad_FullRotation'};

% downVarList = strcat(limbList, '_DOWN_NEW');
 downVarList = strcat(limbList, '_down_cam');

limbVarListX = strcat(limbList, '_xPlot_mm');
limbVarListY = strcat(limbList, '_yPlot_mm');
limbVelVarListX = strcat(limbList, '_xVel_Plot_mmPerSec');
limbVelVarListY = strcat(limbList, '_yVel_Plot_mmPerSec');

if smoothing
    frequencyVarList = strcat('smooth_', frequencyVarList);
    phaseVarList = strcat('smooth_', phaseVarList);
    centroidVarList = strcat('smooth_', centroidVarList);
    velVarList = strcat('smooth_', velVarList);
end

%% Detect turns
tic;

% Extract the needed data
vr = newData.smooth_angVel_radPerSec;
g = newData.uniqueFlyTrajID;

trigger_R = newData.yawExtremum > 0;
trigger_L = newData.yawExtremum < 0;


% Ensure that turns with less than the threshold yaw speed are excluded
trigger_R = trigger_R & (vr > vR_threshold);
trigger_L = trigger_L & (vr < -vR_threshold);

% Exclude turns with walking speed less than the desired walking speed
trigger_R = trigger_R & (newData.smooth_forwardSpeed_mmPerSec > vF_threshold);
trigger_L = trigger_L & (newData.smooth_forwardSpeed_mmPerSec > vF_threshold);

fprintf('Detected turns in %f seconds.\n', toc);


%% Extract data

% Define an anonymous function to alias the appropriate slicing function
if excludeIncompleteTraces
    slicefun = @(l, t) getSlicesSimple(newData, l, winlen, t);
else
    slicefun = @(l, t) getSlices(newData, l, winlen, t, true);
end

tic;

% Identifying information
[ ID_R ] = slicefun( idVarList, trigger_R);
[ ID_L ] = slicefun( idVarList, trigger_L);
ID_R = cat(3, ID_R{:});
ID_L = cat(3, ID_L{:});

% Centroid velocities
[ V_R ] = slicefun( velVarList, trigger_R);
[ V_L ] = slicefun( velVarList, trigger_L);
V_R = cat(3, V_R{:});
V_L = cat(3, V_L{:});

% Limb x-coordinates
[ LX_R ] = slicefun( limbVarListX, trigger_R);
[ LX_L ] = slicefun( limbVarListX, trigger_L);
LX_R = cat(3, LX_R{:});
LX_L = cat(3, LX_L{:});

% Limb y-coordinates
[ LY_R ] = slicefun( limbVarListY, trigger_R);
[ LY_L ] = slicefun( limbVarListY, trigger_L);
LY_R = cat(3, LY_R{:});
LY_L = cat(3, LY_L{:});

% Limb up/down
[ D_R ] = slicefun( downVarList, trigger_R);
[ D_L ] = slicefun( downVarList, trigger_L);
D_R = cat(3, D_R{:});
D_L = cat(3, D_L{:});

% COM positions
[ CM_R ] = slicefun( centroidVarList, trigger_R);
[ CM_L ] = slicefun( centroidVarList, trigger_L);
CM_R = cat(3, CM_R{:});
CM_L = cat(3, CM_L{:});

if extractPhases
    % Limb x-velocities
    [ vLX_R ] = slicefun( limbVelVarListX, trigger_R);
    [ vLX_L ] = slicefun( limbVelVarListX, trigger_L);
    vLX_R = cat(3, vLX_R{:});
    vLX_L = cat(3, vLX_L{:});

    % Limb y-velocities
    [ vLY_R ] = slicefun( limbVelVarListY, trigger_R);
    [ vLY_L ] = slicefun( limbVelVarListY, trigger_L);
    vLY_R = cat(3, vLY_R{:});
    vLY_L = cat(3, vLY_L{:});

    % Instantaneous phases
    [ Phi_R ] = slicefun( phaseVarList, trigger_R);
    [ Phi_L ] = slicefun( phaseVarList, trigger_L);
    Phi_R = cat(3, Phi_R{:});
    Phi_L = cat(3, Phi_L{:});

    % Instantaneous frequencies
    [ dPhi_R ] = slicefun( frequencyVarList, trigger_R);
    [ dPhi_L ] = slicefun( frequencyVarList, trigger_L);
    dPhi_R = cat(3, dPhi_R{:});
    dPhi_L = cat(3, dPhi_L{:});

end

fprintf('Extracted timeseries in %f seconds.\n', toc);

%% Fold left and right turns
tic;

% Fold timeseries
LY_L =  LY_L(:, [4,5,6, 1,2,3], :);
LX_L = -LX_L(:, [4,5,6, 1,2,3], :);
D_L  =   D_L(:, [4,5,6, 1,2,3], :);

% Flip the axes such that +x is rightward and +y is forward 
LX_L = -LX_L;
LX_R = -LX_R;
LY_L = -LY_L;
LY_R = -LY_R;

% Flip the signs of the angular and perpendicular velocities such that the
% angular velocity is positive (R) at the time of the detected turn
V_L(:,1,:) = (-1).*V_L(:,1,:);
V_L(:,3,:) = (-1).*V_L(:,3,:);

% Combine X and Y limb coordinate arrays
L_L = cat(2, LX_L, LY_L);
L_R = cat(2, LX_R, LY_R);

% Combine left and right limb coordinate and up/down arrays
L = cat(3, L_R,  L_L);
D = cat(3, D_R, D_L);

% Append an indicator to the ID array indicating whether a turn was
% originally to the right (1) or the left (0)
ID_R = cat(2, ID_R, ones(size(ID_R,1), 1, size(ID_R, 3)));
ID_L = cat(2, ID_L, zeros(size(ID_L,1), 1, size(ID_L, 3)));
ID = cat(3, ID_R, ID_L);
ID = ID(:,[1,2,4,3],:);
V  = cat(3, V_R, V_L);
CM = cat(3, CM_R, CM_L);

if extractPhases
    % Fold timeseries
    vLY_L =  vLY_L(:, [4,5,6, 1,2,3], :);
    vLX_L = -vLX_L(:, [4,5,6, 1,2,3], :);
    dPhi_L = dPhi_L(:, [4,5,6, 1,2,3], :);
    Phi_L  = Phi_L(:,  [4,5,6, 1,2,3], :);

    % Flip the y-axis such that +y is forward
    vLY_L = -vLY_L;
    vLY_R = -vLY_R;

    % Combine X and Y limb velocity arrays
    vL_L = cat(2, vLX_L, vLY_L);
    vL_R = cat(2, vLX_R, vLY_R);

    % Fold
    vL   = cat(3, vL_R,   vL_L);
    dPhi = cat(3, dPhi_R, dPhi_L);
    Phi  = cat(3, Phi_R,  Phi_L);
    CM   = cat(3, CM_R,   CM_L);

    % Convert units
    dPhi = dPhi./(2*pi).*fps;
    Phi = Phi./(2*pi);
end
clearvars *_L *_R;

% Convert units
V(:,1,:) = rad2deg(V(:,1,:));

fprintf('Folded timeseries in %f seconds.\n', toc);

%% Ensure that all incomplete traces are excluded

if excludeIncompleteTraces
    % Determine which variables we have to check for empty values
    if extractPhases
        notNan = ~squeeze(any(any(isnan(cat(2, ID,V,L,D,CM,vL,Phi,dPhi)),2),1));
    else
        notNan = ~squeeze(any(any(isnan(cat(2, ID,V,L,D,CM)),2),1));
    end

    % Remove incomplete timeseries
    ID = ID(:,:,notNan);
    V = V(:,:,notNan);
    L = L(:,:,notNan);
    D = D(:,:,notNan);
    CM = CM(:,:,notNan);

    if extractPhases
        vL = vL(:,:,notNan);
        Phi = Phi(:,:,notNan);
        dPhi = dPhi(:,:,notNan);
    end
end

end
