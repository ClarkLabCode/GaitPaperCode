function newData = appendMovingAverageUpDownLogical_SensitivityUtility( newData, threshold )
% Quick utility for appending the updated SWING/STANCE logical and saving a 
% copy in the same directory as the original files
    
% Append the UP/DOWN logical based on moving average limb velocities
[ newData ] = addMovingAverageUpDown( newData );

% Remove single frame discontinutities
[ newData ] = RemoveSingleFrameUpDownDiscontinuities_local( newData, '_down' );

% Append the UP/DOWN logical based on smoothed Camera Frame Movement
[ newData ] = addSmoothedCameraFrameUpDown_variableThreshold( newData, threshold );

% Remove single frame discontinuities from camera up/down
[ newData ] = RemoveSingleFrameUpDownDiscontinuities_local( newData, '_down_cam' );
    
end

function [ newData ] = RemoveSingleFrameUpDownDiscontinuities_local( newData, suffix )
% Remove single frame discontinuities from the datasets

% Define variable lists
limbList = {'L1','L2','L3','R1','R2','R3'};
% downVarList = strcat(limbList, '_down');
downVarList = strcat(limbList, suffix );

% Extract the needed data
id = newData.uniqueFlyTrajID;
X = newData{:, downVarList};

% Find single-frame discontinutities
idx = (X ~= circshift(X,-1)) & (circshift(X,-1) == circshift(X,1)) & (id == circshift(id,1)) & (id == circshift(id,-1));

% Remove those discontinutities
X(idx &  circshift(X,-1)) = true;
X(idx & ~circshift(X,-1)) = false;

% Append the modified data
newData{:, downVarList} = X;

end