function [ newData ] = addSmoothedCameraFrameUpDown_variableThreshold( newData, threshold )
% This function generates a test alternative UP/DOWN logical

% Define the variables
limbList = {'L1','L2','L3','R1','R2','R3'};
limbYVarList = strcat(limbList, '_yVel_Cam_mmPerSec');
limbXVarList = strcat(limbList, '_xVel_Cam_mmPerSec');
limbVarList = [limbYVarList, limbXVarList];
downList = strcat(limbList, '_down_cam');

% Calculate the smoothed version of the limb velocities in the camera frame
[ vels ] = local_SmoothWithMovingAverageFilter(newData, limbVarList, 5); % Current

% Calculate the magnitude of the velocity vector at each time-point
Y = vels(:,1:6);
X = vels(:,7:12);

D = hypot(X,Y);

% NOTE: This cutoff is somewhat arbitrary, but by visual inspection generates reasonable results
% NOTE: Default for gait paper was 20 --> 3.1 pix/frame
down = D < threshold; % (20 mm/s --> 3.1  pix/frame) 

% Place the data in the array of interest
newData{:,downList} = down;

end

function [ Y ] = local_SmoothWithMovingAverageFilter(newData, varList, span)
% This function generates a smoothed version of a selected set of variables

% Find the groups corresponding to trajectories
[G, ~] = findgroups(newData.uniqueFlyTrajID);

% Format data for splitapply
X = newData{:, varList};

% Smooth each trajectory
[ Y ] = splitapply(@(x) {local_moving_average(x, span)}, X, G);

% Combine data for output
Y = cat(1, Y{:});

end

function [ out ] = local_moving_average(y, span)
% Apply moving average filter along the first dimension of multiple
% dimensional data. Adapted from the 'moving' subfunction of MATLAB's
% smooth.m from the Curve Fitting Toolbox. This can only handle cases in
% which NaN values are concentrated at the beginning and end of
% trajectories

% Allocate
out = nan(size(y));

% Check to make sure the timeseries is long enough
if size(y,1) > span
    
    % Find nan rows
    ynan = any(isnan(y),2);
    
    % Check to make sure that there are some non-nan rows
    if nnz(~ynan) > span
        
        % Remove NaN rows
        y = y(~ynan,:);
        
        % Actual filtering logic, adapted to operate along the first
        % dimension of the array...
        span = floor(span);
        n = size(y,1);
        span = min(span,n);
        width = span-1+mod(span,2); % force it to be odd
        c = filter(ones(width,1)/width,1,y,[],1);
        cbegin = cumsum(y(1:width-2,:),1);
        cbegin = cbegin(1:2:end, :)./(1:2:(width-2))';
        cend = cumsum(y(n:-1:n-width+3,:), 1);
        cend = cend(end:-2:1,:)./(width-2:-2:1)';
        c = [cbegin; c(width:end,:); cend];
        
        % Store the filtered data
        out(~ynan,:) = c;
    end
end
end