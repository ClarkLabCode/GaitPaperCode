function [ newData ] = computeRamerDouglasPeuckerTurns( newData, epsilon, minFrames, smoothing )
%A function to compute turns in trajectories using the
%Ramer-Douglas-Peucker algorithm to simplify the trajectory.

% Note that a counterclockwise angle is NEGATIVE (which corresponds to the
% fly turning LEFT), and a clockwise angle is POSITIVE (which corresponds
% to the fly turning RIGHT). The sign convention thus matches that of the
% angular velocities derived from the orientation angle returned by OpenCV
% FitEllipse.

%% Validate inputs
if ~exist('newData','var') || isempty(newData) || ~isa(newData, 'table')
    error('Input ''newData'' must be a table.');
end
if ~exist('epsilon','var') || isempty(epsilon) || ~isa(epsilon, 'numeric') || (numel(epsilon) ~= 1)
    error('Input ''epsilon'' must be a scalar numeric threshold for the Raymer-Douglas-Peucker algorithm.');
end
if ~exist('minFrames','var') || isempty(minFrames) || ~isa(minFrames, 'numeric') || (numel(minFrames) ~= 1)
    error('Input ''minFrames'' must be a scalar numeric threshold on the minimum trajectory length to consider.');
end
if nargin < 4
    smoothing = false;
end

%% Simplify trajectories using the RDP algorithm

% Find the groups
[G, ~] = findgroups(newData.uniqueFlyTrajID);

% Format data for splitapply
varList = {'xCOM', 'yCOM'};
if smoothing
    varList = strcat('smooth_', varList);
end
X = newData{:, varList};

% Simplify each trajectory using the RDP algorithm
tic;
[ Y ] = splitapply(@(x) {rdpfun(x, epsilon, minFrames)}, X, G);
fprintf('Simplified trajectories using the RDP algorithm in %f seconds.\n', toc);

%% Append RDP data to the data table
RDP = cat(1,Y{:});

newData.RDPx = RDP(:,1);
newData.RDPy = RDP(:,2);
newData.prevRDPx = RDP(:,3);
newData.prevRDPy = RDP(:,4);
newData.nextRDPx = RDP(:,5);
newData.nextRDPy = RDP(:,6);
newData.RDPangle = RDP(:,7);

end

function [ out ] = rdpfun( x, epsilon, minFrames )
%Function to calculate RDP turns for a single trajectory (used with
%splitapply to calculate turns on a per-trajectory basis).

out = nan(size(x,1), 7);
notNan = ~any(isnan(x),2);

if nnz(notNan) >= minFrames
    % Run the RDP algorithm to reduce the number of points
    [ y ] = ramerDouglasPeucker(x(notNan,:), epsilon);
    
    % Find the directional changes
    R = diff([y(1,:); y; y(end,:)]);
    x1 = R(1:end-1,1);
    y1 = R(1:end-1,2);
    x2 = R(2:end,1);
    y2 = R(2:end,2);
    
    % Append an index
    idx = ismember(x,y, 'rows');
    
    % Append previous and next values
    y = [y, [0, 0; y(1:end-1,1:2)], [y(2:end,1:2); 0, 0]];
    
    % Calculate the angles (adjusting the sign to match desired convention)
    T = (-1)*atan2( x1.*y2 - y1.*x2, x1.*x2 + y1.*y2);
    
    % Append the angles and trigger
    out(idx,:) = [y, T];
end

end


