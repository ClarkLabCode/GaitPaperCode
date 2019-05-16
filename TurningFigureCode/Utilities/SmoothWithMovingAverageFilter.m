function [ Y ] = SmoothWithMovingAverageFilter(newData, varList, span)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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