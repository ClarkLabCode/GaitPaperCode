function [S] = getSlicesOneSidedEither(newData, varList, winlen, trigger, side)
% This function grabs one-sided trajectories 

% Minor input validation
if ~isa(varList,'cell') && (isa(varList, 'string') || isa(varList, 'char'))
    varList = {varList};
end
idx = ~cellfun(@(x) any(strcmp(x, newData.Properties.VariableNames)), varList);
if any(idx)
    error('Data table does not contain a variable named %s.\n', varList{idx});
end
clearvars idx;

% Find the groups
[G, ~] = findgroups(newData.uniqueFlyTrajID);

% Get numeric data
X = newData{:, varList};

% Get slices
[ S ] = splitapply(@(x,y) {slicefun(x, y, winlen, side)}, trigger, X, G);
S = cat(1, S{:});

end

function [ slices ] = slicefun(t, X, w, side)

slices = [];
if any(t)
    
    % Find the triggering points
    k = find(t);
    
    if strcmp (side, 'pre')
        % Remove triggering points that are too close to the end of the trajectory
        k((k-w) < 1) = [];
        
        % Check to make sure any triggering points remain
        if any(k)
            
            % Cut out the slices
            slices = arrayfun(@(x) X(x-w:x,:), k, 'UniformOutput', false);
        end
        
    elseif strcmp (side, 'post')
        % Remove triggering points that are too close to the end of the trajectory
        k((k + w) > size(X,1)) = [];
        
        
        % Check to make sure any triggering points remain
        if any(k)
            
            % Cut out the slices
            slices = arrayfun(@(x) X(x:x+w,:), k, 'UniformOutput', false);
            
        end
        
    else
        error('Invalid directional input.');
    end
end

end
