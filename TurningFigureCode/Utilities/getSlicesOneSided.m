function [S] = getSlicesOneSided(newData, varList, winlen, trigger)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

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
[ S ] = splitapply(@(x,y) {slicefun(x, y, winlen)}, trigger, X, G);
S = cat(1, S{:});

end

function [ slices ] = slicefun(t, X, w)

slices = [];
if any(t)
    
    % Find the triggering points
    k = find(t);
    
    % Remove triggering points that are too close to the end of the trajectory
    k((k + w) > size(X,1)) = [];
    
    % Check to make sure any triggering points remain
    if any(k)
        
        % Cut out the slices
        slices = arrayfun(@(x) X(x:x+w,:), k, 'UniformOutput', false);
    end
end

end
