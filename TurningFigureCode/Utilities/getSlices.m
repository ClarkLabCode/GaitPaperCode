function [ S, trigger, seq ] = getSlices( newData, varList, winlen, trigger, allowOverlap )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

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

if ~allowOverlap
    % Remove trigger points that are too close to one another
    [trigger] = splitapply(@(x) {removeTooCloseTriggers(x, winlen)}, trigger, G);
    trigger = cat(1, trigger{:});
end

% Get slices
[ S ] = splitapply(@(t, x) {slicefun(t, x, winlen)}, trigger, X, G);
S = cat(1, S{:});

% For each slice, get the uniqueFlyTrajID, starting frame, and ending frame
if nargout == 3
    X_seq = newData{:,{'uniqueFlyTrajID', 'Frame'}};
    [ seq ] = splitapply(@(t,x) {slicefun(t,x, winlen)}, trigger, X_seq, G);
    seq = cat(1,seq{:});
    seq = cellfun(@(x) [x(:,1), repmat([min(x(:,2)), max(x(:,2))], size(x,1), 1)], seq, 'uniformoutput', false);
end

end

function [ trigger ] = removeTooCloseTriggers(trigger, winlen)
% Remove triggering points that are too close to one another
k = find(trigger);
k = k(find(diff(k)<=winlen)+1);
trigger(k) = 0;
end

function [ S ] = slicefun(trigger, X, winlen)

if any(trigger)
    % Find the triggering points
    k = find(trigger);
    
    % Pad the stimulus with NaN values as needed
    padAfter = max(k) + winlen - size(X,1);
    padBefore = abs(min(k) - winlen) + 1;
    if padAfter > 0
        X = [X; nan(padAfter, size(X,2))];
    end
    if padBefore > 0
        X = [nan(padBefore, size(X,2)); X];
        k = k + padBefore;
    end
    
    % Cut out the slices
    slices = arrayfun(@(x) X(x-winlen:x+winlen,:), k, 'UniformOutput', false);
    S = slices;
else
    S = [];
end

end