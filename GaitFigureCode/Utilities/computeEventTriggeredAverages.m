function [ A, M ] = computeEventTriggeredAverages( newData, trigger, winlen, varList )
%This function computes event triggered averages of any set of variables in
%the data table given a trigger and window length. 

% A --> Cell Array of the Raw Data
% M --> Table of averages

if ~isa(varList,'cell') && (isa(varList, 'string') || isa(varList, 'char')) 
    varList = {varList};
end
if ~isa(newData, 'table')
   error('Data must be input as a table.'); 
end

idx = ~cellfun(@(x) any(strcmp(x, newData.Properties.VariableNames)), varList);
if any(idx)
    error('Data table does not contain a variable named %s.\n', varList{idx});
end
clearvars idx;

% Find the groups
[G, ~] = findgroups(newData.uniqueFlyTrajID);

% Remove trigger points that are too close to the ends of trajectories
[trigger] = splitapply(@(x) {removeTrajectoryEndTriggers(x, winlen)}, trigger, G);
trigger = cat(1, trigger{:});

% Slice the data
[ A ] = varfun(@(y) splitapply(@(x){getSlices(x, winlen)}, [trigger, y], G), newData, 'InputVariables', varList, 'OutputFormat','table');
A = varfun(@(x) cat(2, x{:}), A,'output','cell');

% Compute averages and format output
C = cell2mat(cellfun(@(x) nanmean(x,2), A,'uniformoutput',false));
M = array2table([(-winlen:winlen)', C], 'VariableNames',[{'frame'}, varList]);
end

function [ trigger ] = removeTrajectoryEndTriggers(trigger, winlen)
% Remove triggering points that are too close to the ends of
% trajectories
k = find(trigger);
k = [k(k-winlen <=0); k(k+winlen >= length(trigger)); k(find(diff(k)<=winlen)+1)];
trigger(k) = 0;
end

function [ S ] = getSlices(X, winlen)
trigger = X(:,1);
stim = X(:,2:end);

if any(trigger)
    k = find(trigger);
    
    % Pad the stimulus with NaN values as needed
    if any((k+winlen) > size(stim,1))
        stim = [stim; nan(k((k+winlen) > size(stim,1))+winlen-size(stim,1), size(stim,2))];
    end
    if any((k-winlen) <=0 )
        n = abs(k((k-winlen)<=0)-winlen)+1;
        stim = [nan(n, size(stim,2)); stim];
        k = k+n;
    end
    
    % Cut out the slices
    slices = arrayfun(@(x) stim(x-winlen:x+winlen,:), ...
        k, ...
        'UniformOutput', 0);
    slices = cellfun(@(x) cat(1,x(:)), slices, 'uniformoutput',false);
    S = cat(2, slices{:});
else
    S = [];
end

end