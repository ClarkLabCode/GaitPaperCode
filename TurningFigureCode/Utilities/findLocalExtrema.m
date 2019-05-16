function [ trigger ] = findLocalExtrema(kind, x, groups, winlen, minHeight, minProminence)
%FINDLOCALEXTREMA: A function to find local extrema in the vector x for
%each of the groups defined by the vector groups.

if length(x) ~= length(groups)
    error('The length of the input vector must equal the length of the grouping vector');
end

% Find the groups
[G, ~] = findgroups(groups);

% Account for different kinds of extrema
switch kind
    case 'maxima'
        % Do nothing
    case 'minima'
        x = (-1).*x;
    case 'all'
        x = abs(x);
    otherwise
        error('%s is an invalid kind of extremum.', kind);
end

% Find the extrema
warning('off', 'signal:findpeaks:largeMinPeakHeight');
[trigger] = splitapply(@(a){findextrema(a, winlen, minHeight, minProminence)}, [zeros(length(x),1), x], G);
warning('on', 'signal:findpeaks:largeMinPeakHeight');

% Combine data for output
trigger = cat(1, trigger{:});

end

function [ y ] = findextrema(x, winlen, minHeight, minProminence)
y = logical(x(:,1));
if size(x,1) >= 15
    if any(x(:,2) >= minHeight)
        % Use findpeaks to find extrema
        [~, locs] = findpeaks(x(:,2), 'MinPeakDistance', min(length(x)/2, winlen), 'MinPeakHeight', minHeight, 'MinPeakProminence', minProminence);
        if any(locs)
            y(locs) = true;
        end
    end
end
end