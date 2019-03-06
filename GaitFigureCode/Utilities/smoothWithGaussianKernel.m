function [ Y, w ] = smoothWithGaussianKernel( newData, varList, N, varargin )
%SMOOTHWITHGAUSSIANKERNEL: A function to smooth a set of variables from a
%data table using a gaussian kernel

%% Create the kernel

if (mod(N,1)>0) || ~mod(N,2)
    error('Window length must be an odd integer.');
end
m = (N-1)/2;
n = -m:m;

% Process user input
arga = find(strcmpi(varargin, 'alpha'), 1);
args = find(strcmpi(varargin, 'sigma'), 1);
if isempty(arga) && isempty(args)
    a = 2.5;
elseif isempty(arga) && ~isempty(args)
    a = (N-1)/(2*varargin{args+1});
elseif isempty(args) && ~isempty(arga)
    a = varargin{arga+1};
else
    error('Values for alpha and sigma values cannot be provided simultaneously.');
end

argm = find(strcmpi(varargin, 'method'),1);
if isempty(argm)
    method = 'filter';
else
    method = varargin{argm+1};
    if ~any(strcmpi(method, {'filter','filtfilt'}))
        error('Invalid filtering method: %s.');
    end
end

% Compute the kernel
w = exp(-(1/2)*(a*n/(N/2)).^2);

% Normalize the kernel
w = w ./ sum(w);

%% Smooth the data

% Find the groups corresponding to trajectories
[G, ~] = findgroups(newData.uniqueFlyTrajID);

% Format data for splitapply
X = newData{:, varList};

% Smooth each trajectory
switch method
    case 'filter'
        [ Y ] = splitapply(@(x) {smoothfun_filter(x, w)}, X, G);
    case 'filtfilt'
        [ Y ] = splitapply(@(x) {smoothfun_filtfilt(x, w)}, X, G);
end

% Combine data for output
Y = cat(1, Y{:});

end

function [ out ] = smoothfun_filter(x,w)
% This can only handle cases in which NaN values are concentrated at the
% beginning and end of trajectories
out = nan(size(x));
if size(x,1) >= length(w)
    xnan = any(isnan(x),2);
    if nnz(~xnan) > length(w)
        out(~xnan,:) = filter(w,1,x(~xnan,:));
    end
end
end

function [ out ] = smoothfun_filtfilt(x,w)
% This can only handle cases in which NaN values are concentrated at the
% beginning and end of trajectories
out = nan(size(x));
if size(x,1) >= 3*length(w)-1
    xnan = any(isnan(x),2);
    if nnz(~xnan) > 3*length(w)-1
        out(~xnan,:) = filtfilt(w,1,x(~xnan,:));
    end
end
end
