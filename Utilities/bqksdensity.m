function [ f ] = bqksdensity(x, q, scale)
%BQKSDENSITY Kernel density estimation for circular data using a quartic
%kernel. - JZV
%
% Citations:
%   Fisher, N. I. (1989). Smoothing a sample of circular data. Journal of
%   structural geology, 11(6), 775-778. 

% Validate inputs
if ~isvector(x) || isscalar(x)
    error('Input data must be a vector.');
end
if size(x,1) == 1
    x = x';
end

% Validate query points
if ~isvector(q)
    error('Query points must be input as a vector');
end
if size(q,1) == 1
    q = q';
end
if any(abs(diff(q) - mean(diff(q))) > 1e-10)
    error('Query points must be evenly spaced');
end

% Get the number of data points
n = length(x);

% Get the MLE of the von Mises parameter kappa
k = vmkappa(x);

% Compute the estimate of the smoothing parameter
h = sqrt(7)/sqrt(k)/(n^(1/5));

% Scale the smoothing parameter, if desired
if (nargin > 2) && (scale > 0)
    h = scale*h;
end

% Compute the raw kernel density estimate
d = abs(q - x');
d = min(d, 2*pi-d);
d(d>h) = NaN;
f = 0.9375 * nansum(d,2) / (n*h);

% Adjust the density
f = sqrt(1 + f / max(f)) - 1;

% Normalize the PDF
f = f ./ (sum(f) .* mean(diff(q)));

end

function [ k ] = vmkappa(x)
% Compute the MLE of the von Mises parameter kappa
% Adapted from circ_kappa from the Circular Statistics Toolbox for MATLAB

% Get the number of data points
n = length(x);

% Compute the mean resultant length
r = abs(sum(exp(1i*x)))/n;

% Compute estimate of kappa
if r < 0.53
    k = 2*r + r^3 + 5*r^5/6;
elseif r>=0.53 && r<0.85
    k = -0.4 + 1.39*r + 0.43/(1-r);
else
    k = 1/(r^3 - 4*r^2 + 3*r);
end

if n<15 && n>1
    if k < 2
        k = max(k-2*(n*k)^-1,0);
    else
        k = (n-1)^3*k/(n^3+n);
    end
end

end