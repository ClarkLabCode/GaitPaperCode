function [ curvature ] = fast_2d_curvature(n,t,X,Y)
% Approximate the curvature of a parametric 2-d curve in an nice way by
% approximating x(t) and y(t) using nth order polynomials. Since we are
% trying to fit x(t) and y(t) given a single time vector and many samples,
% we only need one Vandermonde matrix, and can then solve the system all at
% once using mldivide. This way, it takes 0.1 seconds instead of 3 minutes
% for 20,000 timeseries. - JZV

tic;

% Get the number of points
N = size(X,1);

% Combine the data into one matrix
X = [X; Y]';

% Scale the time vector (for numerical stability)
x = (t - mean(t))./std(t);

% Make the vandermonde matrix
vandermonde = repmat(x, 1,n+1).^(n:-1:0);

% Fit the polynomial coefficients
p = vandermonde \ X;
p = p.';

% Compute the first and second derivatives
dp = p(:,1:n) .* (n:-1:1);
ddp = dp(:,1:n-1) .* (n-1:-1:1);

% Quickly evaluate the derivatives at the zero point
% In the words of MATLAB's polyval.m,
% "Make it scream for scalar x.  Polynomial evaluation can be
% implemented as a recursive digital filter."
df = filter(1,[1 0],dp);
df = df(:,n);
ddf = filter(1,[1 0],ddp);
ddf = ddf(:,n-1);

% Separate x and y data
fx = df(1:N);
fy = df(N+1:end);
fxx = ddf(1:N);
fyy = ddf(N+1:end);

% Compute the curvature
curvature = abs(fx.*fyy - fy.*fxx) ./ ((fx.^2+fy.^2).^(3/2));

fprintf('Computed curvature in %f seconds.\n', toc);
end