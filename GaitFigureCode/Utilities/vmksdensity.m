function [f] = vmksdensity(x, q, nu)
%VMKSDENSITY: Kernel density estimation for the von Mises distribution

% Citations:
%   M. Oliveira, R.M. Crujeiras, A. Rodríguez-Casal. A plug-in rule for
%   bandwidth selection in circular density estimation. Computational
%   Statistics and Data Analysis 56 (2012) 3898–3908.
%   https://www.sciencedirect.com/science/article/pii/S0167947312002204.

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

% Get the number of data points
n = length(x);

% Set the bandwidth
if nargin < 3 || isempty(nu)
    %     % Get the MLE of the von Mises parameter kappa
    % k = vmkappa(x);
    
    %     % Least-asymptotic-mean-squared-error estimate
    %     num = 3*n*(k^2)*besseli(2,2*k);
    %     den = 4*sqrt(pi)*(besseli(0,k)^2);
    %     nu = (num/den)^(2/5);
    
    % Silverman's normal bandwidth rule-of-thumb
    sd = sqrt(2*(1-abs(nanmean(exp(1i*x)))));
    nu = (n^(1/5))/(0.9*sd);
end

% Compute the kernel density estimate of the PDF
f = sum(exp(nu*cos(q-x')),2) / (n*2*pi*besseli(0,nu));

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

