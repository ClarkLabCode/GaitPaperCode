function [ R, t_ms ] = PlotVelocityAutocorrelation(newData, maxTau, smoothing, fps, corder, plotci, upsampling)

% References:
%   Marple, S. L. (1999). Estimating group delay and phase delay via
%   discrete-time" analytic" cross-correlation. IEEE transactions on signal
%   processing, 47(9), 2604-2607.
%

if ~exist('maxTau','var') || isempty(maxTau)
    maxTau = 150;
end

if ~exist('smoothing','var') || isempty(smoothing)
    smoothing = false;
end

if ~exist('fps','var') || isempty(fps)
    fps = 150;
end

if ~exist('corder','var') || isempty(corder)
    corder = linspecer(3);
end

if ~exist('plotci','var') || isempty(plotci)
    plotci = false;
end

if ~exist('upsampling','var') || isempty(upsampling)
    upsampling = 1;
end

varList = {'angVel_radPerSec', 'forwardSpeed_mmPerSec','translationalSpeed_mmPerSec'};
legendStr = {'v_{r}', 'v_{||}','v_{\perp}'};

if smoothing
    varList = strcat('smooth_', varList);
end

% Get a vector of delays
t = (0:(1/upsampling):maxTau)';
t_ms = t / (fps/1000);

% Set parameters for bootstrapping
ci_level = 0.01;
nb = 1000;

%% Compute autocorrelation functions

% Extract the needed data
X = newData{:, varList};
id = newData.uniqueFlyTrajID;

% Remove missing values
notNan = ~any(isnan(X),2);
X = X(notNan,:);
id = id(notNan);
g = findgroups(id);

% Compute autocorrelation functions
tic;
C = splitapply(@(x) {acf(x, maxTau, upsampling)}, X, g);
C = cellfun(@(x) reshape(x, 1, []), C, 'uni', false);
C = cell2mat(C);
fprintf('Computed autocorrelation functions in %f seconds.\n', toc);

% Compute the average of the divergence statistic (over trajectories)
R = reshape(nanmean(C), [], length(varList));

%% Compute bootstrapped confidence intervals

if plotci
    % Open a parallel pool
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        poolobj = parpool('local');
    end
    
    % Compute bootstrapped 100*(1-ci_level) % confidence intervals
    tic;
    opts = statset('useparallel', true);
    R_ci = bootci(nb, {@nanmean, C}, 'alpha', ci_level, 'options', opts);
    R_ci_lower = reshape(R_ci(1,:), [], length(varList));
    R_ci_upper = reshape(R_ci(2,:), [], length(varList));
    fprintf('Computed %d%% confidence intervals in %f seconds.\n', 100*(1-ci_level), toc);
end

%% Plot autocorrelation functions

figure('Position',[200,500,1000,1000],'WindowStyle','docked');
hold on;
set(gca, 'colororder', corder);
if plotci
    xx = repmat(t_ms, 1, length(varList));
    PlotAsymmetricErrorPatch(xx, R, R_ci_lower, R_ci_upper, corder);
else
    plot(t_ms, R, '-', 'linewidth', 2);
end
axis('square');
yticks(-1:0.5:1);
xlim([0, maxTau]/ (fps/1000));
xlabel('time (ms)');
ylabel('R(t)');
ConfAxis();
legend(legendStr);

end

function [ g ] = acf(x, n, m)
k = size(x,2);
X = fft(x, 2*n + 1, 1);
G = X .* conj(X);
G = [G(1,:); 2*G(2:n,:); G(n+1,:); zeros(2*m*n-n-1,k)];
g = real(ifft(G)) * m;
g = g(1:(n*m+1),:) ./ g(1,:);
end