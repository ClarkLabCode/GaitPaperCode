function PlotAverageModulationAsFunctionOfCentroidKinematics(V, X, ID, labelStr, legendStr, clim, vfInterval, mean_type )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here


lineColor = lines(size(X,2));

if nargin >= 8 && ~isempty(mean_type) && strcmpi(mean_type, 'circular')
    accumfun = @(x) rad2deg(atan(nansum(sin(mod(x,2*pi)))./nansum(cos(mod(x,2*pi)))));
else
    accumfun = @nanmean;
end

fontSize = 16;

nboot = 1000;
bootAlpha = 0.05;

tl = vfInterval(1);
tu = vfInterval(2);

%% Define bins

vrMax = 425;
if min(V(:,1)) >= 0
    vrBinEdges = 0:20:vrMax;
else
    vrBinEdges = -vrMax:20:vrMax;
end
vfBinEdges = 0:0.5:32;

vpBinEdges = -12:0.5:12;

vfBinCenters = (vfBinEdges(1:end-1))';
vrBinCenters = (vrBinEdges(1:end-1))';
vpBinCenters = (vpBinEdges(1:end-1))';

%% Discretize the centroid dynamics

vrDiscr = discretize( V(:,1), vrBinEdges );
vfDiscr = discretize( V(:,2), vfBinEdges );
vpDiscr = discretize( V(:,3), vpBinEdges );

% Check for zero values
vfDiscr(vfDiscr==0) = NaN;
vrDiscr(vrDiscr==0) = NaN;
vpDiscr(vpDiscr==0) = NaN;

% Remove NaN values
notNan = ~any(isnan([vrDiscr, vfDiscr, vpDiscr, ID]),2);
vrDiscr = vrDiscr(notNan);
vfDiscr = vfDiscr(notNan);
vpDiscr = vpDiscr(notNan);
V = V(notNan,:);
X = X(notNan,:);
ID = ID(notNan);

% Get the expected number of indices
nr = length(vrBinCenters);
nf = length(vfBinCenters);
np = length(vpBinCenters);
ni = length(unique(ID));

% Find groups
vrDiscr = findgroups(vrDiscr);
vfDiscr = findgroups(vfDiscr);
vpDiscr = findgroups(vpDiscr);
ID = findgroups(ID);

% Check for empty bins
if (length(unique(vrDiscr)) ~= nr)
    error('There are empty v_{r} bins.');
end
if (length(unique(vfDiscr)) ~= nf)
    error('There are empty v_{||} bins.');
end
if (length(unique(vpDiscr)) ~= np)
    error('There are empty v_{\perp} bins.');
end
if (length(unique(ID)) ~= ni)
    error('There are empty IDs.');
end

%% Plot the singly conditioned means
clearvars sp;

% Compute stats over IDs, conditioned on yaw
if size(X,2)==1
    s = accumarray([vrDiscr, ID], X, [nr,ni], accumfun);
else
    s = nan(nr,ni, size(X,2));
    for ind = 1:size(X,2)
        s(:,:,ind) = accumarray([vrDiscr, ID], X(:,ind), [nr,ni], accumfun);
    end
end

% Compute mean and CIs over IDs
a = squeeze(nanmean(s,2));
ci = bootci(nboot, {@(f) nanmean(f, 1), permute(s, [2,1,3])}, 'alpha', bootAlpha);

MakeFigure;
sp(1) = subplot(1,3,1);
PlotAsymmetricErrorPatch(vrBinCenters, a, squeeze(ci(1,:,:,:)), squeeze(ci(2,:,:,:)), lineColor);
axis('square');
xlabel('v_{r} (\circ{\cdot}s^{-1})');
legend(legendStr);
ylabel(labelStr);
ConfAxis('fontSize', fontSize);
xlim([min(vrBinCenters), max(vrBinCenters)]);
ylim([min(clim) max(clim)]);

% Compute stats over IDs, conditioned on forward velocity
if size(X,2)==1
    s = accumarray([vfDiscr, ID], X, [nf,ni], accumfun);
else
    s = nan(nf,ni, size(X,2));
    for ind = 1:size(X,2)
        s(:,:,ind) = accumarray([vfDiscr, ID], X(:,ind), [nf,ni], accumfun);
    end
end

% Compute mean and SEM over IDs
a = squeeze(nanmean(s,2));
ci = bootci(nboot, {@(f) nanmean(f, 1), permute(s, [2,1,3])}, 'alpha', bootAlpha);
sp(2) = subplot(1,3,2);
% PlotErrorPatch(vfBinCenters, a, da, lineColor);
PlotAsymmetricErrorPatch(vfBinCenters, a, squeeze(ci(1,:,:,:)), squeeze(ci(2,:,:,:)), lineColor);
axis('square');
xlabel('v_{||} (mm{\cdot}s^{-1})');
ylabel(labelStr);
ConfAxis('fontSize', fontSize);
xlim([min(vfBinCenters), max(vfBinCenters)]);
ylim([min(clim) max(clim)]);

% Compute stats over IDs, conditioned on lateral velocity
if size(X,2)==1
    s = accumarray([vpDiscr, ID], X, [np,ni], accumfun);
else
    s = nan(np,ni, size(X,2));
    for ind = 1:size(X,2)
        s(:,:,ind) = accumarray([vpDiscr, ID], X(:,ind), [np,ni], accumfun);
    end
end

% Compute mean and SEM over IDs
a = squeeze(nanmean(s,2));
ci = bootci(nboot, {@(f) nanmean(f, 1), permute(s, [2,1,3])}, 'alpha', bootAlpha);

sp(3) = subplot(1,3,3);
PlotAsymmetricErrorPatch(vpBinCenters, a, squeeze(ci(1,:,:,:)), squeeze(ci(2,:,:,:)), lineColor);
axis('square');
xlabel('v_{\perp} (mm{\cdot}s^{-1})');
ylabel(labelStr);
ConfAxis('fontSize', fontSize);
xlim([min(vpBinCenters), max(vpBinCenters)]);
ylim([min(clim) max(clim)]);

%% Plot the mean as a function of yaw, conditioned on forward walking
clearvars sp;


idx = (V(:,2) > tl) & (V(:,2) < tu);
titleStr = sprintf('%d < v_{||} < %d mm{\\cdot}s^{-1}; N = %d frames', tl, tu, nnz(idx));

% Compute stats over IDs, conditioned on yaw
if size(X,2)==1
    s = accumarray([vrDiscr(idx), ID(idx)], X(idx,:), [nr,ni], accumfun);
else
    s = nan(nr,ni, size(X,2));
    for ind = 1:size(X,2)
        s(:,:,ind) = accumarray([vrDiscr(idx), ID(idx)], X(idx,ind), [nr,ni], accumfun);
    end
end

% Compute mean and SEM over IDs
[~,noTurn] = min(abs(vrBinCenters));

if strcmp(func2str(accumfun), 'nanmean')
    sn = s(noTurn,:,:);
    s = (s - sn)./sn;
end

a = squeeze(nanmean(s,2));
ci = bootci(nboot, {@(f) nanmean(f, 1), permute(s, [2,1,3])}, 'alpha', bootAlpha);

MakeFigure;
PlotAsymmetricErrorPatch(vrBinCenters, a, squeeze(ci(1,:,:,:)), squeeze(ci(2,:,:,:)), lineColor);
axis('square');
xlabel('v_{r} (\circ{\cdot}s^{-1})');
legend(legendStr, 'location', 'eastoutside');

if strcmp(func2str(accumfun), 'nanmean')
    ylabel(sprintf('fractional change relative to no turning in %s',labelStr));
else
    ylabel(sprintf('%s',labelStr));
end
ConfAxis('fontSize', fontSize);
xlim([min(vrBinCenters), max(vrBinCenters)]);
title(titleStr);

end

