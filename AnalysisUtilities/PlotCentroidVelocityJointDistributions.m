function PlotCentroidVelocityJointDistributions( newData, speedBinEdges, vfBinEdges, vpBinEdges, vrBinEdges, cmap, logscale, linewidth, smoothing)
%PlotCentroidVelocityComponentDistributions: A function to visualize
%joint distributions of the centroid dynamical variables. 

%% Validate inputs

if ~exist('speedBinEdges','var') || isempty(speedBinEdges)
    speedBinEdges = linspace(0, 50, 200);
end
speedBinCenters = speedBinEdges(1:end-1) + diff(speedBinEdges)/2;

if ~exist('vrBinEdges','var') || isempty(vrBinEdges)
    vrBinEdges = linspace(-300, 300, 100);
end
vrBinCenters = vrBinEdges(1:end-1) + diff(vrBinEdges)/2;

if ~exist('vfBinEdges','var') || isempty(vfBinEdges)
    vfBinEdges = linspace(-5, 35, 100);
end
vfBinCenters = vfBinEdges(1:end-1) + diff(vfBinEdges)/2;

if ~exist('vpBinEdges','var') || isempty(vpBinEdges)
    vpBinEdges = linspace(-10, 10, 100);
end
vpBinCenters = vpBinEdges(1:end-1) + diff(vpBinEdges)/2;

if ~exist('cmap','var') || isempty(cmap)
    cmap = viridis(256);
end

if ~exist('logscale','var') || isempty(logscale)
    logscale = true;
end

if ~exist('linewidth','var') || isempty(linewidth)
    linewidth = 2;
end

if ~exist('smoothing','var') || isempty(smoothing)
    smoothing = true;
end

% Set the number of contour lines
numLvl = 5;

%% Extract the data

% Define the list of dynamical variables
varList = {'linVel_mmPerSec','angVel_radPerSec','forwardSpeed_mmPerSec','translationalSpeed_mmPerSec'};

% Check whether smoothed variables are desired
if smoothing
   varList = strcat('smooth_', varList);    
end

% Extract the data
V = newData{:, varList};

% Remove NaN values
isNotNanIdx = ~any(isnan(V),2);
V = V(isNotNanIdx,:);

% Split up the velocity components
% (and convert the angular velocity to deg/s from rad/s)
speed = V(:,1);
vr    = rad2deg(V(:,2));
vf    = V(:,3);
vp    = V(:,4);

clearvars V;

%% Calculate the single-variable histograms

% Speed
N_speed = histcounts(speed, speedBinEdges, 'normalization','pdf');

% Angular velocity
N_ang = histcounts(vr, vrBinEdges, 'normalization','pdf');

% Forward velocity
N_fwd = histcounts(vf, vfBinEdges, 'normalization','pdf');

% Perpendicular velocity
N_perp = histcounts(vp, vpBinEdges, 'normalization','pdf');

% Perpendicular velocity versus forward velocity
N = histcounts2(vf, vp, vfBinEdges, vpBinEdges, 'normalization','pdf');

figure('Position',[200,500,1000,1000],'WindowStyle','docked');
if logscale
    hold on;
    imagesc(vfBinCenters, vpBinCenters, log10(N)');
    [xList, yList] = meshgrid(vfBinCenters, vpBinCenters);
    contour(xList, yList, log10(N)',numLvl,'edgecolor', 'k');
else
    hold on;
    imagesc(vfBinCenters, vpBinCenters, N');
    [xList, yList] = meshgrid(vfBinCenters, vpBinCenters);
    contour(xList, yList, N',numLvl,'edgecolor', 'k');

end

axis('xy','square','tight');

colormap(cmap);
cbar = colorbar;
cbar.Location = 'northoutside';

xlabel('v_{||} (mm/s)');
ylabel('v_{\perp} (mm/s)');

if logscale
   ylabel(cbar, 'log_{10} P( v_{||}, v_{\perp} )');
%    caxis([-6.5, -2.5]);
else
    ylabel(cbar, 'P( v_{||}, v_{\perp} )'); 
end

ConfAxis;

% Angular velocity versus forward velocity
N = histcounts2(vf, vr, vfBinEdges, vrBinEdges, 'normalization','pdf');

figure('Position',[200,500,1000,1000],'WindowStyle','docked');
if logscale
    hold on;
    imagesc(vfBinCenters, vrBinCenters, log10(N)');
    [xList, yList] = meshgrid(vfBinCenters, vrBinCenters);
    contour(xList, yList, log10(N)',numLvl,'edgecolor', 'k');
else
    hold on;
    imagesc(vfBinCenters, vrBinCenters, (N)');
    [xList, yList] = meshgrid(vfBinCenters, vrBinCenters);
    contour(xList, yList, (N)',numLvl,'edgecolor', 'k');    
end

axis('xy','square','tight');

colormap(cmap);
cbar = colorbar;
cbar.Location = 'northoutside';

xlabel('v_{||} (mm/s)');
ylabel('v_{r} (\circ/s)');

if logscale
   ylabel(cbar, 'log_{10} P( v_{||}, v_{r} )');
%    caxis([-6.5, -2.5]);
else
    ylabel(cbar, 'P( v_{||}, v_{r} )'); 
end

ConfAxis;

% Angular velocity versus perpendicular velocity
N = histcounts2(vp, vr, vpBinEdges, vrBinEdges, 'normalization','pdf');

figure('Position',[200,500,1000,1000],'WindowStyle','docked');
if logscale
    hold on;
    imagesc(vpBinCenters, vrBinCenters, log10(N)'); 
    [xList, yList] = meshgrid(vpBinCenters, vrBinCenters);
    contour(xList, yList, log10(N)',numLvl,'edgecolor', 'k');
else
    hold on;
    imagesc(vpBinCenters, vrBinCenters, (N)');
    [xList, yList] = meshgrid(vpBinCenters, vrBinCenters);
    contour(xList, yList, (N)',numLvl,'edgecolor', 'k');
end

axis('xy','square','tight');

colormap(cmap);
cbar = colorbar;
cbar.Location = 'northoutside';

xlabel('v_{\perp} (mm/s)');
ylabel('v_{r} (\circ/s)');

if logscale
   ylabel(cbar, 'log_{10} P( v_{\perp}, v_{r} )');
else
    ylabel(cbar, 'P( v_{\perp}, v_{r} )'); 
end

ConfAxis;


end