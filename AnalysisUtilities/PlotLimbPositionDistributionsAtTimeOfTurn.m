function PlotLimbPositionDistributionsAtTimeOfTurn(newData)

% Set the minimum and maximum forward speed
vfMin = 0.5;
vfMax = Inf;
vrMin = 0;
vrMax = Inf;

%% Extract the needed data

vfMin = 1;
limbList = {'L1','L2','L3','R1','R2','R3'};
limbVarListX = strcat(limbList, '_xPlot_mm');
limbVarListY = strcat(limbList, '_yPlot_mm');

% Define an indexing vector
vf = newData.smooth_forwardSpeed_mmPerSec;
vr = abs(rad2deg(newData.smooth_angVel_radPerSec));
idx = (vf > vfMin) & (vf < vfMax) & (vr > vrMin) & (vr < vrMax) & (newData.yawExtremum ~= 0);

% Extract the necessary variables
turnDirection = newData.yawExtremum(idx);
L = newData{idx, [limbVarListX, limbVarListY]};
ID = newData.videoID(idx);

% Remove NaN values
notNan = ~any(isnan(L),2);
L = L(notNan,:);
ID = ID(notNan,:);
turnDirection = turnDirection(notNan);

% Fold left and right turns
L(turnDirection<0,:) = L(turnDirection<0,[4,5,6,1,2,3,10,11,12,7,8,9]);
L(turnDirection<0, 1:6) = - L(turnDirection<0, 1:6);

% Make forward and rightward positive
L = -L;

[ ~, cmp, ~ ] = MakeTurningColormaps();

% Define limb list
limbList = {'O1','O2','O3','I1','I2','I3'};

%% Marginal distributions

xEdges = (-3:0.05:3)';
yEdges = (-3:0.05:3)';

xCenters = xEdges(1:end-1) + diff(xEdges)/2;
yCenters = yEdges(1:end-1) + diff(yEdges)/2;

nx = length(xCenters);
ny = length(yCenters);

fx = nan(nx, 6);
fy = nan(ny, 6);
for ind = 1:6
   fx(:,ind) = histcounts(L(:,ind), xEdges, 'normalization','pdf');
   fy(:,ind) = histcounts(L(:,ind+6), yEdges, 'normalization','pdf');
end

MakeFigure;
plot(xCenters, fx, 'linewidth', 2);
xlabel('r_{\perp} (mm)');
ylabel('pdf (mm^{-1})');
legend(limbList);
axis('square');
ConfAxis('fontSize', 16);
title('limb positions at yaw extremum');
ylim([0 4]);
yticks(0:1:4);

MakeFigure;
plot(yCenters, fy, 'linewidth', 2);
xlabel('r_{||} (mm)');
ylabel('pdf (mm^{-1})');
legend(limbList);
axis('square');
ConfAxis('fontSize', 16);
title('limb positions at yaw extremum');
ylim([0 2]);
yticks(0:1/2:2);

%% Joint distributions

fxy = nan(nx,ny,6);
for ind = 1:6
    fxy(:,:,ind) = histcounts2(L(:,ind), L(:,ind+6), xEdges, yEdges, 'normalization','pdf');
end

MakeFigure;
contour(xCenters, yCenters, sum(fxy,3)', 20, 'linewidth', 2);

axis('equal');
cbar = colorbar;
ylabel(cbar, 'pdf (mm^{-2})');
xlabel('r_{\perp} (mm)');
ylabel('r_{||} (mm)');
xlim([-3 3]);
ylim([-3 3]);
ConfAxis('fontSize', 16);
title('limb positions at yaw extremum');
colormap(cmp);
caxis([0 6]);
cbar.Ticks = 0:2:6;

end