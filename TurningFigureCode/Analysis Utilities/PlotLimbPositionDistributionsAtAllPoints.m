function PlotLimbPositionDistributionsAtAllPoints(newData)

%% Extract the needed data

vfMin = 1;
limbList = {'L1','L2','L3','R1','R2','R3'};
limbVarListX = strcat(limbList, '_xPlot_mm');
limbVarListY = strcat(limbList, '_yPlot_mm');
L = newData{:, [limbVarListX, limbVarListY]};
vf = newData.smooth_forwardSpeed_mmPerSec;
L = L(vf>vfMin,:);

% Make forward and rightward positive
L = -L;

[ ~, cmp, ~ ] = MakeTurningColormaps();

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
title(sprintf('limb positions given v_{||} > %d mm/s', vfMin));
ylim([0 4]);
yticks(0:1:4);

MakeFigure;
plot(yCenters, fy, 'linewidth', 2);
xlabel('r_{||} (mm)');
ylabel('pdf (mm^{-1})');
legend(limbList);
axis('square');
ConfAxis('fontSize', 16);
title(sprintf('limb positions given v_{||} > %d mm/s', vfMin));
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
title(sprintf('limb positions given v_{||} > %d mm/s', vfMin));
colormap(cmp);
caxis([0 6]);
cbar.Ticks = 0:2:6;
end