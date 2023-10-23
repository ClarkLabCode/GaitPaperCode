%% Convert limb UMAP into cylindrical coordinates

[th, rho, z] = cart2pol(limbUMAP(:,2), limbUMAP(:,3), limbUMAP(:,1));

% Convert angular coordinate into cycles
th = mod(th + pi, 2*pi)/(2*pi);

[ cmpBlueRed, cmpRed, cmpBlue ] = MakeTurningColormaps();

%% Marginal distributions of UMAP cylindrical coordinates

% Phase
q = (0:0.01:1)';
f = bqksdensity(2*pi*th, 2*pi*q);
MakeFigure;
plot(q,f, 'linewidth', 2);
xlabel('UMAP phase (cycles)');
ylabel('pdf (1/cycles)');
axis('square');
ConfAxis('fontSize', 16);
ylim([0 0.2]);

% Radius
q = (0:0.01:3)';
f = ksdensity(rho, q, 'bandwidth', 1/16);
MakeFigure;
plot(q,f, 'linewidth', 2);
xlabel('UMAP radius (arb. units)');
ylabel('pdf (1/arb. units)');
axis('square');
ConfAxis('fontSize', 16);

% z
q = (-6:0.01:8)';
f = ksdensity(z, q, 'bandwidth', 1/16);
MakeFigure;
plot(q,f, 'linewidth', 2);
xlabel('UMAP z (arb. units)');
ylabel('pdf (1/arb. units)');
axis('square');
ConfAxis('fontSize', 16);

%% Joint distribution of UMAP radius and z

binEdgesR = (0:0.05:3)';
binEdgesZ = (-4:0.05:3)';
binCentersR = binEdgesR(1:end-1) + diff(binEdgesR)/2;
binCentersZ = binEdgesZ(1:end-1) + diff(binEdgesZ)/2;

f = histcounts2(z, rho, binEdgesZ, binEdgesR, 'normalization','pdf');
MakeFigure;
imagesc(binCentersZ, binCentersR,f');
hold on;
contour(binCentersZ, binCentersR', f', 0:1/2:4, 'EdgeColor','k', 'LineWidth', 2);
axis('xy', 'equal','tight');
xlabel('UMAP z (arb. units)');
ylabel('UMAP radius (arb. units)');
cbar = colorbar;
ylabel(cbar, 'pdf (1/arb. units^2)');
yticks(0:3);
caxis([0 2]);
cbar.Ticks = 0:4;
colormap(cmpRed);
ConfAxis('fontSize', 16);

%% Joint distribution of UMAP and L2 phases

binEdges = (0:0.01:1)';
binCenters = binEdges(1:end-1) + diff(binEdges)/2;
f = histcounts2(th, squeeze(Phi(t==0,2,:)), binEdges, binEdges, 'normalization','pdf');
MakeFigure;

imagesc(binCenters,binCenters,f');
hold on;
contour(binCenters, binCenters', f', 0:2:10, 'EdgeColor','k', 'LineWidth', 2);
axis('xy', 'square','tight');
xlabel('UMAP phase (cycles)');
ylabel('L_2 phase (cycles)');
cbar = colorbar;
colormap(cmpRed);
ylabel(cbar, 'pdf (1/cycles^2)');
caxis([0 10]);
ConfAxis('fontSize', 16);

%% Joint distribution of UMAP z and mean frequency

meanFreq = squeeze(mean(mean(dPhi(:,:,:),2),1));

binEdgesF = (0:0.25:25)';
binEdgesZ = (-4:0.05:3)';
binCentersF = binEdgesF(1:end-1) + diff(binEdgesF)/2;
binCentersZ = binEdgesZ(1:end-1) + diff(binEdgesZ)/2;

f = histcounts2(z, meanFreq, binEdgesZ, binEdgesF, 'normalization','pdf');
MakeFigure;
imagesc(binCentersZ, binCentersF,f');
hold on;
contour(binCentersZ, binCentersF', f', 0:0.1:0.3, 'EdgeColor','k', 'LineWidth', 2);
axis('xy', 'square','tight');
xlabel('UMAP z (arb. units)');
ylabel('mean frequency (Hz)');
cbar = colorbar;
ylabel(cbar, 'pdf (1/arb. units/Hz)');
caxis([0 0.3]);
colormap(cmpRed);
ConfAxis('fontSize', 16);