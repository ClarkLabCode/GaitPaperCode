%% Show limb data covariance matrix and first two PCs

C = cov(X_L);
[ cmpBlueRed, cmpRed, cmpBlue ] = MakeTurningColormaps();
MakeFigure;
imagesc(C);
axis('xy','square','tight', 'off');
cbar = colorbar;
caxis([-1,1]);
colormap(cmpBlueRed);
ylabel(cbar, 'Covariance (standardized units)');
cbar.Ticks = [-1,0,1];
ConfAxis('fontSize', 16);

%% Decompose data using PCA

[coeff,score,latent] = pca(X_L);

%% Plot PC spectrum

MakeFigure;
plot(latent(1:25) / sum(latent), '-ok', 'MarkerFaceColor','k', 'MarkerEdgeColor','None','Linewidth', 2, 'MarkerSize', 10);
axis('square');
ylim([0 0.2]);
yticks(0:0.1:0.2);
xlabel('Principal component');
ylabel('Fraction of total variance');
ConfAxis('fontSize', 16);

%% Plot first four PCs

for ind = 1:2
    MakeFigure;
    imagesc(t_ms, 1:12, reshape(coeff(:,ind), length(t_ms), 12)');
    axis('xy','square','tight');
    xlabel('Time (ms)');
    yticks(1:12);
    yticklabels([strcat({'L1','L2','L3','R1','R2','R3'}, '_x'), strcat({'L1','L2','L3','R1','R2','R3'}, '_y')]);
    ConfAxis('fontSize', 16);
    cbar = colorbar;
    caxis([-0.2,0.2]);
    cbar.Ticks = [-0.2, 0, 0.2];
    colormap(cmpBlueRed);
    title(sprintf('PC %d', ind));
    ylabel(cbar, 'weighting (standardized units)');
end

%% Scatter the first two PCs

meanFreq = squeeze(mean(mean(dPhi(:,:,:),2),1));
MakeFigure;
scattern(score(:, 1:2),10,meanFreq, 'filled');
axis('square');
cbar = colorbar();
ylabel(cbar, 'mean frequency (Hz)');
caxis([0 15]);
colormap(cmpRed);
xlabel('PC_1 (arb. units)');
ylabel('PC_2 (arb. units)');
zlabel('PC_3 (arb. units)');
ConfAxis('fontSize', 14);

MakeFigure;
scattern(score(:,1:2), 10, squeeze(L(t==0,8,:)), 'filled');
axis('square');
caxis([-1,1]);
colormap(cmpBlueRed);
xlabel('PC_1 (arb. units)');
ylabel('PC_2 (arb. units)');
zlabel('PC_3 (arb. units)');
cbar = colorbar;
ylabel(cbar, 'L2_{||} (mm)');
ConfAxis('fontSize', 14);


MakeFigure;
scattern(score(:,1:2), 10, squeeze(L(t==0,11,:)), 'filled');
axis('square');
caxis([-1,1]);
colormap(cmpBlueRed);
xlabel('PC_1 (arb. units)');
ylabel('PC_2 (arb. units)');
zlabel('PC_3 (arb. units)');
cbar = colorbar;
ylabel(cbar, 'L3_{||} (mm)');
ConfAxis('fontSize', 14);

%% Scatter the first three PCs

meanFreq = squeeze(mean(mean(dPhi(:,:,:),2),1));
MakeFigure;
scattern(score(:, 1:3),10,meanFreq, 'filled');
axis('square');
cbar = colorbar();
ylabel(cbar, 'mean frequency (Hz)');
caxis([0 15]);
colormap(cmpRed);
xlabel('PC_1 (arb. units)');
ylabel('PC_2 (arb. units)');
zlabel('PC_3 (arb. units)');
ConfAxis('fontSize', 14);

MakeFigure;
scattern(score(:,1:3), 10, squeeze(L(t==0,8,:)), 'filled');
axis('square');
caxis([-1,1]);
colormap(cmpBlueRed);
xlabel('PC_1 (arb. units)');
ylabel('PC_2 (arb. units)');
zlabel('PC_3 (arb. units)');
cbar = colorbar;
ylabel(cbar, 'L2_{||} (mm)');
ConfAxis('fontSize', 14);

MakeFigure;
scattern(score(:,1:3), 10, squeeze(L(t==0,11,:)), 'filled');
axis('square');
caxis([-1,1]);
colormap(cmpBlueRed);
xlabel('PC_1 (arb. units)');
ylabel('PC_2 (arb. units)');
zlabel('PC_3 (arb. units)');
cbar = colorbar;
ylabel(cbar, 'L3_{||} (mm)');
ConfAxis('fontSize', 14);

%% Plot PC1-2 radius versus mean frequency

r = hypot(score(:,1), score(:,2));

meanFreq = squeeze(mean(mean(dPhi(:,:,:),2),1));

binEdgesF = (0:0.25:25)';
binEdgesR = (0:0.2:20)';
binCentersF = binEdgesF(1:end-1) + diff(binEdgesF)/2;
binCentersR = binEdgesR(1:end-1) + diff(binEdgesR)/2;

f = histcounts2(r, meanFreq, binEdgesR, binEdgesF, 'normalization','pdf');
MakeFigure;
imagesc(binCentersR, binCentersF,f');
hold on;
contour(binCentersR, binCentersF', f', 0:0.01:0.03, 'EdgeColor','k', 'LineWidth', 2);
axis('xy', 'square','tight');
xlabel('PC 1-2 radius (arb. units)');
ylabel('Mean frequency (Hz)');
cbar = colorbar;
ylabel(cbar, 'pdf (1/arb. units/Hz)');
caxis([0 0.03]);
colormap(cmpRed);
ConfAxis('fontSize', 16);

%% Plot PC3-4 radius versus mean frequency

r = hypot(score(:,3), score(:,4));

meanFreq = squeeze(mean(mean(dPhi(:,:,:),2),1));

binEdgesF = (0:0.25:25)';
binEdgesR = (0:0.2:14)';
binCentersF = binEdgesF(1:end-1) + diff(binEdgesF)/2;
binCentersR = binEdgesR(1:end-1) + diff(binEdgesR)/2;

f = histcounts2(r, meanFreq, binEdgesR, binEdgesF, 'normalization','pdf');
MakeFigure;
imagesc(binCentersR, binCentersF,f');
hold on;
contour(binCentersR, binCentersF', f', 0:0.01:0.06, 'EdgeColor','k', 'LineWidth', 2);
axis('xy', 'square','tight');
xlabel('PC 3-4 radius (arb. units)');
ylabel('Mean frequency (Hz)');
cbar = colorbar;
ylabel(cbar, 'pdf (1/arb. units/Hz)');
caxis([0 0.06]);
colormap(cmpRed);
ConfAxis('fontSize', 16);
