
[ cmpBlueRed, cmpRed, cmpBlue ] = MakeTurningColormaps();

centroidLimits = [-8 8];

viewport = [-30 20];

cmpFeetDown = 0.8*jet(7);
cmpFeetDown = cmpFeetDown(2:7,:);

%% Scatter centroid embedding by centroid velocity components at time 0

MakeFigure;
scattern(centroidUMAP, 10, squeeze(V(t==0,1,:)), 'filled');
axis('equal');
caxis([-500 500]);
colormap(cmpBlueRed);
cbar = colorbar();
ylabel(cbar, 'v_{r}(0) (\circ/s)');
xlabel('z_1 (arb. units)');
ylabel('z_2 (arb. units)');
ConfAxis('fontSize', 14);
xlim(centroidLimits);
ylim(centroidLimits);

MakeFigure;
scattern(centroidUMAP, 10, squeeze(V(t==0,2,:)), 'filled');
axis('equal');
caxis([0 30]);
colormap(cmpRed);
xlabel('z_1 (arb. units)');
ylabel('z_2 (arb. units)');
cbar = colorbar;
ylabel(cbar, 'v_{||}(0) (mm/s)');
ConfAxis('fontSize', 14);
xlim(centroidLimits);
ylim(centroidLimits);

MakeFigure;
scattern(centroidUMAP, 10, squeeze(V(t==0,3,:)), 'filled');
axis('equal');
caxis([-10 10]);
colormap(cmpBlueRed);
xlabel('z_1 (arb. units)');
ylabel('z_2 (arb. units)');
cbar = colorbar;
ylabel(cbar, 'v_{\perp}(0) (mm/s)');
ConfAxis('fontSize', 14);
xlim(centroidLimits);
ylim(centroidLimits);

MakeFigure;
scattern(centroidUMAP, 10, log10(curvature), 'filled');
axis('equal');
caxis([-3 -1]);
colormap(cmpRed);
xlabel('z_1 (arb. units)');
ylabel('z_2 (arb. units)');
cbar = colorbar;
cbar.Ticks = -5:1;
cbar.TickLabels = num2str((-5:1)', '10^{%d}');
ylabel(cbar, 'path curvature (1/mm)');
ConfAxis('fontSize', 14);
xlim(centroidLimits);
ylim(centroidLimits);

%% Label centriod embedding by whether a point is a turn

MakeFigure;
scattern(centroidUMAP(isTurn,:), 20, isExtremum(isTurn), 'filled');
axis('equal');
colormap([0 0 0; 1 0 0]);
xlim(centroidLimits);
ylim(centroidLimits);
cbar = colorbar;
caxis([-1 1]);
cbar.Ticks = [-1/2 1/2];
cbar.TickLabels = {'leftward yaw extremum','rightward yaw extremum'};
xlabel('z_1 (arb. units)');
ylabel('z_2 (arb. units)');
ConfAxis('fontSize', 14);

%% Scatter centroid embedding colored by limb position


for ind = 1:6
    MakeFigure;
    scattern(centroidUMAP, 10, squeeze(Lcentered(t==0,ind,:)), 'filled');
    axis('equal');
    caxis([-0.5 0.5]);
    xlabel('z_1 (arb. units)');
    ylabel('z_2 (arb. units)');
    title('centroid kinematics');
    cbar = colorbar;
    ylabel(cbar, sprintf('%s_{\\perp}(0) (mm)', limbList{ind}));
    cbar.Ticks = -1:0.5:1;
    
    colormap(cmpBlueRed);
    ConfAxis('fontSize', 14);

end

for ind = 1:6
    MakeFigure;
    scattern(centroidUMAP, 10, squeeze(Lcentered(t==0,6+ind,:)), 'filled');
    axis('equal');
    caxis([-1 1]);
    xlabel('z_1 (arb. units)');
    ylabel('z_2 (arb. units)');  
    cbar = colorbar;
    ylabel(cbar, sprintf('%s_{||}(0) (mm)', limbList{ind}));
    cbar.Ticks = -1:0.5:1;
    
    colormap(cmpBlueRed);
    title('centroid kinematics');
    ConfAxis('fontSize', 14);

end

%% Scatter centroid embedding by the number of feet down

MakeFigure;
scattern(centroidUMAP, 10, squeeze(sum(D(t==0,:,:),2)), 'filled');
axis('equal');

cbar = colorbar();
ylabel(cbar, 'total feet down at time 0');
caxis([1 6]);
cbar.Ticks = (1+(6-1)/(2*6)):((6-1)/6):6;
cbar.TickLabels = 1:6;

colormap(cmpFeetDown);
xlabel('z_1 (arb. units)');
ylabel('z_2 (arb. units)');
zlabel('z_3 (arb. units)');
ConfAxis('fontSize', 14);

%% Color limb embedding by distance from camera

% c = [  -23.9313  -43.9128   20.9979];
c = [  -16.2361  -46.4646   23.5763];
% c = [ -51.279346112382285  29.227930036008868  34.199211639809512];

d = sum((limbUMAP(:,1:3) - c).^2,2);
d = (d-min(d))/(max(d)-min(d));
MakeFigure;
scattern(limbUMAP(:,1:3), 10, d, 'filled');
axis('equal');
m = 4;
view(viewport);
xticks(-10:2:10);
yticks(-10:2:10);
zticks(-10:2:10);

colormap(cmpBlue);
xl = xlim; yl = ylim; zl = zlim;

%% Scatter limb embedding by centroid velocity components at time 0 (assumes 3D)

MakeFigure;
scattern(limbUMAP(:, 1:3), 10, squeeze(V(t==0,1,:)), 'filled');
axis('equal');
caxis([-500 500]);
colormap(cmpBlueRed);
cbar = colorbar();
ylabel(cbar, 'v_{r}(0) (\circ s^{-1})');
xlabel('u_1 (arb. units)');
ylabel('u_2 (arb. units)');
zlabel('u_3 (arb. units)');
ConfAxis('fontSize', 14);

view(viewport);
xticks(-10:2:10);
yticks(-10:2:10);
zticks(-10:2:10);
xlim(xl); ylim(yl); zlim(zl);

MakeFigure;
scattern(limbUMAP(:, 1:3), 10, squeeze(V(t==0,2,:)), 'filled');
axis('equal');
caxis([0 30]);
colormap(cmpRed);
xlabel('u_1 (arb. units)');
ylabel('u_2 (arb. units)');
zlabel('u_3 (arb. units)');
cbar = colorbar;
ylabel(cbar, 'v_{||}(0) (mm s^{-1})');
ConfAxis('fontSize', 14);

view(viewport);
xticks(-10:2:10);
yticks(-10:2:10);
zticks(-10:2:10);

MakeFigure;
scattern(limbUMAP(:, 1:3), 10, squeeze(V(t==0,3,:)), 'filled');
axis('equal');
caxis([-10 10]);
colormap(cmpBlueRed);
xlabel('u_1 (arb. units)');
ylabel('u_2 (arb. units)');
zlabel('u_3 (arb. units)');
cbar = colorbar;
ylabel(cbar, 'v_{\perp}(0) (mm s^{-1})');
ConfAxis('fontSize', 14);

view(viewport);
xticks(-10:2:10);
yticks(-10:2:10);
zticks(-10:2:10);

MakeFigure;
scattern(limbUMAP(:, 1:3), 10, log10(curvature), 'filled');
axis('equal');
caxis([-3 -1]);
colormap(cmpRed);
xlabel('u_1 (arb. units)');
ylabel('u_2 (arb. units)');
zlabel('u_3 (arb. units)');
cbar = colorbar;
cbar.Ticks = -5:1;
cbar.TickLabels = num2str((-5:1)', '10^{%d}');
ylabel(cbar, 'path curvature (1/mm)');
ConfAxis('fontSize', 14);

view(viewport);
xticks(-10:2:10);
yticks(-10:2:10);
zticks(-10:2:10);

%% Scatter non-phase dimensions of limb embedding by centroid movements (4d)

MakeFigure;
scattern(limbUMAP(:, [1,4]), 10, squeeze(V(t==0,1,:)), 'filled');
axis('square');
caxis([-500 500]);
colormap(cmpBlueRed);
cbar = colorbar();
ylabel(cbar, 'v_{r}(0) (\circ s^{-1})');
xlabel('u_1 (arb. units)');
ylabel('u_4 (arb. units)');
ConfAxis('fontSize', 14);

MakeFigure;
scattern(limbUMAP(:, [1,4]), 10, squeeze(V(t==0,2,:)), 'filled');
axis('square');
caxis([0 30]);
colormap(cmpRed);
xlabel('u_1 (arb. units)');
ylabel('u_4 (arb. units)');
cbar = colorbar;
ylabel(cbar, 'v_{||}(0) (mm s^{-1})');
ConfAxis('fontSize', 14);

MakeFigure;
scattern(limbUMAP(:, [1,4]), 10, squeeze(V(t==0,3,:)), 'filled');
axis('square');
caxis([-10 10]);
colormap(cmpBlueRed);
xlabel('u_1 (arb. units)');
ylabel('u_4 (arb. units)');
cbar = colorbar;
ylabel(cbar, 'v_{\perp}(0) (mm s^{-1})');
ConfAxis('fontSize', 14);

MakeFigure;
scattern(limbUMAP(:, [1,4]), 10, log10(curvature), 'filled');
axis('square');
caxis([-3 -1]);
colormap(cmpRed);
xlabel('u_1 (arb. units)');
ylabel('u_4 (arb. units)');
cbar = colorbar;
cbar.Ticks = -5:1;
cbar.TickLabels = num2str((-5:1)', '10^{%d}');
ylabel(cbar, 'path curvature (1/mm)');
ConfAxis('fontSize', 14);

%% Scatter limb embedding colored by limb position

for ind = 1:6
    MakeFigure;
    scattern(limbUMAP(:, 1:3), 10, squeeze(Lcentered(t==0,ind,:)), 'filled');
    axis('equal');
    caxis([-0.5 0.5]);
    xlabel('u_1 (arb. units)');
    ylabel('u_2 (arb. units)');
    zlabel('u_3 (arb. units)');
    
    cbar = colorbar;
    ylabel(cbar, sprintf('%s_{\\perp}(0) (mm)', limbList{ind}));
    ConfAxis('fontSize', 14);
    colormap(cmpBlueRed);
    cbar.Ticks = -1:0.5:1;
    
    view(viewport);
    xticks(-10:2:10);
    yticks(-10:2:10);
    zticks(-10:2:10);
end

for ind = 1:6
    MakeFigure;
    scattern(limbUMAP(:, 1:3), 10, squeeze(Lcentered(t==0,6+ind,:)), 'filled');
    axis('equal');
    caxis([-1 1]);
    xlabel('u_1 (arb. units)');
    ylabel('u_2 (arb. units)');
    zlabel('u_3 (arb. units)');
    
    cbar = colorbar;
    ylabel(cbar, sprintf('%s_{||}(0) (mm)', limbList{ind}));
    ConfAxis('fontSize', 14);
    cbar.Ticks = -1:0.5:1;
    colormap(cmpBlueRed);
    
    view(viewport);
    xticks(-10:2:10);
    yticks(-10:2:10);
    zticks(-10:2:10);
end

%% Scatter limb embedding colored by limb position (non-phase dimensions)

for ind = 1:6
    MakeFigure;
    scattern(limbUMAP(:, [1,4]), 10, squeeze(Lcentered(t==0,ind,:)), 'filled');
    caxis([-0.5 0.5]);
    axis('square');
    xlabel('u_1 (arb. units)');
    ylabel('u_4 (arb. units)');
    
    cbar = colorbar;
    ylabel(cbar, sprintf('%s_{\\perp}(0) (mm)', limbList{ind}));
    ConfAxis('fontSize', 14);
    colormap(cmpBlueRed);
    cbar.Ticks = -1:0.5:1;
end

for ind = 1:6
    MakeFigure;
    scattern(limbUMAP(:, [1,4]), 10, squeeze(Lcentered(t==0,6+ind,:)), 'filled');
    caxis([-1 1]);
    axis('square');
    xlabel('u_1 (arb. units)');
    ylabel('u_4 (arb. units)');
    
    cbar = colorbar;
    ylabel(cbar, sprintf('%s_{||}(0) (mm)', limbList{ind}));
    ConfAxis('fontSize', 14);
    cbar.Ticks = -1:0.5:1;
    colormap(cmpBlueRed);
end

%% Scatter limb embedding colored by limb phase

% Make a circular colormap
hmap(1:2^16,1) = linspace(0,1,2^16);
hmap(:,[2 3]) = 0.8;
huemap = hsv2rgb(hmap);

for ind = 1:6
    MakeFigure;
    scattern(limbUMAP(:, 1:3), 10, squeeze(Phi(t==0,ind,:)), 'filled');
    axis('equal');
    caxis([0 1]);
    xlabel('u_1 (arb. units)');
    ylabel('u_2 (arb. units)');
    zlabel('u_3 (arb. units)');
    
    cbar = colorbar;
    ylabel(cbar, sprintf('\\phi_{%s}(0) (cycles)', limbList{ind}));
    ConfAxis('fontSize', 14);
    cbar.Ticks = 0:0.5:1;
    colormap(huemap);
    
    view(viewport);
    xticks(-10:2:10);
    yticks(-10:2:10);
    zticks(-10:2:10);
end

%% Scatter limb emebedding, colored by mean frequency

meanFreq = squeeze(mean(mean(dPhi(:,:,:),2),1));

MakeFigure;
scattern(limbUMAP(:, 1:3),10,meanFreq, 'filled');
axis('equal');

cbar = colorbar();
ylabel(cbar, 'mean frequency (Hz)');

caxis([0 15]);
colormap(cmpRed);
xlabel('u_1 (arb. units)');
ylabel('u_2 (arb. units)');
if size(limbUMAP,2) ==3, zlabel('u_3 (arb. units)'); end;
ConfAxis('fontSize', 14);

view(viewport);
xticks(-10:2:10);
yticks(-10:2:10);
zticks(-10:2:10);

%% Scatter limb embedding, colored by total feet down

MakeFigure;
scattern(limbUMAP(:,1:3), 10, squeeze(sum(D(t==0,:,:),2)), 'filled');
axis('equal');

cbar = colorbar();
ylabel(cbar, 'total feet down at time 0');
caxis([1 6]);
cbar.Ticks = (1+(6-1)/(2*6)):((6-1)/6):6;
cbar.TickLabels = 1:6;


colormap(cmpFeetDown);
xlabel('u_1 (arb. units)');
ylabel('u_2 (arb. units)');
if size(limbUMAP,2) ==3, zlabel('u_3 (arb. units)'); end
ConfAxis('fontSize', 14);

view(viewport);
xticks(-10:2:10);
yticks(-10:2:10);
zticks(-10:2:10);
