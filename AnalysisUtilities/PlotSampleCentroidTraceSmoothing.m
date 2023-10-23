function PlotSampleCentroidTraceSmoothing( newData )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

mmPerPix = 0.043;
fps = 150;

winlen = 30;

% Extract data
[ ID, V, ~, ~, CM ] = ExtractFoldedTurningData( newData, winlen, true, true);

%% Plot random sample

vind = find(squeeze(nanmean(V(:,2,:),1))>15);
p = randperm(size(vind,1),1);
idx = vind(p(1));

cmsmoothed = CM(:,:,idx);
v = V(:,:,idx);
id = ID(:,:,idx);

tableIdx = ismember(newData{:, {'uniqueFlyTrajID','Frame'}}, id(:, 1:2), 'rows');
cmraw = newData{tableIdx, {'xCOM','yCOM'}};
theta = newData{tableIdx, {'Orient_Rad_FullRotation','smooth_Orient_Rad_FullRotation'}};

theta = pi/2-theta;

thetaraw = theta(:,1);
thetasmoothed = theta(:,2);

arrowInterval = 3;

figure('Position',[200,500,500,700],'WindowStyle','docked'); hold on;
plot(cmraw(:,1), cmraw(:,2), '-k', 'linewidth', 2);
hold on;
plot(cmsmoothed(:,1), cmsmoothed(:,2),'-r', 'linewidth', 2);
axis('equal');


plot([min(xlim), min(xlim) + 1/mmPerPix], [min(ylim), min(ylim)], '-k', 'linewidth', 2);
text(min(xlim) + 0.5/mmPerPix, min(ylim)-0.1/mmPerPix, '1 mm', 'FontSize', 20, 'horizontalalignment','center');

figure('Position',[200,500,500,700],'WindowStyle','docked'); hold on;
plot(cmraw(:,1), cmraw(:,2), '-k', 'linewidth', 2);
hold on;
quiver(cmraw(1:arrowInterval:end,1), cmraw(1:arrowInterval:end,2), cos(thetaraw(1:arrowInterval:end)), sin(thetaraw(1:arrowInterval:end)), 'k', 'linewidth', 1, 'autoscalefactor', 1, 'ShowArrowHead','off');
plot(cmsmoothed(:,1), cmsmoothed(:,2),'-r', 'linewidth', 2);
quiver(cmsmoothed(1:arrowInterval:end,1), cmsmoothed(1:arrowInterval:end,2),  cos(thetasmoothed(1:arrowInterval:end)), sin(thetasmoothed(1:arrowInterval:end)), 'r', 'linewidth', 1, 'ShowArrowHead','off');

axis('equal');

plot([min(xlim), min(xlim) + 1/mmPerPix], [min(ylim), min(ylim)], '-k', 'linewidth', 2);
text(min(xlim) + 0.5/mmPerPix, min(ylim)-0.1/mmPerPix, '1 mm', 'FontSize', 20, 'horizontalalignment','center');

end

