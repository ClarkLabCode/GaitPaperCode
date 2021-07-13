% PlotTetrapodExample.m: This script gathers the data for plotting a
% single tripod example from the dataset

% NOTE: This is pulled from the IsoD1 masked dataset

%% Clip to the desired trajectory

% Define the ID of interest
exampleID = 6097;

% Define the endpoints of the trajectory
frameStart = 555884;
frameEnd = 555984;

% Isolate the trajectory of interest
singleFly = newData(newData.uniqueFlyTrajID == exampleID,:);

% Clip the trajectory to the desired length
idx = (singleFly.Frame > frameStart & singleFly.Frame <= frameEnd);
singleFly = singleFly(idx,:);

%% Plot the result

% % Define the colors used
% colors = linspecer(6);
% colors_all = linspecer(9);

% Get the trajectory length
traj_length = size(singleFly,1);

% Calculate the time variable
time = (([1:traj_length]-1)/150)*1000; % in milliseconds. 150 fps is the frame rate

% Make a new figure
MakeFigure;

% Plot the limb relative phases
subplot(3,1,1);
title('Tetrapod Gait');
hold on;
% plot(time, mod(singleFly.InstantaneousPhase_L1y-singleFly.InstantaneousPhase_R1y, (2*pi))/(2*pi),...
%     'linestyle', '-', 'LineWidth', 2, 'marker', 'none', 'color', colors_all(1,:));
% plot(time, mod(singleFly.InstantaneousPhase_L2y-singleFly.InstantaneousPhase_R2y, (2*pi))/(2*pi),...
%     'linestyle', '-', 'LineWidth', 2, 'marker', 'none', 'color', colors_all(2,:));
% plot(time, mod(singleFly.InstantaneousPhase_L3y-singleFly.InstantaneousPhase_R3y, (2*pi))/(2*pi),...
%     'linestyle', '-', 'LineWidth', 2, 'marker', 'none', 'color', colors_all(3,:));
plot(time, mod(singleFly.smooth_InstantaneousPhase_L1y-singleFly.smooth_InstantaneousPhase_R1y, (2*pi))/(2*pi),...
    'linestyle', '-', 'LineWidth', 2, 'marker', 'none', 'color', colors_all(1,:));
plot(time, mod(singleFly.smooth_InstantaneousPhase_L2y-singleFly.smooth_InstantaneousPhase_R2y, (2*pi))/(2*pi),...
    'linestyle', '-', 'LineWidth', 2, 'marker', 'none', 'color', colors_all(2,:));
plot(time, mod(singleFly.smooth_InstantaneousPhase_L3y-singleFly.smooth_InstantaneousPhase_R3y, (2*pi))/(2*pi),...
    'linestyle', '-', 'LineWidth', 2, 'marker', 'none', 'color', colors_all(3,:));

xlabel('Time (ms)');
ylabel('Relative Phase (cycles)');
legend({'L1-R1','L2-R2','L3-L3'});
ConfAxis;

% Plot the forward body velocity
subplot(3,1,2);
hold on;
% plot(time, singleFly.translationalSpeed, 'linestyle', '-', 'LineWidth', 2, 'marker', 'none', 'color', colors_all(1,:));
plot(time, singleFly.smooth_translationalSpeed_mmPerSec, 'linestyle', '-', 'LineWidth', 2, 'marker', 'none', 'color', colors_all(1,:));
xlabel('Time (ms)');
ylabel('v_{\perp} (mm/s)');
ConfAxis;

% % Plot the y-position data below
% subplot(3,1,1);
% hold on;
% plot(time, singleFly.L1_yPlot_mm, 'linestyle', '-', 'LineWidth', 2, 'marker', 'none', 'color', colors(1,:));
% plot(time, singleFly.L2_yPlot_mm, 'linestyle', '-', 'LineWidth', 2, 'marker', 'none', 'color', colors(2,:));
% plot(time, singleFly.L3_yPlot_mm, 'linestyle', '-', 'LineWidth', 2, 'marker', 'none', 'color', colors(3,:));
% plot(time, singleFly.R1_yPlot_mm, 'linestyle', '-', 'LineWidth', 2, 'marker', 'none', 'color', colors(4,:));
% plot(time, singleFly.R2_yPlot_mm, 'linestyle', '-', 'LineWidth', 2, 'marker', 'none', 'color', colors(5,:));
% plot(time, singleFly.R3_yPlot_mm, 'linestyle', '-', 'LineWidth', 2, 'marker', 'none', 'color', colors(6,:));
% title('Tetrapod Gait');
% 
% % Set the axis limits
% xlim([0 traj_length/150*1000]);
% ylim([-3 3]);
% 
% % Hide the y-axis
% set(gca, 'ytick', []);
% 
% % Invert the y-axis
% set(gca,'Ydir','reverse')
% 
% % Provide an x-axis label
% xlabel('Time (msec)');
% 
% % Provide a y-axis label
% ylabel('Position (mm)');
% 
% % Add a legend
% legend({'L1','L2','L3','R1','R2','R3'});
% 
% % Configure the axes
% ConfAxis;

% Plot the step plot below (down_cam)
subplot(3,1,3);
step = true(6,traj_length);

% Populate the row with the current data
step(1,:) = singleFly.R3_down_cam;
step(2,:) = singleFly.R2_down_cam;
step(3,:) = singleFly.R1_down_cam;
step(4,:) = singleFly.L3_down_cam;
step(5,:) = singleFly.L2_down_cam;
step(6,:) = singleFly.L1_down_cam;

% Set the tick marks on the y-axis
set(gca, 'ytick', [1:6]);
set(gca, 'yticklabels', {'R3','R2','R1','L3','L2','L1'});

% Invert the y-axis
set(gca,'Ydir','reverse');

% Hide the x-axis
set(gca, 'xtick', []);

% Visualize the image
hold on;
imagesc(step);
colormap('gray');

% Provide an x-axis label
xlabel('Time (msec)');

% Configure the axes
axis('tight');
ConfAxis;
set(gca,'linewidth',6);

% % Plot the step plot below (down)
% subplot(3,1,3);
% step = true(6,traj_length);
% 
% % Populate the row with the current data
% step(1,:) = singleFly.R3_down;
% step(2,:) = singleFly.R2_down;
% step(3,:) = singleFly.R1_down;
% step(4,:) = singleFly.L3_down;
% step(5,:) = singleFly.L2_down;
% step(6,:) = singleFly.L1_down;
% 
% % Set the tick marks on the y-axis
% set(gca, 'ytick', [1:6]);
% set(gca, 'yticklabels', {'R3','R2','R1','L3','L2','L1'});
% 
% % Invert the y-axis
% set(gca,'Ydir','reverse');
% 
% % Hide the x-axis
% set(gca, 'xtick', []);
% 
% % Visualize the image
% hold on;
% imagesc(step);
% colormap('gray');
% 
% % Provide an x-axis label
% xlabel('Time (msec)');
% 
% % Configure the axes
% axis('tight');
% ConfAxis;
% set(gca,'linewidth',6);