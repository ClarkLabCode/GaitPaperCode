function calcGaitIndex( newData )
% This function computes the tripod index on each of the example
% trajectories extracted in the gait paper.

% Gait Index
winlen = 8;
% +1 for tripod
% -1 for tetrapod
% 0 for non-canonical

% Tripod Index - Percentage of frames in a video that display leg
% combinations defined by the tripod gait

% Tetrapod Index - Percentage of frames in a video that display leg
% combinations defined by the tetrapod gait


%% Tripod-like example

% Define the ID of interest
exampleID = 584;

% Define the endpoints of the trajectory
frameStart = 356136;
frameEnd = 356236;

% Get the desired data
[ singleFly ] = getData(newData, exampleID, frameStart, frameEnd, winlen);

% Plot the fly
plotExampleFly(singleFly, 'Tripod'); 


%% Tetrapod-like example

% Define the ID of interest
exampleID = 6097;

% Define the endpoints of the trajectory
frameStart = 555884;
frameEnd = 555984;

% Get the desired data
[ singleFly ] = getData(newData, exampleID, frameStart, frameEnd, winlen);

% Plot the fly
plotExampleFly(singleFly, 'Tetrapod-like'); 


%% Non-canonical example

% Define the ID of interest
exampleID = 827;

% Define the endpoints of the trajectory
frameStart = 475745;
frameEnd = 475845;

% Get the desired data
[ singleFly ] = getData(newData, exampleID, frameStart, frameEnd, winlen);

% Plot the fly
plotExampleFly(singleFly, 'Non-canonical'); 

end

function [ score ] = getScore( numFeet )
% Convert number of feet down into the gait index value

score = zeros(size(numFeet));

score(numFeet == 3) = 1;
score(numFeet == 4) = -1;

end

function [ triScore, tetraScore ] = getTriTetraPercentages( numFeet )

triScore = zeros(size(numFeet));
tetraScore = zeros(size(numFeet));

triScore(numFeet == 3) = 1;
tetraScore(numFeet == 4) = 1;

end

function [ singleFly ] = getData(newData, exampleID, frameStart, frameEnd, winlen)

% Isolate the trajectory of interest
singleFly = newData(newData.uniqueFlyTrajID == exampleID,:);

% Compute the number of feet down
varList = {'L1_down_cam','L2_down_cam','L3_down_cam',...
           'R1_down_cam','R2_down_cam','R3_down_cam'};
singleFly.numFeet = sum(singleFly{:,varList},2);

% Score the number of feet in the trajectory
[ singleFly.score ] = getScore( singleFly.numFeet );

% Score the percentage tripod and percentage tetrapod
[ singleFly.triScore, singleFly.tetraScore ] = getTriTetraPercentages( singleFly.numFeet );

% Compute the moving tripod and tetrapod percentages over the desired window
singleFly.TripodIndex = movmean(singleFly.triScore, winlen);
singleFly.TetrapodIndex = movmean(singleFly.tetraScore, winlen);

% Compute the moving gait index over the desired window
singleFly.GaitIndex = movmean(singleFly.score, winlen);

% Clip the trajectory to the desired length
idx = (singleFly.Frame > frameStart & singleFly.Frame <= frameEnd);
singleFly = singleFly(idx,:);


end

function plotExampleFly(singleFly, titleText)

% Define the colormap
colors = linspecer(6);

% Get the trajectory length
traj_length = size(singleFly,1);

% Calculate the time variable
time = (([1:traj_length]-1)/150)*1000; % in milliseconds. 150 fps is the frame rate

% Make a new figure
MakeFigure;

% Plot the y-position data below
subplot(4,1,1);
hold on;
plot(time, singleFly.L1_yPlot_mm, 'linestyle', '-', 'LineWidth', 2, 'marker', 'none', 'color', colors(1,:));
plot(time, singleFly.L2_yPlot_mm, 'linestyle', '-', 'LineWidth', 2, 'marker', 'none', 'color', colors(2,:));
plot(time, singleFly.L3_yPlot_mm, 'linestyle', '-', 'LineWidth', 2, 'marker', 'none', 'color', colors(3,:));
plot(time, singleFly.R1_yPlot_mm, 'linestyle', '-', 'LineWidth', 2, 'marker', 'none', 'color', colors(4,:));
plot(time, singleFly.R2_yPlot_mm, 'linestyle', '-', 'LineWidth', 2, 'marker', 'none', 'color', colors(5,:));
plot(time, singleFly.R3_yPlot_mm, 'linestyle', '-', 'LineWidth', 2, 'marker', 'none', 'color', colors(6,:));
title(titleText);

% Set the axis limits
xlim([0 traj_length/150*1000]);
ylim([-3 3]);

% % Hide the y-axis
% set(gca, 'ytick', []);

% Invert the y-axis
set(gca,'Ydir','reverse')

% Provide an x-axis label
xlabel('Time (ms)');

% Provide a y-axis label
ylabel('Position (mm)');

% Add a legend
legend({'L1','L2','L3','R1','R2','R3'});

% Configure the axes
ConfAxis;

% Plot the step plot below (down_cam)
subplot(4,1,2);
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
xlabel('Time (ms)');

% Configure the axes
axis('tight');
ConfAxis;
set(gca,'linewidth',6);

% Plot the gait index
subplot(4,1,3);
plot(time, singleFly.GaitIndex, 'linestyle', '-', 'LineWidth', 2, 'marker', 'none', 'color', 'k');
% Set the axis limits
xlim([0 traj_length/150*1000]);
ylim([-1 1]);
% Provide an x-axis label
xlabel('Time (ms)');
% Provide a y-axis label
ylabel('Gait Index');
% Configure the axes
ConfAxis;

% Plot the tripod and tetrapod percentages
subplot(4,1,4);
hold on;
plot(time, singleFly.TripodIndex, 'linestyle', '-', 'LineWidth', 2, 'marker', 'none', 'color', colors(1,:));
plot(time, singleFly.TetrapodIndex, 'linestyle', '-', 'LineWidth', 2, 'marker', 'none', 'color', colors(2,:));
% Add a legend
legend({'Tripod Index', 'Tetrapod Index'});
% Set the axis limits
xlim([0 traj_length/150*1000]);
ylim([0 1]);
% Provide an x-axis label
xlabel('Time (ms)');
% Provide a y-axis label
ylabel('Fraction of Frames');
% Configure the axes
ConfAxis;


end