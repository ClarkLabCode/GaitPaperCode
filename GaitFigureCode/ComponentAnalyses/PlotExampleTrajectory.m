function singleFly = PlotExampleTrajectory( newData, trigger_type, exp_type, lineThickness, cmap, picked_ID )
% This function splits all the trajectories into to buckets, those that
% stop and those that don't stop and plots these two groups separately in
% two subplots.

%% Cut individual trajectories from the full dataset

if strcmp(trigger_type, 'Optogenetic')
    % % % FOR OPTOGENETIC CASE % % %
    trigger = newData.head_hit;
    
elseif strcmp(trigger_type, 'Visual')
    % % % FOR VISUAL STIMULUS CASE % % %
    stim_dur = 960;
    trigger = newData.stimulusOnset & (newData.stimulusVelocity_x == stim_dur | newData.stimulusVelocity_y == stim_dur);
    
elseif strcmp(trigger_type, 'Control')
    % % % FOR RANDOM TRIGGER CONTROL COMPARISON % % %
    trigger = zeros(height(newData),1);
    if strcmp(exp_type, 'Moonwalker')
        n = 3000; % Moonwalker Opto Stimlus
    elseif strcmp(exp_type, 'Visual')
        n = 2000; % Visual Stimulus
    end
    trigger(1:n) = 1;
    idx = randperm(length(trigger));
    trigger = trigger(idx);

else
    error('Invalid trigger type.');
    
end

% Define all the inputs for cropping trajectories
winlen = 50;
% varList = newData.Properties.VariableNames;
varList = {'forwardSpeed_mmPerSec',...
    'L1_StanceDur_mmPerSec', 'L2_StanceDur_mmPerSec', 'L3_StanceDur_mmPerSec', ...
    'R1_StanceDur_mmPerSec', 'R2_StanceDur_mmPerSec', 'R3_StanceDur_mmPerSec',...
    'L1_StanceStart', 'L2_StanceStart', 'L3_StanceStart',...
    'R1_StanceStart', 'R2_StanceStart', 'R3_StanceStart',...
    'L1_down_cam', 'L2_down_cam', 'L3_down_cam', ...
    'R1_down_cam', 'R2_down_cam', 'R3_down_cam', ...
    'L1_yPlot_mm', 'L2_yPlot_mm', 'L3_yPlot_mm', ...
    'R1_yPlot_mm', 'R2_yPlot_mm', 'R3_yPlot_mm', ...
    'L1_xPlot_mm', 'L2_xPlot_mm', 'L3_xPlot_mm', ...
    'R1_xPlot_mm', 'R2_xPlot_mm', 'R3_xPlot_mm', ...
    'translationalSpeed_mmPerSec', 'angVel_radPerSec', 'linVel_mmPerSec'};

% Extract the individual trajectories in a cell array 
[ A, ~ ] = computeEventTriggeredAverages( newData, trigger, winlen, varList );

%% Remove all trajectories where the fly is walking below a threshold walking speed pre-stimulus

pre_speed = A{1}(1:winlen-1,:);
removeIdx = mean(pre_speed,1) <= 5; 

for n = 1:length(A)
    A{n}(:,removeIdx) = [];
end

%% Remove all trajectories where the fly is fully stopped in the last frames of the trajectory

% Selection Criteria
post_speed = A{1}(winlen+2:end,:);
removeIdx = any(post_speed <= 0,1);

for n = 1:length(A)
    A{n}(:,removeIdx) = [];
end

%% Gather the data from the example fly into a single data table

singleFlyVarList = {'forwardSpeed_mmPerSec',...
    'L1_yPlot_mm', 'L2_yPlot_mm', 'L3_yPlot_mm', ...
    'R1_yPlot_mm', 'R2_yPlot_mm', 'R3_yPlot_mm', ...
    'L1_xPlot_mm', 'L2_xPlot_mm', 'L3_xPlot_mm', ...
    'R1_xPlot_mm', 'R2_xPlot_mm', 'R3_xPlot_mm', ...
    'translationalSpeed_mmPerSec', 'angVel_radPerSec','linVel_mmPerSec'};
locs = [1,20:34];

for idx = 1:length(locs)
    X(:,idx) = A{locs(idx)}(:,picked_ID);
end

singleFly = array2table(X,'VariableNames',singleFlyVarList);

%% Define the selected trajectory

n = picked_ID;

% Define the time variable
time = (([-winlen:winlen])/150)*1000; % in milliseconds. 150 fps is the frame rate
    
% Create a figure
makeFigure;

% Plot the y-position data below
subplot(3,1,1);
hold on;
l1 = plot(time, A{20}(:,n), 'linestyle', '-', 'LineWidth', lineThickness, 'marker', 'none', 'color', cmap(1,:));
l2 = plot(time, A{21}(:,n), 'linestyle', '-', 'LineWidth', lineThickness, 'marker', 'none', 'color', cmap(2,:));
l3 = plot(time, A{22}(:,n), 'linestyle', '-', 'LineWidth', lineThickness, 'marker', 'none', 'color', cmap(3,:));
r1 = plot(time, A{23}(:,n), 'linestyle', '-', 'LineWidth', lineThickness, 'marker', 'none', 'color', cmap(4,:));
r2 = plot(time, A{24}(:,n), 'linestyle', '-', 'LineWidth', lineThickness, 'marker', 'none', 'color', cmap(5,:));
r3 = plot(time, A{25}(:,n), 'linestyle', '-', 'LineWidth', lineThickness, 'marker', 'none', 'color', cmap(6,:));
title(num2str(n));

% Configure the axes
xlim([time(1) time(end)]);
ylim([-3 3]);
% set(gca, 'ytick', []);
set(gca,'Ydir','reverse')
xlabel('Time (ms)');
ylabel('Position (mm)');
PlotConstLine(0,2);
% Add a legend
hleglines = [l1(1), l2(1), l3(1), r1(1), r2(1), r3(1)];
legend(hleglines, {'L1','L2','L3','R1','R2','R3'});
% Configure the axes
ConfAxis;

% Plot the body velocity components
subplot(3,1,2);
hold on;
plot(time, A{1}(:,n), 'linestyle', '-', 'LineWidth', lineThickness, 'marker', 'none', 'color', 'b');
xlabel('Time (msec)');
ylabel({'Body Speed', '(mm/s)'});
xlim([time(1) time(end)]);
ConfAxis;
% Add a line indicating when the stimulus turns on
PlotConstLine(0,2);
PlotConstLine(0,1);

% Plot the step plot below (down_cam)
subplot(3,1,3);
step = true(6,winlen*2+1);
% Populate the row with the current data
step(1,:) = A{19}(:,n)';
step(2,:) = A{18}(:,n)';
step(3,:) = A{17}(:,n)';
step(4,:) = A{16}(:,n)';
step(5,:) = A{15}(:,n)';
step(6,:) = A{14}(:,n)';
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

end