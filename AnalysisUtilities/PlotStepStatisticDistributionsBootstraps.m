function [ swing_counts, stance_counts ] = PlotStepStatisticDistributionsBootstraps( newData, num_bootstraps, vfedges, corder_all, f, gof )
% This function plots various statistics as a function of forward speed.
% NOTE: This function assumes that
% computeSwingStanceAmplitudesAndDurations.m has already been run on the
% dataset.

if ~exist('corder_all','var') || isempty(corder_all)
    corder_all = linspecer(6);
end

%% Get the data in order

% Define a limb list
limbList = {'L1','L2','L3','R1','R2','R3'};

% Define some variables for gathering the data
swingVelVarList = {'_SwingAvgVf_mmPerSec'};
stanceVelVarList = {'_StanceAvgVf_mmPerSec'};

% Define some logical variables for isolating events
swingVarList = strcat(limbList, '_SwingStart');
stanceVarList = strcat(limbList, '_StanceStart');

swingStarts = newData{:,swingVarList};
stanceStarts = newData{:,stanceVarList};

%% Partition newData into individual videos for bootstrap analysis

% Define the number of videos
num_videos = length(unique(newData.videoID));

tic;
% Loop through all the videos
for n = 1:num_videos
   cur_data{n} = newData(newData.videoID == n,:);
end
toc;

%% Perform bootstrap analyses (Swing analyses)

% Define the relevant variables
calcVarList = {'_CamStepAmplitude_mm', '_EgoStepAmplitude_mm', '_SwingDur_mmPerSec'};

% Pre-allocate the variable for storing the bootstrap confidence values
ci_swings = nan(6,2,length(vfedges)-1,length(calcVarList));
ci_swings_median = nan(6,2,length(vfedges)-1,length(calcVarList));

tic;
% Loop through each of the limbs
for a = 1:6
    
    % Compute the input data for the current limb and variable combination
    for n = 1:num_videos
        % Define the variables used for this combination of limb and variable type
        curVar = strcat(limbList{a},calcVarList);
        curBoolVar = strcat(limbList{a},'_SwingStart');
        
        % Compute the variable that selects the rows
        curBool = cur_data{n}{:,[curBoolVar]};
        
        % Select out the rows of interest
        temp_steps = cur_data{n}{curBool, [strcat(limbList(a),swingVelVarList), curVar]};
        
        % Filter out NaN values
        temp_steps(isnan(temp_steps(:,2)),:) = [];
        
        % Discretize velocity into bins
        G = discretize(temp_steps(:,1),vfedges);
        
        % Remove the entries that do not fall into any of the bins
        temp_steps(isnan(G),:) = [];
        G(isnan(G)) = [];
        
        % Put the data in a variable for use with bootci
        cur_steps{n} = temp_steps(:,2:end);
        cur_groups{n} = G;
    end
    
    % Convert from columns to rows
    cur_steps = cur_steps';
    cur_groups = cur_groups';
    
    % Define the bootstrapping function
    bootfun = @(a,b)(splitapply(@(x)nanmean(x), cell2mat(a), cell2mat(b)));
    
    % Perform the bootstrapping with
    ci_swings(a,:,:,:) = bootci(num_bootstraps, bootfun, cur_steps, cur_groups);
    
    clear cur_steps cur_groups
    
end
toc;

% Move around the data in ci_swings so that it is organized by analysis
lower_ci_swings = squeeze(ci_swings(:,1,:,:));
upper_ci_swings = squeeze(ci_swings(:,2,:,:));

%% Perform bootstrap analyses (Stance analyses)

% Define the relevant variables
calcVarList = {'_StanceDur_mmPerSec', '_Period_mmPerSec'};

% Pre-allocate the variable for storing the bootstrap confidence values
ci_stances = nan(6,2,length(vfedges)-1,length(calcVarList));

tic;
% Loop through each of the limbs
for a = 1:6
    
    % Compute the input data for the current limb and variable combination
    for n = 1:num_videos
        % Define the variables used for this combination of limb and variable type
        curVar = strcat(limbList{a},calcVarList);
        curBoolVar = strcat(limbList{a},'_StanceStart');
        
        % Compute the variable that selects the rows
        curBool = cur_data{n}{:,[curBoolVar]};
        
        % Select out the rows of interest
        temp_steps = cur_data{n}{curBool, [strcat(limbList(a),stanceVelVarList), curVar]};
        
        % Filter out NaN values
        temp_steps(isnan(temp_steps(:,2)),:) = [];
        
        % Discretize velocity into bins
        G = discretize(temp_steps(:,1),vfedges);
        
        % Remove the entries that do not fall into any of the bins
        temp_steps(isnan(G),:) = [];
        G(isnan(G)) = [];
        
        % Put the data in a variable for use with bootci
        cur_steps{n} = temp_steps(:,2:end);
        cur_groups{n} = G;
    end
    
    % Convert from columns to rows
    cur_steps = cur_steps';
    cur_groups = cur_groups';
    
    % Define the bootstrapping function
    bootfun = @(a,b)(splitapply(@(x)nanmean(x), cell2mat(a), cell2mat(b)));
    
    % Perform the bootstrapping with
    ci_stances(a,:,:,:) = bootci(num_bootstraps, bootfun, cur_steps, cur_groups);
    
    clear cur_steps cur_groups
    
end
toc;

% Move around the data in ci_swings so that it is organized by analysis
lower_ci_stances = squeeze(ci_stances(:,1,:,:));
upper_ci_stances = squeeze(ci_stances(:,2,:,:));

%% Create a plot for the Camera Step Amplitude

% Calculate the camera frame step amplitudes
[steps_mean, ~, ~, ~] = calculateStepStatistics(newData, limbList, ...
    swingVelVarList, '_CamStepAmplitude_mm', swingStarts, vfedges);
    
% Plot the camera frame step amplitudes
vfBins = repmat((vfedges(1:end-1)+vfedges(2:end))'./2, 1,6);

MakeFigure;
PlotConfidenceIntervalWithErrorPatch( vfBins, steps_mean, lower_ci_swings(:,:,1)', upper_ci_swings(:,:,1)', corder_all);
legend(limbList);
ylabel({'Camera Frame Step Amplitude','mm'});
ylim([0 2.5]);
xlim([0 35]);
xlabel({'Forward Speed','mm/s'});
ConfAxis;

% %% Create a plot for the Ego Step Amplitude
% 
% % Calculate the camera frame step amplitudes
% [steps_mean, ~, ~] = calculateStepStatistics(newData, limbList, ...
%     velVarList, '_EgoStepAmplitude_mm', swingStarts, vfedges);
%     
% % Plot the camera frame step amplitudes
% vfBins = repmat((vfedges(1:end-1)+vfedges(2:end))'./2, 1,6);
% 
% MakeFigure;
% PlotConfidenceIntervalWithErrorPatch( vfBins, steps_mean, lower_ci_swings(:,:,2)', upper_ci_swings(:,:,2)', corder_all);
% legend(limbList);
% ylabel({'Egocentric Step Amplitude','mm'});
% ylim([0 1.5]);
% xlabel({'Forward Speed','mm/s'});
% ConfAxis;

%% Create a plot for the Swing Duration

% Calculate the camera frame step amplitudes
[steps_mean, ~, steps_count, ~] = calculateStepStatistics(newData, limbList, ...
    swingVelVarList, '_SwingDur_mmPerSec', swingStarts, vfedges);
    
% Calculate the total number of swings for each limb
swing_counts = sum(steps_count)';

% Plot the camera frame step amplitudes
vfBins = repmat((vfedges(1:end-1)+vfedges(2:end))'./2, 1,6);

MakeFigure;
PlotConfidenceIntervalWithErrorPatch( vfBins, steps_mean, lower_ci_swings(:,:,3)', upper_ci_swings(:,:,3)', corder_all);
legend(limbList);
ylabel({'Swing Duration','ms'});
ylim([0 50]);
xlim([0 35]);
xlabel({'Forward Speed','mm/s'});
ConfAxis;

%% Create a plot for the Stance Duration

% Calculate the camera frame step amplitudes
[steps_mean, ~, steps_count, ~] = calculateStepStatistics(newData, limbList, ...
    stanceVelVarList, '_StanceDur_mmPerSec', stanceStarts, vfedges);
    
% Calculate the total number of stances for each limb
stance_counts = sum(steps_count)';

% Plot the camera frame step amplitudes
vfBins = repmat((vfedges(1:end-1)+vfedges(2:end))'./2, 1,6);

MakeFigure;
PlotConfidenceIntervalWithErrorPatch( vfBins, steps_mean, lower_ci_stances(:,:,1)', upper_ci_stances(:,:,1)', corder_all);
plot(vfBins(:,1), f.a .* vfBins(:,1).^f.b, 'color', [.5,.5,.5], 'linestyle', '--', 'linewidth', 1);
legend(limbList);
ylabel({'Stance Duration','ms'});
ylim([0 250]);
xlim([0 35]);
xlabel({'Forward Speed','mm/s'});
title(strcat('r^2=', num2str(gof.rsquare)));
ConfAxis;

% %% Create a plot for the Stepping Period
% 
% % Calculate the camera frame step amplitudes
% [steps_mean, ~, ~] = calculateStepStatistics(newData, limbList, ...
%     velVarList, '_Period_mmPerSec', stanceStarts, vfedges);
%     
% % Plot the camera frame step amplitudes
% vfBins = repmat((vfedges(1:end-1)+vfedges(2:end))'./2, 1,6);
% 
% MakeFigure;
% PlotConfidenceIntervalWithErrorPatch( vfBins, steps_mean, lower_ci_stances(:,:,2)', upper_ci_stances(:,:,2)', corder_all);
% legend(limbList);
% ylabel({'Step Period','ms'});
% ylim([0 1000]);
% % ylim([0 800]);
% % ylim([0 1600]);
% xlabel({'Forward Speed','mm/s'});
% ConfAxis;

end

%% Local function for calculating 

function [steps_mean, steps_sem, steps_count, steps_median] = calculateStepStatistics(newData, limbList, ...
    velVarList, var_name, logicals, vfedges)

for i =1:6
    
    % Get the velocity data and associated variable
    steps = newData{logicals(:,i), [strcat(limbList(i), velVarList),strcat(limbList(i),var_name)]};
    
    % Filter out NaN values
    steps(isnan(steps(:,2)),:) = [];
    
    % Discretize velocity into bins
    G = discretize(steps(:,1),vfedges);
    
    % Remove the entries that do not fall into any of the bins
    steps(isnan(G),:) = [];
    G(isnan(G)) = [];
    
    % Bin by walking speed to get the mean and standard deviation    
    steps_mean(:,i) = splitapply(@mean,steps(:,2), G);
    steps_median(:,i) = splitapply(@median,steps(:,2), G);
    steps_sem(:,i) = splitapply(@(x) std(x)/(length(x)^(1/2)),steps(:,2), G);
    steps_count(:,i) = splitapply(@(x) length(x),steps(:,2), G);
     
end

end