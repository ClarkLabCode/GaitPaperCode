function [ p, h, stats, data, time, meanSpeed, walk_ci, normMeanSpeed, normWalk_ci, postDurations_pdf, post_pdf_ci, postDurations_cdf, post_cdf_ci, post_stances ] = PlotOptoStanceIncreaseAllLimbs( newData, trigger_type, exp_type, lineThickness, suppressPlots, num_bootstraps )
% This function splits all the trajectories into to buckets, those that
% stop and those that don't stop and plots these two groups separately in
% two subplots.

%% Cut individual trajectories from the full dataset

if strcmp(trigger_type, 'Optogenetic')
    % % % FOR OPTOGENETIC CASE % % %
    trigger = newData.head_hit;
    
elseif strcmp(trigger_type, 'Visual')
    % % % FOR VISUAL STIMULUS CASE % % %
    stim_dur = [720, 960, 1440];
    trigger = newData.stimulusOnset & ...
        (newData.stimulusVelocity_x == stim_dur(1) | newData.stimulusVelocity_y == stim_dur(1) | ...
         newData.stimulusVelocity_x == stim_dur(2) | newData.stimulusVelocity_y == stim_dur(2) | ...     
         newData.stimulusVelocity_x == stim_dur(3) | newData.stimulusVelocity_y == stim_dur(3));
    
elseif strcmp(trigger_type, 'Control')
    % % % FOR RANDOM TRIGGER CONTROL COMPARISON % % %
    trigger = zeros(height(newData),1);
    if strcmp(exp_type, 'Moonwalker')
        n = 3000; % Moonwalker Opto Stimlus
    elseif strcmp(exp_type, 'Visual')
        n = 4000; % Visual Stimulus
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
    'R1_down_cam', 'R2_down_cam', 'R3_down_cam'};

% Extract the individual trajectories in a cell array 
[ A, ~ ] = computeEventTriggeredAverages( newData, trigger, winlen, varList );

%% Remove all trajectories where the fly is walking below a threshold walking speed pre-stimulus

% Define the range of frames that we use for prestimulus averaging
preStimRange = 36:winlen-1;

pre_speed = A{1}(preStimRange,:);
removeIdx = mean(pre_speed,1) <= 5; 

for n = 1:length(A)
    A{n}(:,removeIdx) = [];
end

%% Remove all trajectories where the fly is fully stopped in the last frames of the trajectory

% % Selection Criteria
% post_speed = A{1}(end-20:end,:);
% removeIdx = any(post_speed <= 1,1);

% Selection Criteria
post_speed = A{1}(winlen+2:end,:);
removeIdx = any(post_speed <= 0.1, 1);

for n = 1:length(A)
    A{n}(:,removeIdx) = [];
end

%% Get the pre + post stimulus stance duration for each trajectory

% Define the limb ordering
limbList = {'L1','L2','L3','R1','R2','R3'};

% Get the number of trajectories in this data
num_ex = size(A{1},2);

% Initialize variables for storing the location and duration data (pre + post)
preArray = 1:winlen;
pre_stanceStart = NaN(length(preArray), num_ex);
pre_durs = NaN(length(preArray), num_ex);

postArray = winlen+1:(winlen*2+1);
post_stanceStart = NaN(length(postArray), num_ex);
post_durs = NaN(length(postArray), num_ex);

idx = 1;
for n = 2:7
    
    % Get the last non-NaN stance duration before trigger in each trajectory
    pre_stanceStart(:,:,idx) = logical(A{n+6}(preArray, :)); %L1 --> R3
    pre_durs(:,:,idx) = A{n}(preArray, :);
    
    % Get the first non-NaN stance duration after trigger in each trajectory
    post_stanceStart(:,:,idx) = logical(A{n+6}(postArray, :));
    post_durs(:,:,idx) = A{n}(postArray, :);
    
    % Increment the limb counter
    idx = idx + 1;
    
end

% Initialize the variables for storing the pre and post activation data
pre_stances = nan(num_ex, 6);
post_stances = nan(num_ex, 6);

% Populate our two matrice with the stance durations for each limb in each example
for n = 1:num_ex % NOTE: Each one of the _durs variables is the same size
    for m = 1:6
        
        % Get the pre-stimulus stance duration
        curStances = pre_durs(logical(pre_stanceStart(:,n,m)),n,m);
        try
            pre_stances(n,m) = curStances(end-1);
        catch
            pre_stances(n,m) = NaN;
        end
         
        % Get the post-stimulus stance duration
        curStances = post_durs(logical(post_stanceStart(:,n,m)),n,m);
        try
            post_stances(n,m) = curStances(1);
        catch
            post_stances(n,m) = NaN;
        end
        
    end
end

%% Identify trajectories in which the fly was either not walking prestimulus or stopped post stimulus based on stance durations

% Get an index to select the trajectories that were kept in the analysis
keepExIdx = ~any(isnan([pre_stances, post_stances]),2);

% Remove the nan values from the pre_stance and post_stance data
pre_stances = pre_stances(keepExIdx,:);
post_stances = post_stances(keepExIdx,:);

% Remove any rows that contain entries larger than 250 ms (Suggesting a misidentified step length in this context, given the restrictions on stopping)
keepExIdx = ~any(([pre_stances, post_stances] > 250),2);
pre_stances = pre_stances(keepExIdx,:);
post_stances = post_stances(keepExIdx,:);

% Collapse all the data into two column matrices after removing any rows that correspond to a removed trajectory
data = [pre_stances(:), post_stances(:)];

%% Bootstrap the body velocity and stance durations over trajectories to estimate error

% Define the bin sizes
% dur_binEdges = [-1000/150:2000/150:500]; % 2 frame bins
dur_binEdges = [-500/150:1000/150:250]; % 1 frame bins
dur_binCenters = (dur_binEdges(1:end-1) + dur_binEdges(2:end))/2; 

% Bootstrap the pre-stimulus stance durations
bootfun = @(x) ksdensity(x(:), dur_binCenters);
% pre_pdf_ci = bootci(num_bootstraps, {bootfun, pre_stances}, 'type', 'percentile');
pre_pdf_ci = bootci(num_bootstraps, {bootfun, pre_stances}, 'type', 'bca');

% Boostrap the post-stimulus stance durations
bootfun = @(x) ksdensity(x(:), dur_binCenters);
% post_pdf_ci = bootci(num_bootstraps, {bootfun, post_stances}, 'type', 'percentile');
post_pdf_ci = bootci(num_bootstraps, {bootfun, post_stances}, 'type', 'bca');

% % Bootstrap the pre-stimulus stance durations
bootfun = @(x) ksdensity(x(:), dur_binCenters, 'Function', 'cdf');
% pre_cdf_ci = bootci(num_bootstraps, {bootfun, pre_stances}, 'type', 'percentile');
pre_cdf_ci = bootci(num_bootstraps, {bootfun, pre_stances}, 'type', 'bca');

% Boostrap the post-stimulus stance durations
bootfun = @(x) ksdensity(x(:), dur_binCenters, 'Function', 'cdf');
% post_cdf_ci = bootci(num_bootstraps, {bootfun, post_stances}, 'type', 'percentile');
post_cdf_ci = bootci(num_bootstraps, {bootfun, post_stances}, 'type', 'bca');

% Bootstrap the mean trajectories
forwardSpeed = A{1}(:,keepExIdx);
bootfun = @(x)mean(x,1)';
walk_ci = bootci(num_bootstraps, bootfun, forwardSpeed');

% Bootstrap the normalized mean trajectories
forwardSpeed = A{1}(:,keepExIdx);
bootfun = @(x)mean(x./mean(x(:,preStimRange),2),1);
normWalk_ci = bootci(num_bootstraps, bootfun, forwardSpeed');

%% Test the difference between the pre and post stimulus durations using a sign rank test

% Wilcoxon signed rank test
[p,h,stats] = signrank(data(:,1),data(:,2));

%% Plot the mean walking speed for the population (CI ERROR BARS)

% Number of examples included in the analysis
num_ex = sum(keepExIdx);

% Calculate the mean for the unnormalized traces
meanSpeed = mean(A{1}(:,keepExIdx),2);

% Define the time variable
time = [-winlen:winlen]'*1000/150;

if ~suppressPlots
    
    % Unnormalized Walking Speed
    makeFigure;
    PlotConfidenceIntervalWithErrorPatch( time, meanSpeed, walk_ci(1,:)', walk_ci(2,:)');
    title({trigger_type, strcat('Mean Forward Speed (Bootstrapped) n=',num2str(num_ex))});
    xlabel('Time (ms)');
    ylabel('Forward Speed (mm/s)');
    ylim([-5 25]);
    PlotConstLine(0,2);
    ConfAxis;
    
end

% Calculate the mean for the normalized data
preSpeed = mean(A{1}(preStimRange,keepExIdx),1);
speed = A{1}(:,keepExIdx);
normSpeed = speed ./ repmat(preSpeed,winlen*2+1,1);
normMeanSpeed = mean(normSpeed,2);

if ~suppressPlots
    makeFigure;
    PlotConfidenceIntervalWithErrorPatch( time, normMeanSpeed, normWalk_ci(1,:)', normWalk_ci(2,:)');
    title({trigger_type, strcat('Normalized Mean Forward (Bootstrapped) Speed n=',num2str(num_ex))});
    xlabel('Time (ms)');
    ylabel('Forward Speed (mm/s)');
    ylim([-1 2]);
    PlotConstLine(0,2);
    ConfAxis;
end

%% Plot the mean walking speed for the population (SEM ERROR BARS)

% Number of examples included in the analysis
num_ex = sum(keepExIdx);

% Calculate the mean and SEM for the unnormalized traces
meanSpeed = mean(A{1}(:,keepExIdx),2);
semSpeed = std(A{1}(:,keepExIdx),[],2)./(num_ex^(1/2));

% Define the time variable
time = [-winlen:winlen]'*1000/150;

if ~suppressPlots
    
    % Unnormalized Walking Speed
    makeFigure;
    PlotXvsY(time,meanSpeed,'error',semSpeed);
    title({trigger_type, strcat('Mean Forward Speed (SEM) n=',num2str(num_ex))});
    xlabel('Time (ms)');
    ylabel('Forward Speed (mm/s)');
    ylim([-5 25]);
    PlotConstLine(0,2);
    ConfAxis;
    
end

% Calculate the mean and SEM for the normalized data
preSpeed = mean(A{1}(preStimRange,keepExIdx),1);
speed = A{1}(:,keepExIdx);
normSpeed = speed ./ repmat(preSpeed,winlen*2+1,1);
normMeanSpeed = mean(normSpeed,2);
normSemSpeed = std(normSpeed,[],2)./(num_ex^(1/2));

if ~suppressPlots
    makeFigure;
    PlotXvsY(time,normMeanSpeed,'error',normSemSpeed);
    title({trigger_type, strcat('Normalized Mean Forward Speed (SEM) n=',num2str(num_ex))});
    xlabel('Time (ms)');
    ylabel('Forward Speed (mm/s)');
    ylim([-1 2]);
    PlotConstLine(0,2);
    ConfAxis;
end

%% Create a histogram plot for each of the pre and post stimulus stances

% Get the pdf
preDurations_pdf = ksdensity(pre_stances(:), dur_binCenters);
postDurations_pdf = ksdensity(post_stances(:), dur_binCenters);

if ~suppressPlots
    % Plot
    makeFigure;
    hold on;
    PlotConfidenceIntervalWithErrorPatch( [dur_binCenters;dur_binCenters]', [preDurations_pdf;postDurations_pdf]', [pre_pdf_ci(1,:);post_pdf_ci(1,:)]', [pre_pdf_ci(2,:);post_pdf_ci(2,:)]');
    title({trigger_type, strcat('Stance Duration Distributions - All Limbs n=',num2str(num_ex))});
    xlabel('Stance Duration');
    ylabel('pdf');
%     ylim([0 .02]);
    legend({'Pre Stimulus','Post Stimulus'});
    ConfAxis;
end

%% Create a cdf plot for each of the pre and post stimulus stances

% Get the cdf
preDurations_cdf = ksdensity(pre_stances(:), dur_binCenters, 'Function', 'cdf');
postDurations_cdf = ksdensity(post_stances(:), dur_binCenters, 'Function', 'cdf');

if ~suppressPlots
    % Plot
    makeFigure;
    hold on;
    PlotConfidenceIntervalWithErrorPatch([dur_binCenters;dur_binCenters]', [preDurations_cdf;postDurations_cdf]', [pre_cdf_ci(1,:);post_cdf_ci(1,:)]', [pre_cdf_ci(2,:);post_cdf_ci(2,:)]');
    title({trigger_type, strcat('Stance Duration Distributions - All Limbs n=',num2str(num_ex))});
    xlabel('Stance Duration');
    ylabel('cdf');
    ylim([0 1]);
    legend({'Pre Stimulus','Post Stimulus'});
    ConfAxis;
end

%% Create a scatter plot of all datapoints

if ~suppressPlots
    makeFigure;
    s = scatter(data(:,1), data(:,2), 'filled');
    s.MarkerFaceAlpha = .2;
    hold on;
    plot([1:1000],[1:1000],'linewidth', 2, 'color', [.5,.5,.5], 'linestyle', '--');
    xlabel('Pre Stimulus Stance Durations (ms)');
    ylabel('Post Stimulus Stance Durations (ms)');
    title({trigger_type});
    ConfAxis;
    axis equal;
    xlim([0 1000]);
    ylim([0 1000]);
end
end