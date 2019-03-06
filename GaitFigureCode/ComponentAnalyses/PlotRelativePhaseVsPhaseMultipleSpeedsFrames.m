function PlotRelativePhaseVsPhaseMultipleSpeedsFrames( newData, speed_edges, corder, smoothing )
% This function calculates the relative phase as a function of limb phase for several speed conditions 

%% Define the speed edges

% Create the plot legend
legend_labels = arrayfun(@(x,y) strcat(num2str(x),'-',num2str(y),' mm/s'), ...
    speed_edges(1:end-1), speed_edges(2:end), 'uniformoutput',false);

% Define the phase edges
phase_edges = 0:.1:1;
phase_centers = (phase_edges(1:end-1)+phase_edges(2:end))/2;
phases = repmat(phase_centers',1,length(speed_edges)-1);

% Define the variable list
if smoothing
    varList = {'smooth_forwardSpeed_mmPerSec', 'smooth_InstantaneousPhase_L2y', 'smooth_InstantaneousPhase_R2y'};
    vf = newData.smooth_forwardSpeed_mmPerSec;
else
    varList = {'forwardSpeed_mmPerSec', 'InstantaneousPhase_L2y', 'InstantaneousPhase_R2y'};
    vf = newData.forwardSpeed_mmPerSec;
end

%% Loop through and calculate the timeseries variables

% Preallocate some variables
num = zeros(size(phases,1),size(phases,2));
mean_rel = zeros(size(phases,1),size(phases,2));
sem_rel = zeros(size(phases,1),size(phases,2));
std_rel = zeros(size(phases,1),size(phases,2));

% Loop through the timeseries data
for n = 1:length(speed_edges)-1
   
    
    cur_data = newData{vf >= speed_edges(n) & vf < speed_edges(n+1), varList};

    % Compute the relative phase (L2-R2)
    cur_rel = mod(cur_data(:,2) - cur_data(:,3), (2*pi))./(2*pi);

    % Convert the L2 phase to cycles
    cur_L2 = mod(cur_data(:,2),(2*pi)) / (2*pi);
    
    % Discretize the slow data removing any NaN values
    grp(:) = discretize(cur_L2(:), phase_edges);
    nan_idx = isnan(grp);
    cur_rel= cur_rel(:);
    cur_rel(nan_idx) = [];
    grp(nan_idx) = [];
    grp = grp';
    
    % Calculate the mean and SEM for the data
    cur_num = splitapply(@(x)sum(~isnan(x)),cur_rel,grp);
    mean_cur_rel = splitapply(@(x)nanmean(x),cur_rel,grp);
    std_cur_rel = splitapply(@(x)nanstd(x),cur_rel,grp);
    sem_cur_rel = splitapply(@(x)nanstd(x)./((sum(~isnan(x)))^(1/2)),cur_rel,grp);

    % Deposit the data into a single variable
    num(:,n) = cur_num;
    mean_rel(:,n) = mean_cur_rel;
    sem_rel(:,n) = sem_cur_rel;
    std_rel(:,n) = std_cur_rel;
    
    % Clear some intermediate variables
    clear grp cur_data cur_rel cur_L2
    
end

%% Create the figure
makeFigure;
PlotXvsY(phases, mean_rel, 'error', std_rel, 'color', corder); % STD
xlabel('\phi_{L2}');
ylabel('\phi_{L2}-\phi_{R2}');
title('STD Error');
ConfAxis;
PlotConstLine(.5,1);
legend(legend_labels);
xlim([0 1]);
ylim([0 1]);
axis square;

makeFigure;
PlotXvsY(phases, mean_rel, 'error', sem_rel, 'color', corder); % SEM
xlabel('\phi_{L2}');
ylabel('\phi_{L2}-\phi_{R2}');
title('SEM Error');
ConfAxis;
PlotConstLine(.5,1);
legend(legend_labels);
xlim([0 1]);
ylim([.4 .6]);
axis square;

end