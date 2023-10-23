function [p_ctrl, h_ctrl, stats_ctrl, p_exp, h_exp, stats_exp, h, p, ks2stat] = PlotVisualStimulusGaitFigures( newData, lineThickness, suppressPlots, num_bootstraps )
% This function generates all the plots associated with the Visual Stimulus data
% for the gait paper

%% Generate the relevant data

% Generate the data associated with the control condition
[ p_ctrl, h_ctrl, stats_ctrl, data_ctrl, time_ctrl, meanSpeed_ctrl, errorSpeed_ci_ctrl, normMeanSpeed_ctrl, errorNormSpeed_ci_ctrl, ...
    postDurations_pdf_ctrl, post_pdf_ci_ctrl, postDurations_cdf_ctrl, post_cdf_ci_ctrl, post_stances_ctrl ] ...
    = PlotOptoStanceIncreaseAllLimbs( newData, 'Control', 'Visual', 1, suppressPlots, num_bootstraps );

% Generate the data associated with the experimental condition
[ p_exp, h_exp, stats_exp, data_exp, time_exp, meanSpeed_exp, errorSpeed_ci_exp, normMeanSpeed_exp, errorNormSpeed_ci_exp, ...
    postDurations_pdf_exp, post_pdf_ci_exp, postDurations_cdf_exp, post_cdf_ci_exp, post_stances_exp ] ...
    = PlotOptoStanceIncreaseAllLimbs( newData, 'Visual', 'Visual', 1, suppressPlots, num_bootstraps );

%% Truncate the time-series data prior to plotting (Remove the extra pre-stimulus frames)

% Define the starting frame for the plots
start_idx = 36;

% Truncate the control time-series
time_ctrl = time_ctrl(start_idx:end);
meanSpeed_ctrl = meanSpeed_ctrl(start_idx:end);
errorSpeed_ci_ctrl = errorSpeed_ci_ctrl(:,start_idx:end);
normMeanSpeed_ctrl = normMeanSpeed_ctrl(start_idx:end); 
errorNormSpeed_ci_ctrl = errorNormSpeed_ci_ctrl(:,start_idx:end);

% Truncate the experimental time-series
time_exp = time_exp(start_idx:end);
meanSpeed_exp = meanSpeed_exp(start_idx:end);
errorSpeed_ci_exp = errorSpeed_ci_exp(:,start_idx:end);
normMeanSpeed_exp = normMeanSpeed_exp(start_idx:end); 
errorNormSpeed_ci_exp = errorNormSpeed_ci_exp(:,start_idx:end);

%% Plot the control and experimental mean walking speed time-series

% Define the legend labels
legend_labels{1} = strcat('Random Trigger n=', num2str(size(data_ctrl, 1)/6));
legend_labels{2} = strcat('Experiment n=', num2str(size(data_exp, 1)/6));

% Unnormalized Walking Speed
MakeFigure;
PlotConfidenceIntervalWithErrorPatch( [time_ctrl, time_exp], [meanSpeed_ctrl, meanSpeed_exp], [errorSpeed_ci_ctrl(1,:);errorSpeed_ci_exp(1,:)]', [errorSpeed_ci_ctrl(2,:);errorSpeed_ci_exp(2,:)]');
title({'Visual Stimulus', 'Unnormalized Walking Speed'});
xlabel('Time (ms)');
ylabel('Forward Speed (mm/s)');
ylim([-5 25]);
xlim([min(time_ctrl), max(time_ctrl)]);
PlotConstLine(0,2);
legend(legend_labels);
ConfAxis;

%% Plot the control and experimental normalized walking speed time-series

% Define the legend labels
legend_labels{1} = strcat('Random Trigger n=', num2str(size(data_ctrl, 1)/6));
legend_labels{2} = strcat('Experiment n=', num2str(size(data_exp, 1)/6));

% Normalized Walking Speed
MakeFigure;
PlotConfidenceIntervalWithErrorPatch( [time_ctrl, time_exp], [normMeanSpeed_ctrl, normMeanSpeed_exp], [errorNormSpeed_ci_ctrl(1,:);errorNormSpeed_ci_exp(1,:)]', [errorNormSpeed_ci_ctrl(2,:);errorNormSpeed_ci_exp(2,:)]');
title({'Visual Stimulus', 'Normalized Walking Speed'});
xlabel('Time (ms)');
ylabel('Forward Speed (mm/s)');
ylim([0 2]);
xlim([min(time_ctrl), max(time_ctrl)]);
PlotConstLine(0,2);
legend(legend_labels);
ConfAxis;

%% Plot the control and experimental post stimulus stance distributions (pdf)

% Define the colormap
cmap = lines(2);

% Define the bin edges
dur_binEdges = [-500/150:1000/150:250]; % 1 frame bins
dur_binCenters = (dur_binEdges(1:end-1) + dur_binEdges(2:end))/2;

% Define the legend labels
legend_labels{1} = strcat('Random Trigger n=', num2str(size(data_ctrl, 1)));
legend_labels{2} = strcat('Experiment n=', num2str(size(data_exp, 1)));

% Plot
MakeFigure;
hold on;
PlotConfidenceIntervalWithErrorPatch( [dur_binCenters;dur_binCenters]', [postDurations_pdf_ctrl;postDurations_pdf_exp]', ...
    [post_pdf_ci_ctrl(1,:);post_pdf_ci_exp(1,:)]', [post_pdf_ci_ctrl(2,:);post_pdf_ci_exp(2,:)]', cmap);
title({'Visual','Stance Duration Distributions - All Limbs'});
xlabel('Stance Duration');
ylabel('pdf');
% ylim([0 .02]);
legend(legend_labels); 
ConfAxis;

%% Plot the control and experimental post stimulus stance distributions (cdf)

% Perform a KS test on the stance durations
[h, p, ks2stat] = kstest2(post_stances_ctrl(:), post_stances_exp(:));

% Define the legend labels
legend_labels{1} = strcat('Random Trigger n=', num2str(size(data_ctrl, 1)));
legend_labels{2} = strcat('Experiment n=', num2str(size(data_exp, 1)));

% Plot
MakeFigure;
hold on;
PlotConfidenceIntervalWithErrorPatch( [dur_binCenters;dur_binCenters]', [postDurations_cdf_ctrl;postDurations_cdf_exp]', ...
    [post_cdf_ci_ctrl(1,:);post_cdf_ci_exp(1,:)]', [post_cdf_ci_ctrl(2,:);post_cdf_ci_exp(2,:)]', cmap);
plot([mean(data_ctrl(:,2)); mean(data_ctrl(:,2))],[0 1],'--','Color', cmap(1,:), 'linewidth', lineThickness);
plot([mean(data_exp(:,2));  mean(data_exp(:,2))], [0 1],'--','Color', cmap(2,:), 'linewidth', lineThickness);
title({'Visual','Stance Duration Distributions - All Limbs',sprintf('p = %.4e',p),sprintf('D = %.4f',ks2stat)});
xlabel('Stance Duration');
ylabel('cdf');
% ylim([0 1]);
legend(legend_labels); 
ConfAxis;

end