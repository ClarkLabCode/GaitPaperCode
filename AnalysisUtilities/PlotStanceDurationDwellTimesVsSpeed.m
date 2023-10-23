function PlotStanceDurationDwellTimesVsSpeed( newData, vfBinEdges, corder, limbVarType )
% The probability of each stance duration's dwell time as a function of
% forward speed.

%% NOTE: Define the desired stance configurations for analysis 

% Plot the canoncial 5 feet down combinations
% [011111, 111011] --> Forelimb Swing (31,59)
% [101111, 111101] --> Midlimb Swing (47,61)
% [110111, 111110] --> Hindlimb Swing (55,62)

% Tetrapod Mid-Hind Pairings
% [1,1,0,1,0,1] --> 53
% [1,0,1,1,1,0] --> 46
% [110101, 101110] --> Mid-hind combo (53, 46)

% Tetrapod Fore-Mid Pairings
% [0,1,1,1,0,1] --> 29
% [1,0,1,0,1,1] --> 43
% [011101, 101011] --> Fore-mid combo (29, 43)

% Tetrapod Fore-Hind Pairings
% [0,1,1,1,1,0] --> 30
% [1,1,0,0,1,1] --> 51
% [011110, 110011] --> Fore-hind combo (30, 51)

% Tripod Pairings
% [0,1,0,1,0,1] --> 21
% [1,0,1,0,1,0] --> 42
% [010101, 101010] --> Tripod states (21, 42)

%% Validate inputs

if ~exist('vfBinEdges','var') || isempty(vfBinEdges)
    vfBinEdges = [15, 25];
end

if ~exist('corder','var') || isempty(corder)
%     corder = linspecer(7,'sequential');
    corder = jet(7) .* 0.8;
    
end

%% Set some parameters for the analyses and plots

% Put all the desired states into a single list
% stateList = [53, 46, 29, 43, 30, 51, 21, 42]; % No pentapod
% stateList = [31, 59, 47, 61, 55, 62, 53, 46, 29, 43, 30, 51, 21, 42]; % 5,4,3 Feet down combinations
numFeet = [3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5]; % 3,4,5 Feet down combinations
stateList = [21, 42, 53, 46, 29, 43, 30, 51, 31, 59, 47, 61, 55, 62]; % 3,4,5 Feet down combinations

% Set some additional parameters for the plots
durBinEdges = [.5:1:7.5] .* 100/15; % Individual frames in msec

%% Calculate the necessary variables

% Convert the set of footfall logicals into a number between 0-63 corresponding to each gait
flyID = newData.uniqueFlyTrajID;

if strcmp(limbVarType, 'Egocentric')
    % Egocentric Frame Binaries
    states = uint16(newData.L1_down * (2^5) + newData.L2_down * (2^4) + newData.L3_down * (2^3) + ...
        newData.R1_down * (2^2) + newData.R2_down * (2^1) + newData.R3_down * (2^0));
elseif strcmp(limbVarType, 'Camera')
    % Camera Frame Binaries
    states = uint16(newData.L1_down_cam * (2^5) + newData.L2_down_cam * (2^4) + newData.L3_down_cam * (2^3) + ...
        newData.R1_down_cam * (2^2) + newData.R2_down_cam * (2^1) + newData.R3_down_cam * (2^0));
else
    error('Invalid limb variable type.');
end

prevStates = states(1:end-2);
nextStates = states(3:end);
states = states(2:end-1);

% Get the centroid velocities associated with each of these stances
vel_parallel = newData.forwardSpeed_mmPerSec;
vel_parallel = vel_parallel(2:end-1);

% Define a grouping variable for each trajectory in the dataset
G = findgroups(newData.uniqueFlyTrajID);
G = G(2:end-1);

% NOTE: Kill all the transitions between unique trajectories
trajChange = flyID(1:end-1) ~= flyID(2:end);
trajChange = trajChange(1:end-1);
flyID = flyID(2:end-1);

for n = 1:length(stateList)
        
    % Get all the locations
    startLocs = splitapply(@(x,y) {find(x==stateList(n) & x ~= y)}, states, prevStates, G);
    endLocs = splitapply(@(x,y) {find(x == stateList(n) & x ~= y)}, states, nextStates, G);
    vel_par_cell = splitapply(@(x) {x}, vel_parallel, G); 
    
    % Get the miniumum length for each of these events
    lengths = cellfun(@(x,y){min([length(x),length(y)])}, startLocs, endLocs);
    
    % Get the durations associated with each of these state events
    durs = cellfun(@(len,srt,stp) get_durations_local(len,srt,stp), lengths, startLocs, endLocs, 'uniformoutput', false); 
    
    % Get the mean parallel body velocity for each of these state events
    vels = cellfun(@(vel,srt,stp,len) get_body_vel_local(vel,srt,stp,len), vel_par_cell, startLocs, endLocs, lengths, 'uniformoutput', false); 
    
    % Concatenate all the cells in durs into a single list
    idx = cellfun('isempty',durs);
    durs = durs(~idx);
    durs = cat(1, durs{:});
    
    % Convert the durations from frames into milliseconds
    durs = durs .* (100/15); 
    duration{n} = durs;
    
    % Concatenate all the cells in vels into a single list
    idx = cellfun('isempty',vels);
    vels = vels(~idx);
    vels = cat(1, vels{:});
    velocity{n} = vels;
    
%     % TEST: Plot a histogram of the distribution of step durations in frames
%     MakeFigure;
%     histogram(durs);
%     title(state_labels{n});
%     xlabel('Duration (Frames)');
%     ylabel('Counts');
%     ConfAxis;
    
end

%% Calculate the input data for the plots

% Calculate the data for populating the charts (All walking speeds)
numStances = length(duration);
for n = 1:numStances
    N(:,n) = histcounts(duration{n}, durBinEdges,'Normalization', 'pdf');
    Ncdf(:,n) = histcounts(duration{n}, durBinEdges,'Normalization', 'cdf');
end

% Calculate the data for populating the charts (Conditioned on walking speeds)
numSpeeds = length(vfBinEdges) + 1;
speeds = [-inf, vfBinEdges, inf];
for m = 1:numSpeeds
    for n = 1:numStances
        curDur = duration{n};
        curVel = velocity{n};
        curDur = curDur(speeds(m) < curVel & curVel <= speeds(m+1));
        N_binned(:, n, m) = histcounts(curDur, durBinEdges,'Normalization', 'pdf');
        Ncdf_binned(:, n, m) = histcounts(curDur, durBinEdges,'Normalization', 'cdf');
    end
end

%% PDF Plots

% Define the time axis for the plots
time = round((durBinEdges(1:end-1)+durBinEdges(2:end))/2,2);

% Plot the distribution of all selected stance configurations in a single plot (All walking speeds)
MakeFigure;
for n = 1:numStances
    hold on;
    plot(time, N(:,n), 'linewidth', 2, 'color', corder(numFeet(n)+1,:));
end
configureHistogramPlotAxes('Stance Duration Dwell Times', 'all');

title_text = {strcat(num2str(0),'< v_{||} <=', num2str(vfBinEdges(1))), ...
              strcat(num2str(vfBinEdges(1)),'< v_{||} <=', num2str(vfBinEdges(2))), ...
              strcat(num2str(vfBinEdges(2)),'< v_{||}')};

% Plot the distribution of all selected stance configurations in a single plot (Conditioned on walking speeds) 
MakeFigure;
for m = 1:numSpeeds
    subplot(1,3,m);
    for n = 1:numStances
        hold on;
        plot(time, N_binned(:, n, m), 'linewidth', 2, 'color', corder(numFeet(n)+1,:));
    end    
    configureHistogramPlotAxes(title_text{m}, 'partial');
end

%% CDF Plots

% Define the time axis for the plots
time = round((durBinEdges(1:end-1)+durBinEdges(2:end))/2,2);

% Append 0,0 to all the cdf plots
time = [0, time];
Ncdf = [zeros(1,14);Ncdf];
Ncdf_binned = cat(1,zeros(1,14,3),Ncdf_binned);

% Plot the distribution of all selected stance configurations in a single plot (All walking speeds)
MakeFigure;
for n = 1:numStances
    hold on;
    plot(time, Ncdf(:,n), 'linewidth', 2, 'color', corder(numFeet(n)+1,:));
end
configureHistogramPlotAxes('Stance Duration Dwell Times', 'cdf');

title_text = {strcat(num2str(0),'< v_{||} <=', num2str(vfBinEdges(1))), ...
              strcat(num2str(vfBinEdges(1)),'< v_{||} <=', num2str(vfBinEdges(2))), ...
              strcat(num2str(vfBinEdges(2)),'< v_{||}')};

% Plot the distribution of all selected stance configurations in a single plot (Conditioned on walking speeds) 
MakeFigure;
for m = 1:numSpeeds
    subplot(1,3,m);
    for n = 1:numStances
        hold on;
        plot(time, Ncdf_binned(:, n, m), 'linewidth', 2, 'color', corder(numFeet(n)+1,:));
    end    
    configureHistogramPlotAxes(title_text{m}, 'cdf');
end

end

%% Local functions

% Local function for converting the stance starts and stops into durations
% in number of frames
function dur = get_durations_local(len, srt, stp)

try
    if srt(1) <= stp(1)
        dur = [(stp(1:len)-srt(1:len))+1];
    else
        dur = [(stp(2:len)-srt(1:len-1))+1];
    end
catch
    dur = [];
end

end

% Calculate the mean body velocity over the current step
function vels = get_body_vel_local(vel, srt, stp, len)

try
    if srt(1) <= stp(1)
        vels = arrayfun(@(a,b) sum(vel(a:b))/length(a:b), srt(1:len), stp(1:len));
    else
        vels = arrayfun(@(a,b) sum(vel(a:b))/length(a:b), srt(1:len-1), stp(2:len));
    end
catch
    vels = [];
end

end

% Local function to configure the histogram plots
function configureHistogramPlotAxes(title_text, plotType)


xlabel('Duration (ms)');
switch plotType
    case 'all'
        ylabel('P(Stance Configuration) ms^{-1}');
        xlim([0 50]);
        ylim([0 .15]);
    case 'partial'
        ylabel({'P(Stance Configuration) | v_{||} ) ms^{-1}'});
        xlim([0 50]);
        ylim([0 .15]);
    case 'cdf'
        ylabel('cdf');
        xlim([0 50]);
        ylim([0 1]);
end

title(title_text);
ConfAxis;
axis square;

end