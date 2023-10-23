function PlotStanceDurationDensity( newData, logscale, clims, cmap, plot_type )
%This function plots the stance duration as a function of forward speed as
%a joint density

% NOTE: This function assumes that
% computeSwingStanceAmplitudesAndDurations.m has already been run on the
% dataset.

%% Validate the inputs

if ~exist('logscale','var') || isempty(logscale)
    logscale = true;
end

if ~exist('clims','var') || isempty(clims)
    clims = [];
end

if ~exist('cmap','var') || isempty(cmap)
    cmap = viridis(256);
end

if ~exist('plot_type','var') || isempty(plot_type)
    plot_type = 'joint';
end

% Set the number of contour lines
numLvl = 3;

%% Get the data in order

% Define a limb list
limbList = {'L1','L2','L3','R1','R2','R3'};

% Define some variables for gathering the data
velVarList = {'_StanceAvgVf_mmPerSec'};

% Define some logical variables for isolating events
swingVarList = strcat(limbList, '_SwingStart');
stanceVarList = strcat(limbList, '_StanceStart');

swingStarts = newData{:,swingVarList};
stanceStarts = newData{:,stanceVarList};

% Define some velocity bin edges
vfedges = [5:1:35];

% Define some stance bin edges
stanceEdges = [-500/150:(1000/150):250];

% Define the centers for all the bins
vfCenters = (vfedges(1:end-1) + vfedges(2:end))./2;
stanceCenters = (stanceEdges(1:end-1) + stanceEdges(2:end))./2;

% Get all the stance duration data in one array
[ steps_all ] = gatherStanceDurationVelocityData(newData, limbList, velVarList, '_StanceDur_mmPerSec', stanceStarts);

%% Calculate the joint distribution

% Calculate the data
if strcmp(plot_type, 'joint')
    N = histcounts2(steps_all(:,1), steps_all(:,2), vfedges, stanceEdges, 'normalization','pdf');
    
elseif strcmp(plot_type, 'conditional')
    
    N = zeros(length(vfedges)-1, length(stanceEdges)-1);
    
    for n = 1:length(vfedges)-1
        keep = steps_all(:,1) >= vfedges(n) & steps_all(:,1) < vfedges(n+1);
        N(n,:) = histcounts(steps_all(keep,2), stanceEdges, 'normalization','pdf');
    end
else
    error('Invalid Plot Type.');
end
    
% Plot the data (Joint Distribution)
MakeFigure;

if logscale
    hold on;
    imagesc(vfCenters, stanceCenters, log10(N)');
    [xList, yList] = meshgrid(vfCenters, stanceCenters);
    contour(xList, yList, log10(N)',numLvl,'edgecolor', 'k');

else
    hold on;
    imagesc(vfCenters, stanceCenters, N');
    [xList, yList] = meshgrid(vfCenters, stanceCenters);
    contour(xList, yList, N',numLvl,'edgecolor', 'k');

end
axis('xy','square','tight');

colormap(cmap);
cbar = colorbar('southoutside');

if ~isempty(clims)
    caxis(clims);
end

if logscale
    a = cbar.Ticks;
    cbar.TickLabels = cellstr(num2str(a', '10^{%0.1f}'));
end

xlabel('Forward Speed (mm/s)');
ylabel('Stance Duration (ms)');
ylim([0 250]);
if strcmp(plot_type,'joint')
    title('P(v_f,t_{stance})');
else
    title('P(t_{stance} | v_f)');
end
ConfAxis;

%% Identify trajectories in each of the desired stance duration ranges

% Plot the data (again)
MakeFigure;

if logscale
    hold on;
    imagesc(vfCenters, stanceCenters, log10(N)');
    [xList, yList] = meshgrid(vfCenters, stanceCenters);
    contour(xList, yList, log10(N)',numLvl,'edgecolor', 'k');
else
    hold on;
    imagesc(vfCenters, stanceCenters, N');
    [xList, yList] = meshgrid(vfCenters, stanceCenters);
    contour(xList, yList, N',numLvl,'edgecolor', 'k');
end
axis('xy','square','tight');

colormap(cmap);
cbar = colorbar('southoutside');

if ~isempty(clims)
    caxis(clims);
end

if logscale
    a = cbar.Ticks;
    cbar.TickLabels = cellstr(num2str(a', '10^{%0.1f}'));
end

xlabel('Forward Speed (mm/s)');
ylabel('Stance Duration (ms)');
ylim([0 250]);
if strcmp(plot_type,'joint')
    title('P(v_f,t_{stance})');
else
    title('P(t_{stance} | v_f)');
end
ConfAxis;
hold on;

%% Annotate the dataset with the tripod trajectory

% Load the dataset
load('/Users/bdeangelis/Desktop/Datasets/_CURRENT_DATASETS/IsoD1-Masked-1000PCs/ExampleFlies_Current/Tripod-gait.mat', 'singleFly');

% Gather the step data from the tripod trajectory
stanceStarts = singleFly{:,stanceVarList};
[ tripod_steps ] = gatherStanceDurationVelocityData(singleFly, limbList, velVarList, '_StanceDur_mmPerSec', stanceStarts);

plot(tripod_steps(:,1), tripod_steps(:,2), 'linestyle','none', 'marker', '.', 'markersize', 12, 'color', 'r');


%% Annotate the plot with the tetrapod trajectory

% Load the dataset
load('/Users/bdeangelis/Desktop/Datasets/_CURRENT_DATASETS/IsoD1-Masked-1000PCs/ExampleFlies_Current/Tetrapod-gait.mat', 'singleFly');

% Gather the step data from the tripod trajectory
stanceStarts = singleFly{:,stanceVarList};
[ tetrapod_steps ] = gatherStanceDurationVelocityData(singleFly, limbList, velVarList, '_StanceDur_mmPerSec', stanceStarts);

plot(tetrapod_steps(:,1), tetrapod_steps(:,2), 'linestyle','none', 'marker', '.', 'markersize', 12, 'color', 'g');


%% Annotate the plot with the non-canonical trajectory

% Load the dataset
load('/Users/bdeangelis/Desktop/Datasets/_CURRENT_DATASETS/IsoD1-Masked-1000PCs/ExampleFlies_Current/Non-canonical-gait.mat', 'singleFly');

% Gather the step data from the tripod trajectory
stanceStarts = singleFly{:,stanceVarList};
[ nongait_steps ] = gatherStanceDurationVelocityData(singleFly, limbList, velVarList, '_StanceDur_mmPerSec', stanceStarts);

plot(nongait_steps(:,1), nongait_steps(:,2), 'linestyle','none', 'marker', '.', 'markersize', 12, 'color', 'b');


end

%% Local function for gathering all the data in one array

function [ steps_all ] = gatherStanceDurationVelocityData(newData, limbList, velVarList, var_name, logicals)

steps_all = [];

for i =1:6
    
    % Get the velocity data and associated variable
    steps = newData{logicals(:,i), [strcat(limbList(i), velVarList),strcat(limbList(i),var_name)]};
    
    % Filter out NaN values
    steps(isnan(steps(:,2)),:) = [];
    steps(isnan(steps(:,1)),:) = [];
    
    % Append to a single dataset
    steps_all = [steps_all; steps];
     
end

end

