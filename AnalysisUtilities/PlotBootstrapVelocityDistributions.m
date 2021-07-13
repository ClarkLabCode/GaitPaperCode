function PlotBootstrapVelocityDistributions( newData, num_bootstraps, ...
    speedBinEdges, vrBinEdges, vfBinEdges, vpBinEdges, corder, cmap, linewidth, smoothing )
% This function allows you to bootstrap distributions of the desired
% variables

%% Validate inputs

if ~exist('speedBinEdges','var') || isempty(speedBinEdges)
    speedBinEdges = linspace(0, 35, 70);
end
speedBinCenters = speedBinEdges(1:end-1) + diff(speedBinEdges)/2;

if ~exist('vrBinEdges','var') || isempty(vrBinEdges)
    vrBinEdges = linspace(-300, 300, 100);
end
vrBinCenters = vrBinEdges(1:end-1) + diff(vrBinEdges)/2;

if ~exist('vfBinEdges','var') || isempty(vfBinEdges)
    vfBinEdges = linspace(-5, 35, 100);
end
vfBinCenters = vfBinEdges(1:end-1) + diff(vfBinEdges)/2;

if ~exist('vpBinEdges','var') || isempty(vpBinEdges)
    vpBinEdges = linspace(-10, 10, 100);
end
vpBinCenters = vpBinEdges(1:end-1) + diff(vpBinEdges)/2;

if ~exist('corder','var') || isempty(corder)
    corder = lines(6);
end

if ~exist('cmap','var') || isempty(cmap)
    cmap = viridis(256);
end

if ~exist('linewidth','var') || isempty(linewidth)
    linewidth = 2;
end

if ~exist('smoothing','var') || isempty(smoothing)
    smoothing = true;
end

%% Check if we are using the smoothed version of the data
if smoothing && ~ismember('smooth_forwardSpeed_mmPerSec', newData.Properties.VariableNames)
    % Calculate smoothed versions of the phase and frequency variables and the centroid velocity components
    [ newData ] = filterData( newData, [], [], [], [] );
end

%% Perform the bootstraps on each of the speed variables

% Define the list of dynamical variables
varList = {'linVel_mmPerSec','angVel_radPerSec','forwardSpeed_mmPerSec','translationalSpeed_mmPerSec'};
binEdges ={speedBinEdges, vrBinEdges, vfBinEdges, vpBinEdges};

% Check whether smoothed variables are desired
if smoothing
    varList = strcat('smooth_', varList);
end

% Define the grouping variable for bootstrapping (videos)
G = findgroups(newData.videoID);

% Loop through each of the variables of interest to generate the confidence intervals
for n = 1:length(varList)

    % Get the speed data for the current analysis
    cur_data = splitapply(@(x){x(:)}, newData{:,varList{n}}, G);
    
    if strcmp(varList{n}, 'angVel_radPerSec') || strcmp(varList{n}, 'smooth_angVel_radPerSec')
        % Define the bootstrapping function used in this analysis
        bootfun = @(x) histcounts(rad2deg(cell2mat(x)), binEdges{n}, 'normalization', 'pdf');
        
        % Compute the mean associated with each dataset
        means{n} = histcounts(rad2deg(cell2mat(cur_data)), binEdges{n}, 'normalization','pdf');
    
    else    
        % Define the bootstrapping function used in this analysis
        bootfun = @(x) histcounts(cell2mat(x), binEdges{n}, 'normalization', 'pdf');
        
        % Compute the mean associated with each dataset
        means{n} = histcounts(cell2mat(cur_data), binEdges{n}, 'normalization','pdf');
    
    end

    % Compute the bootstraps
    speed_bounds{n} = bootci(num_bootstraps, bootfun, cur_data);

end

%% Plot all the data

% Linear Velocity
MakeFigure;
PlotConfidenceIntervalWithErrorPatch( speedBinCenters', means{1}', ...
    speed_bounds{1}(2,:)', speed_bounds{1}(1,:)', corder(1,:));
set(gca, 'YScale', 'log');
xlim([min(speedBinEdges), max(speedBinEdges)]);
xlabel('speed (mm/s)');
ylabel('log pdf (s/mm)');
% axis('square');
ConfAxis;

% Angular Velocity
MakeFigure;
PlotConfidenceIntervalWithErrorPatch( vrBinCenters', means{2}', ...
    speed_bounds{2}(2,:)', speed_bounds{2}(1,:)', corder(2,:));
% xlim([min(vrBinEdges), max(vrBinEdges)]);
set(gca, 'YScale', 'log');
xlabel('v_r (\circ/s)');
ylabel('log pdf (s/\circ)');
% axis('square');
ConfAxis;

% Forward Speed
MakeFigure;
PlotConfidenceIntervalWithErrorPatch( vfBinCenters', means{3}', ...
    speed_bounds{3}(2,:)', speed_bounds{3}(1,:)', corder(3,:));
% xlim([min(vfBinEdges), max(vfBinEdges)]);
set(gca, 'YScale', 'log');
xlabel('v_{||} (mm/s)');
ylabel('log pdf (s/mm)');
ConfAxis;

% Lateral Speed
MakeFigure;
PlotConfidenceIntervalWithErrorPatch( vpBinCenters', means{4}', ...
    speed_bounds{4}(2,:)', speed_bounds{4}(1,:)', corder(4,:));
% xlim([min(vpBinEdges), max(vpBinEdges)]);
set(gca, 'YScale', 'log');
xlabel('v_{\perp} (mm/s)');
ylabel('log pdf (s/mm)');
% axis('square');
ConfAxis;

end

