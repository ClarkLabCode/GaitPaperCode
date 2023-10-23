function PlotStanceConfigurationFractions( newData, analysisType, binEdges, corder, faceAlpha, minorLineWidth, smoothing, removeDiscontinuities, limbVarType, suppressExtraPlots )
% PlotStanceConfigurationfractions: A function to plot the fraction of each
% configuration of limbs down in bins defined by centroid dynamics.

%% Validate inputs

if ~exist('analysisType','var') || isempty(analysisType)
    analysisType = 'linear';
end

if ~exist('binEdges','var') || isempty(binEdges)
    switch analysisType
        case 'linear'
            binEdges = 0:0.5:45;
        case 'forward'
            binEdges = 0:0.5:35;
        case 'angular'
            binEdges = -1000:10:1000;
        case 'absangular'
            binEdges = 0:10:1000;
        otherwise
            error('Invalid analysisType: %s', analysisType);
    end
end

if ~exist('corder','var') || isempty(corder)
%     corder = linspecer(7,'sequential');
    corder = jet(7) .* 0.8;
end

if ~exist('faceAlpha', 'var') || isempty(faceAlpha)
    faceAlpha = 0.6;
end

if ~exist('minorLineWidth', 'var') || isempty(minorLineWidth)
    minorLineWidth = 2;
end

if ~exist('smoothing','var') || isempty(smoothing)
    smoothing = false;
end

if ~exist('removeDiscontinuities','var') || isempty(removeDiscontinuities)
    removeDiscontinuities = false;
end

if ~exist('suppressExtraPlots','var') || isempty(suppressExtraPlots)
    suppressExtraPlots = false;
end

%% Bin the data by the desired dynamical variable

% Get the needed data
[ v ] = getVelocityData( newData, analysisType, smoothing );

% Get the centers of the bins
binCenters = ( binEdges(1:end-1) + diff(binEdges)/2 )';

% Compute the distribution of velocities and discretize the velocity
% appropriately
[N, ~, vbin] = histcounts( v, binEdges, 'normalization','probability');
vbin(vbin==0) = NaN;

%% Compute fractions of each stance configuration in each velocity bin

% Set the list of variables to use

if strcmp(limbVarType, 'Egocentric')
    % Egocentric Frame Logical
    varList = {'L1_down', 'L2_down', 'L3_down','R1_down', 'R2_down', 'R3_down'};
elseif strcmp(limbVarType, 'Camera')
    % Camera Frame Logical
    varList = {'L1_down_cam', 'L2_down_cam', 'L3_down_cam','R1_down_cam', 'R2_down_cam', 'R3_down_cam'};
else
    error('Invalid limb variable type.');
end

% Extract the data
X = newData{:, varList};

% Remove single-frame discontinuities in the up/down logicals
if removeDiscontinuities
    X = splitapply(@(x) {removeSingleFrameUpDownDiscontinuities( x )}, X, findgroups(newData.uniqueFlyTrajID));
    X = cell2mat(X);
end

% Compute the total number of feet down
totalFeetDown = sum(X, 2);

% Compute the stance index (re-defined to match newer standards for left limbs before right)
stanceIndex = sum(X .* (2 .^ (0:5)), 2);

% Make the key linking stanceIndex to totalFeetDown
feetDownKey = sum(de2bi([0:63]),2);

% Get the unique groups defined by the discretized velocities
G = findgroups(vbin);

% Compute major edges (total feet down)
totalFeetDownEdges = (0:7) - 0.5;
majorPercentage = splitapply(@(x) {histcounts(x, totalFeetDownEdges, 'normalization','probability')}, totalFeetDown, G);
majorPercentage = cell2mat(majorPercentage);

% Take the cumulative sum
cumulativeMajorPercentage = cumsum(majorPercentage, 2);

% Append zeros to major percentages
cumulativeMajorPercentage = [zeros(size(cumulativeMajorPercentage, 1), 1), cumulativeMajorPercentage];

% Compute minor edges (specific combinations of feet down)
stanceEdges = (0:64) - 0.5;
minorPercentage = splitapply(@(x) {histcounts(x, stanceEdges, 'normalization','probability')}, stanceIndex, G);
minorPercentage = cell2mat(minorPercentage);

% Sort by the number of feet down
[feetDownKey, idx] = sort(feetDownKey, 'ascend');
stanceEdgeKey = (0:63);
stanceEdgeKey = stanceEdgeKey(idx);
minorPercentage = minorPercentage(:, idx);

% Take the cumulative sum
cumulativeMinorPercentage = cumsum(minorPercentage, 2);

%% Plot with major percentages only

legendStr = {'No feet down','1 foot down','2 feet down', '3 feet down', '4 feet down', '5 feet down','6 feet down'};

figure('Position',[200,500,1000,1000],'WindowStyle','docked');
hold on;

% Plot nonsense markers to make legend
for ind=1:7
    plot(NaN, NaN, 'Linewidth', 4, 'Color', corder(ind,:));
end

% Plot major percentages as patches
x = [binCenters; flipud(binCenters)];
for ind = 2:8
    y = [cumulativeMajorPercentage(:, ind-1); flipud( cumulativeMajorPercentage(:, ind) )];
    patch(x, y, 1, 'FaceColor', corder(ind-1,:), 'FaceAlpha', 1, 'EdgeColor', 'None');
end

configureStanceFractionPlotAxes( binCenters, analysisType );

legend( legendStr, 'Location','Northwest' );

%% Plot with minor percentages

figure('Position',[200,500,1000,1000],'WindowStyle','docked');
hold on;

% Plot nonsense markers to make legend
for ind=1:7
    plot(NaN, NaN, 'Linewidth', 4, 'Color', corder(ind,:));
end

% Plot major percentages as patches
x = [binCenters; flipud(binCenters)];
for ind = 2:8
    y = [cumulativeMajorPercentage(:, ind-1); flipud( cumulativeMajorPercentage(:, ind) )];
    patch(x, y, 1, 'FaceColor', corder(ind-1,:), 'FaceAlpha', faceAlpha, 'EdgeColor', 'None');
end

% Plot minor percentages as lines
for ind = 1:7
    plot(binCenters, cumulativeMinorPercentage(:, feetDownKey == (ind - 1)), 'Linewidth', minorLineWidth, 'Color', corder(ind,:) );
end

configureStanceFractionPlotAxes( binCenters, analysisType );

legend( legendStr, 'Location','Northwest' );

%% Plot with major percentages only, with a histogram of the dynamical variable on which binning is performed

if ~suppressExtraPlots
    
    figure('Position',[200,500,1000,1000],'WindowStyle','docked');
    
    s(1) = subplot(2,1,1);
    plot(binCenters, N, 'linewidth', 2);
    s(1).Position(4) = 0.2;
    xticks([]);
    xlim([min(binCenters), max(binCenters)]);
    ylabel('probability');
    ConfAxis;
    
    s(2) = subplot(2,1,2); hold on;
    
    % Plot nonsense markers to make legend
    for ind=1:7
        plot(NaN, NaN, 'Linewidth', 4, 'Color', [corder(ind,:), 1]);
    end
    
    % Plot major percentages as patches
    x = [binCenters; flipud(binCenters)];
    for ind = 2:8
        y = [cumulativeMajorPercentage(:, ind-1); flipud( cumulativeMajorPercentage(:, ind) )];
        patch(x, y, 1, 'FaceColor', corder(ind-1,:), 'FaceAlpha', 1, 'EdgeColor', 'None');
    end
    
    configureStanceFractionPlotAxes( binCenters, analysisType );
    
    s(2).Position(4) = s(2).Position(4)+0.2;
    s(1).Position(2) = s(2).Position(2)+s(2).Position(4)+0.025;
    s(1).Position(3) = s(2).Position(3);
    
    legend( legendStr, 'Location','Northwest' );
    
end

%% Plot with minor percentages, with a histogram of the dynamical variable on which binning is performed

if ~suppressExtraPlots
    
    figure('Position',[200,500,1000,1000],'WindowStyle','docked');
    
    s(1) = subplot(2,1,1);
    plot(binCenters, N, 'linewidth', 2);
    s(1).Position(4) = 0.2;
    xticks([]);
    xlim([min(binCenters), max(binCenters)]);
    ylabel('probability');
    ConfAxis;
    
    s(2) = subplot(2,1,2); hold on;
    
    % Plot nonsense markers to make legend
    for ind=1:7
        plot(NaN, NaN, 'Linewidth', 4, 'Color', [corder(ind,:), 1]);
    end
    
    % Plot major percentages as patches
    x = [binCenters; flipud(binCenters)];
    for ind = 2:8
        y = [cumulativeMajorPercentage(:, ind-1); flipud( cumulativeMajorPercentage(:, ind) )];
        patch(x, y, 1, 'FaceColor', corder(ind-1,:), 'FaceAlpha', faceAlpha, 'EdgeColor', 'None');
    end
    
    % Plot minor percentages as lines
    for ind = 1:7
        plot(binCenters, cumulativeMinorPercentage(:, feetDownKey == (ind - 1)), 'Linewidth', minorLineWidth, 'Color', corder(ind,:) );
    end
    
    configureStanceFractionPlotAxes( binCenters, analysisType );
    
    s(2).Position(4) = s(2).Position(4)+0.2;
    s(1).Position(2) = s(2).Position(2)+s(2).Position(4)+0.025;
    s(1).Position(3) = s(2).Position(3);
    
    legend( legendStr, 'Location','Northwest' );
    
end

%% Plot non-cumulative percentages of each total number of feet down

if ~suppressExtraPlots
    
    figure('Position',[200,500,1000,1000],'WindowStyle','docked');
    hold on;
    set(gca, 'ColorOrder', corder);
    plot(binCenters, majorPercentage, 'linewidth', 2);
    
    configureStanceFractionPlotAxes( binCenters, analysisType );
    
    ylim([0, 1]);
    legend( legendStr );
    
end

%% Plot non-cumulative percentages of individual configurations, for each number of feet down

if ~suppressExtraPlots
    
    seriesLabels = cellstr( dec2bin( stanceEdgeKey ) );
    
    for ind = 1:7
        figure('Position',[200,500,1000,1000],'WindowStyle','docked');
        
        hold on;
        if nnz(feetDownKey == (ind-1)) < 12
            set(gca, 'ColorOrder', linspecer( nnz(feetDownKey == (ind-1)), 'qualitative'));
            
        else
            set(gca, 'ColorOrder', vega20(nnz(feetDownKey == (ind-1))));
        end
        plot(binCenters, minorPercentage(:, feetDownKey == (ind-1) ), 'linewidth', 2 );
        title(legendStr{ind});
        
        configureStanceFractionPlotAxes( binCenters, analysisType );
        
        try
            ylim([0 max( max( minorPercentage(:, feetDownKey == (ind-1) ), [], 2), [], 1)]);
        catch
        end
        
        legend(seriesLabels(feetDownKey == (ind-1)), 'Location','Northwest');
    end
    
end

%% Plot non-cumulative percentages of individual configureations (3,4,5) Colored by number
% NOTE: The coloring here matches the stacked plot so that we can use the same legend
% NOTE: All stances are loaded into the minorPercentage by number of feet
% down. To get a particular stance we need to search based on the binary
% value.

% Get the forward speeds associated with the max percentage of N feet down
max3 = binCenters(find(majorPercentage(:,4) == max(majorPercentage(:,4))));
max4 = binCenters(find(majorPercentage(:,5) == max(majorPercentage(:,5))));
max5 = binCenters(find(majorPercentage(:,6) == max(majorPercentage(:,6))));

MakeFigure;

% Plot the canonical 3 feet down combinations
% [010101, 101010] --> Tripod states (21, 42)
threeList = [21, 42];
for n = 1:length(threeList)
    hold on;
    plot(binCenters, minorPercentage(:,stanceEdgeKey==threeList(n)), 'linewidth', 1, 'color', corder(4,:));
end
plot([max3,max3],[0,.4],'color', corder(4,:), 'linewidth', 1, 'linestyle','--');

% Plot the canonical 4 feet down combinations
% [110101, 101110] --> Mid-hind combo (53, 46)
% [011101, 101011] --> Fore-mid combo (29, 43)
% [011110, 110011] --> Fore-hind combo (30, 51)
fourList = [53,46,29,43,30,51];
for n = 1:length(fourList)
    hold on;
    plot(binCenters, minorPercentage(:,stanceEdgeKey==fourList(n)), 'linewidth', 1, 'color', corder(5,:));
end
plot([max4,max4],[0,.4],'color', corder(5,:), 'linewidth', 1, 'linestyle','--');

% Plot the canoncial 5 feet down combinations
% [011111, 111011] --> Forelimb Swing (31,59)
% [101111, 111101] --> Midlimb Swing (47,61)
% [110111, 111110] --> Hindlimb Swing (55,62)
fiveList = [31,59,47,61,55,62];
for n = 1:length(fiveList)
    hold on;
    plot(binCenters, minorPercentage(:,stanceEdgeKey==fiveList(n)), 'linewidth', 1, 'color', corder(6,:));
end
plot([max5,max5],[0,.4],'color', corder(6,:), 'linewidth', 1, 'linestyle','--');
    
configureStanceFractionPlotAxes( binCenters, analysisType );
ylim([0, .4]);

end

%% Local Functions

function configureStanceFractionPlotAxes( binCenters, analysisType )

switch analysisType
    case 'linear'
        xlabel('speed (mm/s)');
        ylabel('P( n_{down} | speed )');
    case 'forward'
        xlabel('v_{||} (mm/s)');
        ylabel('P( n_{down} | v_{||} )');
    case 'angular'
        xlabel('v_r (\circ/s)');
        ylabel('P( n_{down} | v_{r} )');
    case 'absangular'
        xlabel('|v_r| (\circ/s)');
        ylabel('P( n_{down} | |v_{r}| )');
end


xlim([min(binCenters), max(binCenters)]);
ylim([0 1]);

ConfAxis;
end

function [ x ] = removeSingleFrameUpDownDiscontinuities( x )
% Remove single - frame discontinuities
if size(x,1) > 3
    for ind=2:size(x,1)-1
        for ii=1:size(x,2)
            if (x(ind,ii) ~= x(ind-1,ii)) && (x(ind-1,ii) == x(ind+1,ii))
                x(ind,ii) = x(ind-1,ii);
            end
        end
    end
end
end

function [ v ] = getVelocityData( newData, analysisType, smoothing )
% Extract the appropriate dynamical data given the analysis type

switch analysisType
    case 'linear'
        velVarList = {'linVel_mmPerSec'};
        
        if smoothing
            velVarList = strcat('smooth_', velVarList);
        end
        
        v = newData{:, velVarList};
        
    case 'forward'
        velVarList = {'forwardSpeed_mmPerSec'};
        
        if smoothing
            velVarList = strcat('smooth_', velVarList);
        end
        
        v = newData{:, velVarList};
        
    case 'angular'
        velVarList = {'angVel_radPerSec'};
        
        if smoothing
            velVarList = strcat('smooth_', velVarList);
        end
        
        v = rad2deg( newData{:, velVarList} );
    case 'absangular'
        velVarList = {'angVel_radPerSec'};
        
        if smoothing
            velVarList = strcat('smooth_', velVarList);
        end
        
        v = abs( rad2deg( newData{:, velVarList} ) );
    otherwise
        error('Invalid analysisType: %s', analysisType);
end

end

