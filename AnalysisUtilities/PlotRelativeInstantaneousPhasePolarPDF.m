function PlotRelativeInstantaneousPhasePolarPDF( newData, nbins, corder_pairwise, corder_all, rlims, vfBinEdges, smoothing, suppressExtraPlots )
%PlotRelativeInstantaneousPhaseJointDistributions: A function to plot
%polar histograms of relative instantaneous phase.

%% Validate inputs

if ~exist('nbins','var') || isempty(nbins)
    nbins = 180;
end

if ~exist('corder_pairwise','var') || isempty(corder_pairwise)
    corder_pairwise = lines(7);
end

if ~exist('corder_all','var') || isempty(corder_all)
    corder_all = linspecer(6);
end

if ~exist('rlims','var') || isempty(rlims)
    rlims = [];
end

if ~exist('vfBinEdges','var') || isempty(vfBinEdges)
    vfBinEdges = [0, 10, 20, 30];
end

if ~exist('smoothing','var') || isempty(smoothing)
    smoothing = false;
end

%% Get the data

% Define the list of phases
phaseVarList = {...
    'InstantaneousPhase_L1y','InstantaneousPhase_L2y','InstantaneousPhase_L3y',...
    'InstantaneousPhase_R1y','InstantaneousPhase_R2y','InstantaneousPhase_R3y'};

% Define the forward speed variable name
velVarList = {'forwardSpeed_mmPerSec'};

% Check whether smoothed variables are desired
if smoothing
    phaseVarList = strcat('smooth_', phaseVarList);
    velVarList   = strcat('smooth_', velVarList);
end

% Extract the needed data
Phi = newData{:, phaseVarList};
v   = newData{:, velVarList};

% Remove rows containing NaN values
isNotNanIdx = ~any(isnan([Phi, v]), 2);
Phi = Phi(isNotNanIdx,:);
v   = v(isNotNanIdx);

% Compute the nine pairwise relative instantaneous phases of interest
Phi_rel = [...
    Phi(:,1) - Phi(:,4),...
    Phi(:,2) - Phi(:,5),...
    Phi(:,3) - Phi(:,6),...  
    Phi(:,2) - Phi(:,1),...
    Phi(:,3) - Phi(:,2),...
    Phi(:,3) - Phi(:,1),...
    Phi(:,5) - Phi(:,4),...
    Phi(:,6) - Phi(:,5)...
    Phi(:,6) - Phi(:,4)...
    ];

% Group the relative phases based on symmetry
Phi_rel = [ Phi_rel(:,1), Phi_rel(:,2), Phi_rel(:,3), Phi_rel(:,4), Phi_rel(:,5), Phi_rel(:,6);...
        Phi_rel(:,1), Phi_rel(:,2), Phi_rel(:,3), Phi_rel(:,7), Phi_rel(:,8), Phi_rel(:,9)];

% Wrap the relative phases to [0, 2*pi)
Phi_rel = mod(Phi_rel, 2*pi);

% Get the bin edges
queryAngles = linspace(0, 2*pi, nbins);

% Define the series labels
seriesLabels = {'L1-R1','L2-R2','L3-R3','M-F','H-M','H-F'};

%% Plot the distributions without conditioning on forward walking speed

% Compute the distributions
N = zeros(length(queryAngles), size(Phi_rel, 2));

for ind=1:size(Phi_rel,2)
    N(:,ind) = vmksdensity(Phi_rel(:,ind), queryAngles);
end

figure('Position',[200,500,1000,1000],'WindowStyle','docked');
ax = polaraxes;
hold on;
set(gca, 'ColorOrder', corder_pairwise);

polarplot(queryAngles, N, 'linewidth', 2);

configurePolarHistogramAxes( ax, '\Delta\phi (cycles)','pdf[ \Delta\phi ] (cycles^{-1})', rlims );

legend(seriesLabels, 'FontSize', 20);

%% Plot the distributions, conditioning on forward walking speed

if ~ suppressExtraPlots
    
    % Compute histograms
    N = zeros(length(queryAngles), size(Phi_rel,2), length(vfBinEdges) + 1);
    
    % v < vfBinEdges(1) mm/s
    idx = v < vfBinEdges(1);
    for ind=1:size(Phi_rel,2)
        N(:, ind, 1) = vmksdensity(Phi_rel(idx,ind), queryAngles);
    end
    
    % v > vfBinEdges(end) mm/s
    idx = v > vfBinEdges(end);
    for ind=1:size(Phi_rel,2)
        N(:, ind, end) = vmksdensity(Phi_rel(idx,ind), queryAngles);
    end
    
    for ind = 2:length(vfBinEdges)
        % vfBinEdges(ind-1) <= v < vfBinEdges(ind) mm/s
        idx = (v>=vfBinEdges(ind-1)) & (v<vfBinEdges(ind));
        for ii=1:size(Phi_rel,2)
            N(:, ii, ind) = vmksdensity(Phi_rel(idx,ii), queryAngles);
        end
    end
    
    for ind=1:size(N,3)
        
        figure('Position',[200,500,1000,1000],'WindowStyle','docked');
        ax = polaraxes;
        hold on;
        set(gca, 'ColorOrder', corder_pairwise);
        
        polarplot(queryAngles, squeeze(N(:,:,ind)), 'linewidth', 2);
        
        if ind==1
            rlabel = sprintf('pdf[ \\Delta\\phi | v_{||} < %d mm/s ] (cycles^{-1})', vfBinEdges(ind));
        elseif ind == size(N,3)
            rlabel = sprintf('pdf[ \\Delta\\phi | v_{||} > %d mm/s ] (cycles^{-1})', vfBinEdges(ind-1));
        else
            rlabel = sprintf('pdf[ \\Delta\\phi | %d \\leq v_{||} < %d mm/s ] (cycles^{-1})', vfBinEdges(ind-1), vfBinEdges(ind));
        end
        
        configurePolarHistogramAxes( ax, '\Delta\phi (cycles)', rlabel, rlims );
        
        legend(seriesLabels, 'FontSize', 20);
    end
    
end

%% Plot polar histograms of all pairwise relative phases

if ~ suppressExtraPlots
    
    seriesLabels = {'L1','L2','L3','R1','R2','R3'};
    
    % Compute the histograms
    idx = 1:6;
    N = zeros(length(queryAngles), 5, 6);
    for ind=1:6
        Phi_rel = mod(Phi(:, idx(idx~=ind) ) - Phi(:, ind), 2*pi);
        for ii=1:5
            N(:, ii, ind) = vmksdensity(Phi_rel(:, ii), queryAngles);
        end
    end
    
    % Plot six separate figures
    for ind=1:6
        figure('Position',[200,500,1000,1000],'WindowStyle','docked');
        ax = polaraxes;
        hold on;
        set(gca, 'ColorOrder', corder_all(idx(idx~=ind), :));
        polarplot(queryAngles, squeeze(N(:,:,ind)), 'linewidth', 2);
        
        configurePolarHistogramAxes( ax, sprintf('Phase relative to %s (cycles)', seriesLabels{ind}), 'pdf[ \Delta\phi ]', [0 max(N(:))] );
        
        legend(seriesLabels(idx(idx~=ind)), 'FontSize', 20);
    end
    
    % Plot everything as six subplots in one summary figure
    figure('Position',[200,500,1000,1000],'WindowStyle','docked');
    for ind=1:6
        ax = subplot(2,3,ind,polaraxes);
        hold on;
        set(gca, 'ColorOrder', corder_all(idx(idx~=ind), :));
        
        polarplot(queryAngles, squeeze(N(:,:,ind)), 'linewidth', 2);
        
        %     configurePolarHistogramAxes( ax, sprintf('Phase relative to %s (cycles)', seriesLabels{ind}), 'pdf[ \Delta\phi ]', [0 max(N(:))] );
        ax.RAxis.Limits = [0 max(N(:))];
        ax.ThetaAxis.Label.String = '\Delta\phi (cycles)';
        
        title(seriesLabels{ind});
        
        legend(seriesLabels(idx(idx~=ind)), 'FontSize', 10, 'Location','South');
        
        % Temporary format change
        ax.FontSize = 20;
        ax.Title.FontSize = 20;
        
    end
    
end
end

function configurePolarHistogramAxes( ax, thetalabel, rlabel, rlims )
ax.ThetaAxisUnits = 'radians';
ax.ThetaZeroLocation = 'right';
ax.ThetaDir = 'counterclockwise';
ax.FontSize = 20;
ax.Title.FontSize = 20;
% ax.FontSize = 40;
% ax.Title.FontSize = 40;
ax.LineWidth = 2;
ax.RAxis.Label.String = rlabel;

if ~isempty(rlims)
    ax.RAxis.Limits = rlims;
end

ax.ThetaAxis.Label.String = thetalabel;
ax.ThetaTick = (0:0.1:1).*(2*pi);
ax.ThetaTickLabel = cellstr(num2str((0:0.1:1)', '%.1f'));
ax.RAxisLocation = 0;
end

