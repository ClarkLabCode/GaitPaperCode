function PlotRelativeInstantaneousPhaseVersusCentroidVelocityComponents( newData, vfBinEdges, vrBinEdges, nbins, cmap, logscale, corder, linewidth, smoothing, showDistributions, suppressExtraPlots)
%PlotRelativeInstantaneousPhaseVersusCentroidVelocityComponents: A function
%to plot distributions of nearest-neighbor relative phases conditioned on
%forward velocity or angular velocity, as well as the circular statistics
%of those distributions.

%% Validate inputs

if ~exist('vfBinEdges','var') || isempty(vfBinEdges)
    vfBinEdges = 5:1:35;
end
vfBinCenters = vfBinEdges(1:end-1) + diff(vfBinEdges)/2;

if ~exist('vrBinEdges','var') || isempty(vrBinEdges)
    vrBinEdges = 0:10:350;
end
vrBinCenters = vrBinEdges(1:end-1) + diff(vrBinEdges)/2;

if ~exist('nbins','var') || isempty(nbins)
    nbins = 100;
end
phiBinEdges = 0:(1/nbins):1;
phiBinCenters = phiBinEdges(1:end-1) + diff(phiBinEdges)/2;

if ~exist('cmap','var') || isempty(cmap)
    cmap = viridis(256);
end

if ~exist('logscale','var') || isempty(logscale)
    logscale = true;
end

if ~exist('corder','var') || isempty(corder)
%     corder = cbrewer('qual','Dark2', 7);
    corder = linspecer(9);
end

if ~exist('linewidth','var') || isempty(linewidth)
    linewidth = 2;
end

if ~exist('smoothing','var') || isempty(smoothing)
    smoothing = true;
end

if ~exist('showDistributions','var') || isempty(showDistributions)
    showDistributions = false;
end

% Define list of series labels
seriesLabels = {'L1-R1','L2-R2','L3-R3','M-F','H-M','H-F'};

%% Extract the data

% Define the list of phases
phaseVarList = {...
    'InstantaneousPhase_L1y','InstantaneousPhase_L2y','InstantaneousPhase_L3y',...
    'InstantaneousPhase_R1y','InstantaneousPhase_R2y','InstantaneousPhase_R3y'};

% Define the list of dynamical variables
velVarList = {'forwardSpeed_mmPerSec','angVel_radPerSec'};

% Check whether smoothed variables are desired
if smoothing
    phaseVarList = strcat('smooth_', phaseVarList);
    velVarList   = strcat('smooth_', velVarList);
end

% Extract the needed data
Phi = newData{:, phaseVarList};
V   = newData{:, velVarList};

% Remove rows containing NaN values
isNotNanIdx = ~any(isnan([Phi, V]), 2);
Phi = Phi(isNotNanIdx,:);
V   = V(isNotNanIdx, :);

% Separate the velocity components
vf = V(:,1);
vr = rad2deg( V(:,2));
vr_abs = abs(vr);
clearvars V;

% Discretize the velocity components into bins
vf_discrete = discretize(vf, vfBinEdges);
vr_discrete = discretize(vr_abs, vrBinEdges);
vf_discrete( vf_discrete == 0 ) = NaN;
vr_discrete( vr_discrete == 0 ) = NaN;

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
vf_discrete = [vf_discrete; vf_discrete]; 
vr_discrete = [vr_discrete; vr_discrete]; 
vf = [vf; vf];
vr = [vr; vr];
vr_abs = [vr_abs; vr_abs];

% Wrap the relative phases to [0, 2*pi)
Phi_rel = mod(Phi_rel, 2*pi);

%% Plot circular statistics binned by forward velocity

% Define the data to be included
X = [Phi_rel, vf_discrete];
X = X(~isnan(vf_discrete), :);

% Define an anonymous function to compute the needed values
% statfun = @(x) [circ_mean(x(:, 1:9),[],1), circ_std(x(:, 1:9), [],[], 1), unique(x(:,end))];
statfun = @(x) [circ_mean(x(:, 1:6),[],1), circ_std(x(:, 1:6), [],[], 1), unique(x(:,end))];

% Get the groups
G = findgroups(X(:,end));

% Apply the anonymous function to each velocity bin
M = splitapply(statfun, X, G);

% Extract data from output
v_bin = vfBinCenters(M(:,end)).';
% mean_delta_phi = mod(M(:,1:9), (2*pi)) / (2*pi);
% sd_delta_phi = M(:,10:18) / (2*pi);
mean_delta_phi = mod(M(:,1:6), (2*pi)) / (2*pi);
sd_delta_phi = M(:,7:12) / (2*pi);
clearvars X M;

% Plot it up
figure('Position',[200,500,1000,1000],'WindowStyle','docked');

subplot(1,2,1);
hold on;
set(gca, 'ColorOrder', corder);
plot(v_bin, mean_delta_phi, 'Linewidth', linewidth);
xlabel('v_{||} (mm/s)');
ylabel('E[ \Delta\phi | v_{||} ] (cycles modulo 1)');
legend(seriesLabels);
ConfAxis;
ylim([0 1]);
xlim([vfBinEdges(1) vfBinEdges(end)]);

subplot(1,2,2);
hold on;
set(gca, 'ColorOrder', corder);
plot(v_bin, sd_delta_phi, 'Linewidth', linewidth);
xlabel('v_{||} (mm/s)');
ylabel('s[ \Delta\phi | v_{||} ] (cycles modulo 1)');
legend(seriesLabels);
ConfAxis;
% ylim([0 .2]);
xlim([vfBinEdges(1) vfBinEdges(end)]);

%% Plot circular statistics binned by absolute angular velocity

if ~suppressExtraPlots
    % Define the data to be included
    X = [Phi_rel, vr_discrete];
    X = X(~isnan(vr_discrete), :);
    
    % Define an anonymous function to compute the needed values
    % statfun = @(x) [circ_mean(x(:, 1:7),[],1), circ_std(x(:, 1:7), [],[], 1), unique(x(:,end))];
    statfun = @(x) [circ_mean(x(:, 1:6),[],1), circ_std(x(:, 1:6), [],[], 1), unique(x(:,end))];
    
    % Get the groups
    G = findgroups(X(:,end));
    
    % Apply the anonymous function to each velocity bin
    M = splitapply(statfun, X, G);
    
    % Extract data from output
    v_bin = vrBinCenters(M(:,end)).';
    mean_delta_phi = mod((M(:,1:6)), (2*pi)) / (2*pi);
    sd_delta_phi = M(:,7:12) / (2*pi);
    
    clearvars X M;
    
    % Plot it up
    figure('Position',[200,500,1000,1000],'WindowStyle','docked');
    
    subplot(1,2,1);
    hold on;
    set(gca, 'ColorOrder', corder);
    plot(v_bin, mean_delta_phi, 'Linewidth', linewidth);
    xlabel('v_{r} (\circ/s)');
    ylabel('E[ \Delta\phi | v_{r} ] (cycles modulo 1)');
    legend(seriesLabels);
    ConfAxis;
    ylim([0 1]);
    
    subplot(1,2,2);
    hold on;
    set(gca, 'ColorOrder', corder);
    plot(v_bin, sd_delta_phi, 'Linewidth', linewidth);
    xlabel('v_{r} (\circ/s)');
    ylabel('s[ \Delta\phi | v_{r} ] (cycles modulo 1)');
    legend(seriesLabels);
    ConfAxis;
    
end
%% Circular statistics conditioned on both forward and rotational velocities

notAbsVrBinEdges = [-fliplr(vrBinEdges(2:end)), vrBinEdges];

cmp = cbrewer('div','Spectral', 2^16);

if ~suppressExtraPlots
    
    MakeFigure;
    for i = 1:3
        subplot(1,3,i);
        colorHistogramByZ(vr, vf, Phi_rel(:,i), @(x) abs(circ_mean(x, [], 1))/ (2*pi), notAbsVrBinEdges, vfBinEdges);
        title(seriesLabels{i});
        
        colormap(cmp);
        cbar = colorbar('southoutside');
        ylabel(cbar, sprintf('E[ \\Delta\\phi_{%s} | v_{r}, v_{||} ] (cycles modulo 1)', seriesLabels{i}));
        cbar.Ticks = 0:0.25:1;
        caxis([0 1]);
        
        xlabel('v_{r} (\circ/s)');
        ylabel('v_{||} (mm/s)');
        axis('xy','square','tight');
        ConfAxis('fontSize', 12);
    end
    
    MakeFigure;
    ind = [1 3 2 4];
    for i = 4:6
        subplot(2,2,ind(i-3));
        colorHistogramByZ(vr, vf, Phi_rel(:,i), @(x) abs(circ_mean(x, [], 1))/ (2*pi), notAbsVrBinEdges, vfBinEdges);
        title(seriesLabels{i});
        
        colormap(cmp);
        cbar = colorbar('eastoutside');
        ylabel(cbar, sprintf('E[ \\Delta\\phi_{%s} | v_{r}, v_{||} ] (cycles modulo 1)', seriesLabels{i}));
        cbar.Ticks = 0:0.25:1;
        caxis([0 1]);
        
        xlabel('v_{r} (\circ/s)');
        ylabel('v_{||} (mm/s)');
        axis('xy','square','tight');
        ConfAxis('fontSize', 12);
    end
    
    cmp = cbrewer('seq','Blues', 2^16);
    
    
    
    MakeFigure;
    for i = 1:3
        subplot(1,3,i);
        colorHistogramByZ(vr, vf, Phi_rel(:,i), @(x) abs(circ_std(x, [], [], 1))/ (2*pi), notAbsVrBinEdges, vfBinEdges);
        title(seriesLabels{i});
        
        colormap(cmp);
        cbar = colorbar('southoutside');
        ylabel(cbar, sprintf('sd[ \\Delta\\phi_{%s} | v_{r}, v_{||} ] (cycles modulo 1)', seriesLabels{i}));
        %    cbar.Ticks = 0:0.25:1;
        caxis([0 1/4]);
        
        xlabel('v_{r} (\circ/s)');
        ylabel('v_{||} (mm/s)');
        axis('xy','square','tight');
        ConfAxis('fontSize', 12);
    end
    
    MakeFigure;
    
    ind = [1 3 2 4];
    for i = 4:6
        subplot(2,2,ind(i-3));
        colorHistogramByZ(vr, vf, Phi_rel(:,i), @(x) abs(circ_std(x, [], [], 1))/ (2*pi), notAbsVrBinEdges, vfBinEdges);
        title(seriesLabels{i});
        
        colormap(cmp);
        cbar = colorbar('eastoutside');
        ylabel(cbar, sprintf('sd[ \\Delta\\phi_{%s} | v_{r}, v_{||} ] (cycles modulo 1)', seriesLabels{i}));
        %    cbar.Ticks = 0:0.25:1;
        caxis([0 1/4]);
        
        xlabel('v_{r} (\circ/s)');
        ylabel('v_{||} (mm/s)');
        axis('xy','square','tight');
        ConfAxis('fontSize', 12);
    end
end

%% Plot distributions of relative phases conditioned on forward velocity

if showDistributions
    
    % Update the series labels for nearest neighbor plots
    seriesLabels = {'L1-R1','L2-R2','L3-R3','L2-L1','L3-L2','R2-R1','R3-R2'};
    
    % Extract the needed data
    Phi = newData{:, phaseVarList};
    V   = newData{:, velVarList};
    
    % Remove rows containing NaN values
    isNotNanIdx = ~any(isnan([Phi, V]), 2);
    Phi = Phi(isNotNanIdx,:);
    V   = V(isNotNanIdx, :);
    
    % Separate the velocity components
    vf = V(:,1);
    vr = rad2deg( V(:,2));
    vr_abs = abs(vr);
    clearvars V;
    
    % Discretize the velocity components into bins
    vf_discrete = discretize(vf, vfBinEdges);
    vr_discrete = discretize(vr_abs, vrBinEdges);
    vf_discrete( vf_discrete == 0 ) = NaN;
    vr_discrete( vr_discrete == 0 ) = NaN;

    
    % Compute the seven nearest-neighbor relative instantaneous phases
    Phi_rel = [...
        Phi(:,1) - Phi(:,4),...
        Phi(:,2) - Phi(:,5),...
        Phi(:,3) - Phi(:,6),...
        Phi(:,2) - Phi(:,1),...
        Phi(:,3) - Phi(:,2),...
        Phi(:,5) - Phi(:,4),...
        Phi(:,6) - Phi(:,5)...
        ];
    
    % Wrap the relative phases to [0, 2*pi) and convert units to cycles/s
    Phi_rel = mod(Phi_rel, 2*pi) / (2*pi);
    
    % Compute the distributions
    N = zeros(length(vfBinCenters), length(phiBinCenters), 7);
    
    isNotNanIdx = ~isnan(vf_discrete);
    uniqueBins  = unique(vf_discrete(isNotNanIdx));
    
    % Iterate over limb pairings
    for ind = 1:7
        % Iterate over velocity bins
        for ii = 1:length(uniqueBins)
            idx = isNotNanIdx & (vf_discrete == uniqueBins(ii));
            N(ii, :, ind) = histcounts( Phi_rel(idx, ind), phiBinEdges, 'normalization','pdf');
        end
    end
    
    % Plot the distributions
    for ind = 1:7
        figure('Position',[200,500,1000,1000],'WindowStyle','docked');
        if logscale
            imagesc( vfBinCenters( uniqueBins ), phiBinCenters, log10( squeeze( N(:,:,ind) ))' );
        else
            imagesc( vfBinCenters( uniqueBins ), phiBinCenters, ( squeeze( N(:,:,ind) ))' );
        end     
        axis('xy','square','tight');
        colormap(cmap);
        cbar = colorbar;
        if logscale
            a = cbar.Ticks;
            cbar.TickLabels = cellstr(num2str(a', '10^{%0.1f}'));
        end
        xlabel('v_{||} (mm/s)');
        ylabel(sprintf('\\Delta_{%s}\\phi (cycles modulo 1)', seriesLabels{ind}));
        ylabel(cbar, sprintf('pdf[ \\Delta_{%s}\\phi | v_{||} ] (cycles^{-2})', seriesLabels{ind}));
        ConfAxis;
    end
    
    %% Plot distributions of relative phases conditioned on angular velocity
    
    % Compute the distributions
    N = zeros(length(vrBinCenters), length(phiBinCenters), 7);
    
    isNotNanIdx = ~isnan(vr_discrete);
    uniqueBins  = unique(vr_discrete(isNotNanIdx));
    
    % Iterate over limb pairings
    for ind = 1:7
        % Iterate over velocity bins
        for ii = 1:length(uniqueBins)
            idx = isNotNanIdx & (vr_discrete == uniqueBins(ii));
            N(ii, :, ind) = histcounts( Phi_rel(idx, ind), phiBinEdges, 'normalization','pdf');
        end
    end
    
    % Plot the distributions
    for ind = 1:7
        figure('Position',[200,500,1000,1000],'WindowStyle','docked');
        if logscale
            imagesc( vrBinCenters( uniqueBins ), phiBinCenters, log10( squeeze( N(:,:,ind) ))' );
        else
            imagesc( vrBinCenters( uniqueBins ), phiBinCenters, ( squeeze( N(:,:,ind) ))' );
        end
        axis('xy','square','tight');
        colormap(cmap);
        cbar = colorbar;
        if logscale
            a = cbar.Ticks;
            cbar.TickLabels = cellstr(num2str(a', '10^{%0.1f}'));
        end
        xlabel('v_{r} (\circ/s)');
        ylabel(sprintf('\\Delta_{%s}\\phi (cycles modulo 1)', seriesLabels{ind}));
        ylabel(cbar, sprintf('pdf[ \\Delta_{%s}\\phi | v_{r} ] (cycles^{-2})', seriesLabels{ind}));
        ConfAxis;
    end
end
end

