function PlotRelativeInstantaneousPhaseJointDistributions_nonSymmetrized( newData, nbins, cmap, logscale, clims, vfBinEdges, smoothing, pointColor, textColor, suppressExtraPlots )
%PlotRelativeInstantaneousPhaseJointDistributions: A function to plot joint
%distributions of relative instantaneous phase in false color.

%% Validate inputs

if ~exist('nbins','var') || isempty(nbins)
    nbins = 100;
end

if ~exist('cmap','var') || isempty(cmap)
    cmap = viridis(256);
end

if ~exist('logscale','var') || isempty(logscale)
    logscale = true;
end

if ~exist('clims','var') || isempty(clims)
    clims = [];
end

if ~exist('vfBinEdges','var') || isempty(vfBinEdges)
    vfBinEdges = [0, 10, 20, 30];
end

if ~exist('smoothing','var') || isempty(smoothing)
    smoothing = false;
end

if ~exist('pointColor','var') || isempty(pointColor)
    pointColor = 'r';
end

if ~exist('textColor','var') || isempty(textColor)
    textColor = 'w';
end

plotStyle = 'contour';

if ~exist('textColor','var') || isempty(textColor)
    switch plotStyle
        case 'imagesc'
            textColor = 'w';
            
        case 'contour'
            textColor = 'r';
    end
end

contourOpts = {20, 'EdgeColor', 'none'};


%% Get the data

% Define the list of phases
limbVarList = {'L1','L2','L3','R1','R2','R3'};
phaseVarList = strcat('InstantaneousPhase_', limbVarList, 'y');

% Define the forward speed variable name
velVarList = {'forwardSpeed_mmPerSec'};

% Check whether smoothed variables are desired
if smoothing
    phaseVarList = strcat('smooth_', phaseVarList);
    velVarList   = strcat('smooth_', velVarList);
end

% Extract the needed data
Phi = newData{:, phaseVarList};
vf   = newData{:, velVarList};

% Remove rows containing NaN values
isNotNanIdx = ~any(isnan([Phi, vf]), 2);
Phi = Phi(isNotNanIdx,:);
vf   = vf(isNotNanIdx);

% Compute needed pairings
l2r2 = mod(Phi(:,2) - Phi(:,5), 2*pi)/(2*pi); % L2-R2
l3l1 = mod(Phi(:,3) - Phi(:,1), 2*pi)/(2*pi); % L3-L1

r2l2 = mod(Phi(:,5) - Phi(:,2), 2*pi)/(2*pi); % R2-L2
r3r1 = mod(Phi(:,6) - Phi(:,4), 2*pi)/(2*pi); % R3-R1

% Get the bin edges
binEdges = 0:(1/nbins):1;

% Get the bin centers
binCenters = binEdges(1:end-1) + diff(binEdges)/2;

% Compute the histogram (without symmetrizing)
N_left  = histcounts2(l2r2, l3l1, binEdges, binEdges, 'normalization','pdf');
N_right = histcounts2(r2l2, r3r1, binEdges, binEdges, 'normalization','pdf');

[xx,yy] = meshgrid(binCenters);

if logscale
    N_left = log10(N_left);
    N_right = log10(N_right);
    
    idx = isinf(N_left) | isnan(N_left);
    N_left(idx) = floor(min(N_left(~idx)));
    
    idx = isinf(N_right) | isnan(N_right);
    N_right(idx)=floor(min(N_right(~idx)));
end


%% Plot blank field with canonical gait locations labeled

figure('Position',[200,500,1000,1000],'WindowStyle','docked');

plotCanonicalGaitLocations(pointColor, 'k');

configurePhaseJointDistributionAxisLabels();

if ~isempty(clims)
    caxis(clims);
end

%% Plot distribution with canonical gaits overlaid

figure('Position',[200,500,1000,1000],'WindowStyle','docked');

subplot(1,2,1);

switch plotStyle
    case 'contour'
        contourf(xx,yy,N_left', contourOpts{:});
    case 'imagesc'
        imagesc(binCenters, binCenters, (N_left)');
end

axis('xy','equal','tight');

plotCanonicalGaitLocations(pointColor, textColor);

configurePhaseJointDistributionAxisLabels();

colormap(cmap);
cbar = colorbar('southoutside');

if ~isempty(clims)
    caxis(clims);
end

if logscale
    a = cbar.Ticks;
    cbar.TickLabels = cellstr(num2str(a', '10^{%0.1f}'));
end

ylabel(cbar, 'pdf[ \Delta_{M}\phi, \Delta_{H-F}\phi ] (cycles^{-2})');

title('left');

subplot(1,2,2);
switch plotStyle
    case 'contour'
        contourf(xx,yy,N_right', contourOpts{:});
    case 'imagesc'
        imagesc(binCenters, binCenters, (N_right)');
end

axis('xy','equal','tight');

plotCanonicalGaitLocations(pointColor, textColor);

configurePhaseJointDistributionAxisLabels();

colormap(cmap);
cbar = colorbar('southoutside');

if ~isempty(clims)
    caxis(clims);
end

if logscale
    a = cbar.Ticks;
    cbar.TickLabels = cellstr(num2str(a', '10^{%0.1f}'));
end

ylabel(cbar, 'pdf[ \Delta_{M}\phi, \Delta_{H-F}\phi ] (cycles^{-2})');

title('right');

%% Plot distribution without canonical gaits overlaid

figure('Position',[200,500,1000,1000],'WindowStyle','docked');

subplot(1,2,1);
switch plotStyle
    case 'contour'
        contourf(xx,yy,N_left', contourOpts{:});
    case 'imagesc'
        imagesc(binCenters, binCenters, (N_left)');
end

axis('xy','equal','tight');

configurePhaseJointDistributionAxisLabels();

colormap(cmap);
cbar = colorbar('southoutside');

if ~isempty(clims)
    caxis(clims);
end

if logscale
    a = cbar.Ticks;
    cbar.TickLabels = cellstr(num2str(a', '10^{%0.1f}'));
end

ylabel(cbar, 'pdf[ \Delta_{M}\phi, \Delta_{H-F}\phi ] (cycles^{-2})');

title('left');

subplot(1,2,2);
switch plotStyle
    case 'contour'
        contourf(xx,yy,N_right', contourOpts{:});
    case 'imagesc'
        imagesc(binCenters, binCenters, (N_right)');
end

axis('xy','equal','tight');

configurePhaseJointDistributionAxisLabels();

colormap(cmap);
cbar = colorbar('southoutside');

if ~isempty(clims)
    caxis(clims);
end

if logscale
    a = cbar.Ticks;
    cbar.TickLabels = cellstr(num2str(a', '10^{%0.1f}'));
end

ylabel(cbar, 'pdf[ \Delta_{M}\phi, \Delta_{H-F}\phi ] (cycles^{-2})');

title('right');

%% Compute distributions conditioned on forward walking speed

if ~suppressExtraPlots
    
    % Compute histograms
    N_left = zeros(length(binCenters), length(binCenters), length(vfBinEdges) + 1);
    N_right = zeros(length(binCenters), length(binCenters), length(vfBinEdges) + 1);
    
    % v < vfBinEdges(1) mm/s
    idx = vf < vfBinEdges(1);
    N_left(:,:,1) = histcounts2(l2r2(idx), l3l1(idx), binEdges, binEdges, 'normalization','pdf');
    N_right(:,:,1) = histcounts2(r2l2(idx), r3r1(idx), binEdges, binEdges, 'normalization','pdf');
    
    % v > vfBinEdges(end) mm/s
    idx = vf > vfBinEdges(end);
    N_left(:,:,end) = histcounts2(l2r2(idx), l3l1(idx), binEdges, binEdges, 'normalization','pdf');
    N_right(:,:,end) = histcounts2(r2l2(idx), r3r1(idx), binEdges, binEdges, 'normalization','pdf');
    
    for ind = 2:length(vfBinEdges)
        % vfBinEdges(ind-1) <= v < vfBinEdges(ind) mm/s
        idx = (vf>=vfBinEdges(ind-1)) & (vf<vfBinEdges(ind));
        N_left(:,:,ind) = histcounts2(l2r2(idx), l3l1(idx), binEdges, binEdges, 'normalization','pdf');
        N_right(:,:,ind) = histcounts2(r2l2(idx), r3r1(idx), binEdges, binEdges, 'normalization','pdf');
    end
    
    if logscale
        N_left = log10(N_left);
        N_right = log10(N_right);
        
        idx = isinf(N_left) | isnan(N_left);
        N_left(idx) = floor(min(N_left(~idx)));
        
        idx = isinf(N_right) | isnan(N_right);
        N_right(idx)=floor(min(N_right(~idx)));
    end
    
end
%% Plot distributions conditioned on forward walking speed

if ~suppressExtraPlots
    
    for ind=1:size(N_left,3)
        
        figure('Position',[200,500,1000,1000],'WindowStyle','docked');
        % Left
        subplot(1,2,1);
        
        switch plotStyle
            case 'contour'
                contourf(xx,yy,squeeze(N_left(:,:,ind))', contourOpts{:});
            case 'imagesc'
                imagesc(binCenters, binCenters, squeeze(N_left(:,:,ind))');
        end
        axis('xy','equal','tight');
        
        plotCanonicalGaitLocations(pointColor, textColor);
        
        configurePhaseJointDistributionAxisLabels();
        
        colormap(cmap);
        cbar = colorbar('northoutside');
        
        if ~isempty(clims)
            caxis(clims);
        end
        
        if logscale
            a = cbar.Ticks;
            cbar.TickLabels = cellstr(num2str(a', '10^{%0.1f}'));
        end
        
        if ind == 1
            cbarLabelStr = sprintf('pdf[ \\Delta_{M}\\phi, \\Delta_{H-F}\\phi | v_{||} < %d mm/s ] (cycles^{-2})', vfBinEdges(1));
        elseif ind == size(N_left,3)
            cbarLabelStr = sprintf('pdf[ \\Delta_{M}\\phi, \\Delta_{H-F}\\phi | v_{||} > %d mm/s ] (cycles^{-2})', vfBinEdges(ind-1));
        else
            cbarLabelStr = sprintf('pdf[ \\Delta_{M}\\phi, \\Delta_{H-F}\\phi | %d \\leq v_{||} < %d mm/s ] (cycles^{-2})', vfBinEdges(ind-1), vfBinEdges(ind));
        end
        ylabel(cbar, cbarLabelStr);
        
        title('left');
        
        % Right
        subplot(1,2,2);
        switch plotStyle
            case 'contour'
                contourf(xx,yy,squeeze(N_right(:,:,ind))', contourOpts{:});
            case 'imagesc'
                imagesc(binCenters, binCenters, squeeze(N_right(:,:,ind))');
        end
        axis('xy','equal','tight');
        
        plotCanonicalGaitLocations(pointColor, textColor);
        
        configurePhaseJointDistributionAxisLabels();
        
        colormap(cmap);
        cbar = colorbar('northoutside');
        
        if ~isempty(clims)
            caxis(clims);
        end
        
        if logscale
            a = cbar.Ticks;
            cbar.TickLabels = cellstr(num2str(a', '10^{%0.1f}'));
        end
        
        if ind == 1
            cbarLabelStr = sprintf('pdf[ \\Delta_{M}\\phi, \\Delta_{H-F}\\phi | v_{||} < %d mm/s ] (cycles^{-2})', vfBinEdges(1));
        elseif ind == size(N_left,3)
            cbarLabelStr = sprintf('pdf[ \\Delta_{M}\\phi, \\Delta_{H-F}\\phi | v_{||} > %d mm/s ] (cycles^{-2})', vfBinEdges(ind-1));
        else
            cbarLabelStr = sprintf('pdf[ \\Delta_{M}\\phi, \\Delta_{H-F}\\phi | %d \\leq v_{||} < %d mm/s ] (cycles^{-2})', vfBinEdges(ind-1), vfBinEdges(ind));
        end
        ylabel(cbar, cbarLabelStr);
        title('right');
    end
    
end
end

function plotCanonicalGaitLocations(pointColor, textColor)
if ~((isstring(textColor) || ischar(textColor)) && strcmpi(textColor,'none'))
    pointOpts = {'.', 'Color',pointColor, 'Linewidth', 10, 'MarkerSize', 30};
    textOpts = {'FontSize', 20, 'Color',textColor, 'HorizontalAlignment','Center'};
    
    hold on;
    plot(1/2, 1/3, pointOpts{:});
    text(1/2, 1/3 - 1/48, 'Wave gait', textOpts{:});
    
    plot(1/3, 2/3, pointOpts{:} );
    text(1/3, 2/3 - 1/48, 'Tetrapod gait', textOpts{:} );
    
    plot(2/3, 2/3, pointOpts{:} );
    text(2/3, 2/3 - 1/48, 'Tetrapod gait', textOpts{:} );
    
    plot(1/2, 1, pointOpts{:} );
    text(1/2, 1 - 1/48, 'Tripod gait', textOpts{:} );
end
axis('xy','equal','tight');
end

function configurePhaseJointDistributionAxisLabels()
xlabel('\Delta_{M}\phi (cycles modulo 1)');
ylabel('\Delta_{H-F}\phi (cycles modulo 1)');
ConfAxis;
axis('equal','tight');
xlim([0, 1]);
ylim([0, 1]);
xticks([0, 1/6, 1/3, 1/2, 2/3, 5/6, 1]);
yticks([0, 1/6, 1/3, 1/2, 2/3, 5/6, 1]);
xticklabels({'0', '1/6', '1/3', '1/2', '2/3', '5/6', '1'});
yticklabels({'0', '1/6', '1/3', '1/2', '2/3', '5/6', '1'});
end