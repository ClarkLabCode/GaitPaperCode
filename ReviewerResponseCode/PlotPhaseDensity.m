function PlotPhaseDensity( newData, vEdges, pEdges, cmap, numLvl )
% This function generates distributions of relative phases as a function of 
%forward speed

%% Validate inputs

if ~exist('vEdges','var') || isempty(vEdges)
    vEdges = [0:1:35];
end
vBinCenters = vEdges(1:end-1) + diff(vEdges)/2;

if ~exist('pEdges','var') || isempty(pEdges)
    pEdges = [0:.05:1];
end
pBinCenters = pEdges(1:end-1) + diff(pEdges)/2;

if ~exist('cmap','var') || isempty(cmap)
    cmap = viridis(256);
end

if ~exist('numLvl','var') || isempty(numLvl)
    numLvl = 3;
end

%% Calculate values

% Labels for each of the pairwise relationships
labelList = {'L1-R1', 'L2-R2', 'L3-R3',...
             'L2-L1', 'L3-L2', 'L3-L1',...
             'R2-R1', 'R3-R2', 'R3-R1'};

% Extract the instantaneous phases
varList = {'InstantaneousPhase_L1y', 'InstantaneousPhase_L2y', 'InstantaneousPhase_L3y', ...
    'InstantaneousPhase_R1y', 'InstantaneousPhase_R2y', 'InstantaneousPhase_R3y'};
Phi = newData{:,varList};

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

% Convert relative phases to cycle fractions
Phi_rel = mod(Phi_rel,(2*pi)) ./(2*pi);

% Pull out the velocity data
vel = newData.forwardSpeed_mmPerSec;

% Discretize the velocity into groups
vDisc = discretize(vel, vEdges);

% Allocate an array for storing density estimate
N = nan(length(pBinCenters), length(vBinCenters), 9);

% Compute each of the kernel density estimates
for i = 1:size(Phi_rel,2)
    for indV = 1:length(vEdges)-1
        
        % Kernel Density Estimates
        idx = (vDisc==indV) & ~isnan(Phi_rel(:,i));
        N(:,indV, i) = bqksdensity(2*pi*Phi_rel(idx,i), 2*pi*pBinCenters);     
        
    end
    
end

%% Plot
makeFigure;
for i = 1:size(Phi_rel,2)
    subplot(3,3,i);
    hold on;
    imagesc(vBinCenters, pBinCenters, squeeze(N(:,:,i)));
    [xList, yList] = meshgrid(vBinCenters, pBinCenters);
    contour(xList, yList, squeeze(N(:,:,i)),numLvl,'edgecolor', 'k');
    colormap(cmap);
    cbar = colorbar;
    cbar.Location = 'northoutside';
    title(labelList{i});
    ylabel('P(\Delta \phi | v_{||})');
    axis('xy','square');
    xlabel('v_{||} (mm/s)');
    xticks([10,20,30]);
    yticks([0, .5, 1]);
    ConfAxis;
end

end