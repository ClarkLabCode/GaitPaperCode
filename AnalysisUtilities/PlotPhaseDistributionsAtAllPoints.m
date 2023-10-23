function [fPDF, fCDF, xq, n, Phi, ID] = PlotPhaseDistributionsAtAllPoints(newData, xq, plotCI, smoothing)

% Set the minimum forward speed
vfMin = 1;

% Set the scale factor for kernel density estimation
bqScale = 1;

% Set whether to use smoothed data
if nargin < 4
    smoothing = true;
end

% Set whether to plot confidence intervals
if nargin < 3
    plotCI = false;
end

% Set the number of bootstraps
nboot = 1000;

% Set the color order
corder = lines(6);

%% Extract the data from the table

% Set the list of variables
limbList = {'L1','L2','L3','R1','R2','R3'};
phaseVarList = strcat('InstantaneousPhase_', limbList, 'y');
if smoothing
    phaseVarList = strcat('smooth_', phaseVarList);
end

% Define an indexing vector
idx = (newData.smooth_forwardSpeed_mmPerSec > vfMin);

% Extract the necessary variables
Phi = newData{idx, phaseVarList};
Phi = mod(Phi, 2*pi) / (2*pi);
ID = newData.videoID(idx);

% Remove NaN values
notNan = ~any(isnan(Phi),2);
Phi = Phi(notNan,:);
ID = ID(notNan,:);

%% Compute kernel density estimates of PDFs

% Iterate over limbs
fPDF = nan(length(xq), 6);
for ind = 1:6
    fPDF(:,ind) = bqksdensity(2*pi*Phi(:,ind), 2*pi*xq, bqScale)/(2*pi);
end

% Compute the CDF from the PDF
fCDF = cumsum(fPDF,1) ./ sum(fPDF,1);

% Get the sample size
n = size(Phi,1);

%% If desired, compute bootstrapped confidence intervals on the PDFs

if plotCI
    
    % Get the list of IDs
    idList = unique(ID);
    ciUpperPDF = nan(length(xq),6);
    ciLowerPDF = nan(length(xq),6);
    
    % Iterate over limbs
    for ii = 1:6
        tic;
        % Get the data for each ID
        X = cell(length(idList),1);
        for ind = 1:length(idList)
            X{ind} = Phi(ID == idList(ind),ii);
        end
        
        % Compute 95% CIs by bootstrapping
        ci = bootci(nboot, {@(x) bqksdensity(2*pi*cell2mat(x), 2*pi*xq, bqScale)/(2*pi), X}, 'Options', statset('UseParallel', true));
        ciLowerPDF(:,ii) = ci(1,:);
        ciUpperPDF(:,ii) = ci(2,:);
        fprintf('Computed confidence intervals for limb %d of 6 in %f seconds\n', ii, toc);
    end
    
end

%% Plot PDF in polar coordinates

figure('Position',[200,500,500,700],'WindowStyle','docked');
ax = polaraxes();
polarplot(2*pi*xq,fPDF, 'linewidth', 2);
ax.FontSize = 14;
ax.LineWidth = 2;
ax.RAxisLocation = 0;
ax.RAxis.Label.String = 'pdf';
% ax.RAxisLocation = 0;
ax.ThetaAxisUnits = 'radians';
ax.ThetaTick = 2*pi*(0:0.1:1);
ax.ThetaTickLabel = 0:0.1:1;
ax.ThetaAxis.Label.String = 'phase (cycles)';
ax.RAxis.Label.String = 'pdf (1/cycles)';
legend(limbList);
ax.RAxis.Limits = [0 0.5];
title(sprintf('v_{||} > %d mm/s', vfMin));

%% Plot PDF and CDF in Cartesian coordinates

figure('Position',[200,500,500,700],'WindowStyle','docked');
if plotCI
    PlotAsymmetricErrorPatch(xq, fPDF, ciLowerPDF, ciUpperPDF, corder);
else
    plot(xq, fPDF, 'linewidth', 2);
end
ylabel('pdf (1/cycles)');
xticks(0:0.1:1);
xticklabels(0:0.1:1);
xlabel('phase (cycles)');
legend(limbList, 'location','eastoutside');
axis('square');
ConfAxis('fontSize', 16);
xlim([0 1]);
ylim([0 0.5]);
yticks(0:0.1:0.5);
title(sprintf('v_{||} > %d mm/s', vfMin));

figure('Position',[200,500,500,700],'WindowStyle','docked');
plot(xq, fCDF, 'linewidth', 2);
ylabel('cdf');
xticks(0:0.1:1);
xticklabels(0:0.1:1);
xlabel('phase (cycles)');
legend(limbList, 'location','eastoutside');
axis('square');
ConfAxis('fontSize', 16);
xlim([0 1]);
ylim([0 1]);
yticks(0:0.2:1);
title(sprintf('v_{||} > %d mm/s', vfMin));

end