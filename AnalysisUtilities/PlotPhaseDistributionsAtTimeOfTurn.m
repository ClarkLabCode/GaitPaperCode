function [fPDF, fCDF, xq, n, Phi, ID] = PlotPhaseDistributionsAtTimeOfTurn(newData, xq, plotCI, smoothing)


% Set the minimum and maximum forward speed
vfMin = 0;
vfMax = Inf;
vrMin = 0;
vrMax = Inf;

% Set the scale factor for kernel density estimation
bqScale = 0.5;

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
vf = newData.smooth_forwardSpeed_mmPerSec;
vr = abs(rad2deg(newData.smooth_angVel_radPerSec));
idx = (vf > vfMin) & (vf < vfMax) & (vr > vrMin) & (vr < vrMax) & (newData.yawExtremum ~= 0);

% Extract the necessary variables
turnDirection = newData.yawExtremum(idx);
Phi = newData{idx, phaseVarList};
Phi = mod(Phi, 2*pi) / (2*pi);
ID = newData.videoID(idx);

% Remove NaN values
notNan = ~any(isnan(Phi),2);
Phi = Phi(notNan,:);
ID = ID(notNan,:);
turnDirection = turnDirection(notNan);

% Fold left and right turns
Phi(turnDirection<0,:) = Phi(turnDirection<0,[4,5,6,1,2,3]);

%% Compute kernel density estimates of PDFs

% Iterate over limbs
fPDF = zeros(length(xq), 6);
for ind = 1:6
    fPDF(:,ind) = bqksdensity(2*pi*Phi(:,ind), 2*pi*xq, bqScale);
end

% Compute the CDF from the PDF
fCDF = cumsum(fPDF,1) ./ sum(fPDF,1);

% Get the sample size
n = size(Phi,3);

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
        ci = bootci(nboot, {@(x) bqksdensity(2*pi*cell2mat(x), 2*pi*xq, bqScale), X},'Options', statset('UseParallel', true));
        ciLowerPDF(:,ii) = ci(1,:);
        ciUpperPDF(:,ii) = ci(2,:);
        fprintf('Computed confidence intervals for limb %d of 6 in %f seconds\n', ii, toc);
    end
    
end

%% Plot PDF in polar coordinates

seriesLabels = {'O1','O2','O3','I1','I2','I3'};

figure('Position',[200,500,500,700],'WindowStyle','docked');
ax = polaraxes();
polarplot(2*pi*[xq;xq(1)], [fPDF; fPDF(1,:)], 'linewidth', 2);
ax.FontSize = 14;
ax.LineWidth = 2;
ax.RAxisLocation = 0;
ax.RAxis.Label.String = 'pdf';
% ax.RAxisLocation = 0;
ax.ThetaAxisUnits = 'radians';
ax.ThetaTick = 2*pi*(0:0.1:1);
ax.ThetaTickLabel = 0:0.1:1;
ax.ThetaAxis.Label.String = 'phase at time of turn (cycles)';
ax.RAxis.Label.String = 'pdf (1/cycles)';
legend(seriesLabels);
ax.RAxis.Limits = [0 0.5];

%% Plot the PDF and CDF in Cartesian coordinates

figure('Position',[200,500,500,700],'WindowStyle','docked');
if plotCI
    PlotAsymmetricErrorPatch(xq, fPDF, ciLowerPDF, ciUpperPDF, corder);
else
    plot(xq, fPDF, 'linewidth', 2);
end
ylabel('pdf (1/cycles)');
xlabel('phase at time of turn (cycles)');
legend(seriesLabels, 'location','eastoutside');
axis('square');
ConfAxis('fontSize', 16);
xlim([0 1]);
xticks(0:0.1:1);
ylim([0 0.5]);
yticks(0:0.1:0.5);

figure('Position',[200,500,500,700],'WindowStyle','docked');
plot(xq, fCDF, 'linewidth', 2);
ylabel('cdf');
xticks((0:0.1:1));
xticklabels(0:0.1:1);
xlabel('phase at time of turn (cycles)');
legend(seriesLabels, 'location','eastoutside');
axis('square');
ConfAxis('fontSize', 16);
xlim([0 1]);
ylim([0 1]);
yticks(0:0.2:1);

end
