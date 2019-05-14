function PlotFeetDownVsLimbPhase( newData, nbins, vfBinEdges, cmap, corder, limbVarType)
% This function calculates the number of feet down as a function of limb
% phase and forward walking speed.

%% Validate inputs

if ~exist('nbins','var') || isempty(nbins)
    nbins = 20;
end

if ~exist('vfBinEdges','var') || isempty(vfBinEdges)
    vfBinEdges = [15, 25];
end

if ~exist('cmap','var') || isempty(cmap)
    cmap = viridis(256);
end

if ~exist('corder','var') || isempty(corder)
%     corder = linspecer(7,'sequential');
    corder = jet(7) .* 0.8;
end

%% Calculate the data for each limb

% Calculate the number of feet down during each fly-frame
limbList = {'L1','L2','L3','R1','R2','R3'};

if strcmp(limbVarType, 'Egocentric')
    % Egocentric frame binaries
    downVarList = strcat(limbList,'_down');
elseif strcmp(limbVarType, 'Camera')
    % Camera frame binaries
    downVarList = strcat(limbList,'_down_cam');
else
    error('Invalid limb variable type.');
end

num_feet = sum(newData{:, downVarList},2);

% Get the centroid velocities associated with each of these stances
vel_parallel = newData.forwardSpeed_mmPerSec;

% Bin the data for each limb 
phaseVarList = strcat('InstantaneousPhase_', limbList, 'y');
phase_data = newData{:,phaseVarList};
phase_data = mod(phase_data, (2*pi))/(2*pi);
nameList = {'Forelimbs','Midlimbs','Hindlimbs'};


% Collapse the phase data and the corresponding variables into three columns representing F,M,H limbs rather than each limb individually
phase_data = [phase_data(:,[1,2,3]);phase_data(:,[4,5,6])];
num_feet = [num_feet; num_feet];
vel_parallel = [vel_parallel; vel_parallel];

% Get the bin edges
binEdges_count = .5:1:6.5;
binEdges_phase = 0:(1/nbins):1;
binEdges_vel = 0.5:1:35.5;

% Get the bin centers
binCenters_phase = binEdges_phase(1:end-1) + diff(binEdges_phase)/2;
binCenters_vel = binEdges_vel(1:end-1) + diff(binEdges_vel)/2;

% Compute the histogram (Without Forward Speed Conditioning)
for n = 1:3
    N(:,:,n) = histcounts2(phase_data(:,n), num_feet, binEdges_phase, binEdges_count, 'normalization','pdf');
end

% Compute the histogram (With Forward Speed Conditioning)
numSpeeds = length(vfBinEdges) + 1;
speeds = [-inf, vfBinEdges, inf];
for m = 1:numSpeeds
    for n = 1:3 % The number of limbs
        
        curPhi = phase_data(:,n);
        curPhi = curPhi(speeds(m) < vel_parallel & vel_parallel <= speeds(m+1));
        curFeet = num_feet(speeds(m) < vel_parallel & vel_parallel <= speeds(m+1));
        N_binned(:,:, n, m) = histcounts2(curPhi, curFeet, binEdges_phase, binEdges_count, 'Normalization', 'pdf');
    end
end

% Compute the histogram where the joint density for each of the individual limbs
for curFeet = 0:6 % Each number of feet down condition
    
    cur_phase = phase_data(num_feet == curFeet,:);
    cur_vel = vel_parallel(num_feet == curFeet);
    
    for limbType = 1:3 % Each limb type analyzed
        N_joint(:,:,curFeet+1,limbType) = histcounts2(cur_phase(:,limbType), cur_vel, binEdges_phase, binEdges_vel, 'normalization','count');
    end
end

% Normalize the data for the joint density plot so that each forward walking speed sums to 1 (Across all phases and number of feet down)
% NOTE: We need a normalization constant for each V_parallel slice through
% the cube. So the normalization vector will be of length size(N_joint,2)
phaseBinWidth = binEdges_phase(2) - binEdges_phase(1);
speedBinWidth = binEdges_vel(2) - binEdges_vel(1);
normConst = sum(sum(N_joint,2),1);
norm_joint = N_joint ./ (normConst * phaseBinWidth * speedBinWidth);

%% Plot the joint densities - Just the midlimbs (We calculate all the others too)

% Midlimbs, 0 Feet Down
makeFigure;
hold on;
imagesc(binCenters_vel, binCenters_phase, norm_joint(:,:,1,2));
[xList, yList] = meshgrid(binCenters_vel, binCenters_phase);
contour(xList, yList, norm_joint(:,:,1,2),'edgecolor', 'k');

title({'Midlimbs','0 Feet Down'});
ylabel({'Phase', '(Cycle Fractions)'});
set(gca,'YDir','normal');
xlabel({'V_{||}', '(mm/s)'});
axis square;
colormap(cmap);
cbar = colorbar;
cbar.Location = 'northoutside';
% if ~isempty(clims)
%     caxis(clims);
% end
ylabel(cbar, 'pdf (s * mm^{-1} * cycles^{-1})');
ConfAxis;

% Midlimbs, 1 Foot Down
makeFigure;
hold on;
imagesc(binCenters_vel, binCenters_phase, norm_joint(:,:,2,2));
[xList, yList] = meshgrid(binCenters_vel, binCenters_phase);
contour(xList, yList, norm_joint(:,:,2,2),'edgecolor', 'k');
title({'Midlimbs','1 Foot Down'});
ylabel({'Phase', '(Cycle Fractions)'});
set(gca,'YDir','normal');
xlabel({'V_{||}', '(mm/s)'});
axis square;
colormap(cmap);
cbar = colorbar;
cbar.Location = 'northoutside';
% if ~isempty(clims)
%     caxis(clims);
% end
ylabel(cbar, 'pdf (s * mm^{-1} * cycles^{-1})');
ConfAxis;

% Midlimbs, 2 Feet Down
makeFigure;
hold on;
imagesc(binCenters_vel, binCenters_phase, norm_joint(:,:,3,2));
[xList, yList] = meshgrid(binCenters_vel, binCenters_phase);
contour(xList, yList, norm_joint(:,:,3,2),'edgecolor', 'k');
title({'Midlimbs','2 Feet Down'});
ylabel({'Phase', '(Cycle Fractions)'});
set(gca,'YDir','normal');
xlabel({'V_{||}', '(mm/s)'});
axis square;
colormap(cmap);
cbar = colorbar;
cbar.Location = 'northoutside';
% if ~isempty(clims)
%     caxis(clims);
% end
ylabel(cbar, 'pdf (s * mm^{-1} * cycles^{-1})');
ConfAxis;

% Midlimbs, 3 Feet Down
makeFigure;
hold on;
imagesc(binCenters_vel, binCenters_phase, norm_joint(:,:,4,2));
[xList, yList] = meshgrid(binCenters_vel, binCenters_phase);
contour(xList, yList, norm_joint(:,:,4,2),'edgecolor', 'k');
title({'Midlimbs','3 Feet Down'});
ylabel({'Phase', '(Cycle Fractions)'});
set(gca,'YDir','normal');
xlabel({'V_{||}', '(mm/s)'});
xlim([5 25]);
axis square;
colormap(cmap);
cbar = colorbar;
cbar.Location = 'northoutside';
% if ~isempty(clims)
%     caxis(clims);
% end
ylabel(cbar, 'pdf (s * mm^{-1} * cycles^{-1})');
ConfAxis;

% Midlimbs, 4 Feet Down
makeFigure;
hold on;
imagesc(binCenters_vel, binCenters_phase, norm_joint(:,:,5,2));
[xList, yList] = meshgrid(binCenters_vel, binCenters_phase);
contour(xList, yList, norm_joint(:,:,5,2),'edgecolor', 'k');
title({'Midlimbs','4 Feet Down'});
ylabel({'Phase', '(Cycle Fractions)'});
set(gca,'YDir','normal');
xlabel({'V_{||}', '(mm/s)'});
xlim([5 25]);
axis square;
colormap(cmap);
cbar = colorbar;
cbar.Location = 'northoutside';
% if ~isempty(clims)
%     caxis(clims);
% end
ylabel(cbar, 'pdf (s * mm^{-1} * cycles^{-1})');
ConfAxis;

% Midlimbs, 5 Feet Down
makeFigure;
hold on;
imagesc(binCenters_vel, binCenters_phase, norm_joint(:,:,6,2));
[xList, yList] = meshgrid(binCenters_vel, binCenters_phase);
contour(xList, yList, norm_joint(:,:,6,2),'edgecolor', 'k');
title({'Midlimbs','5 Feet Down'});
ylabel({'Phase', '(Cycle Fractions)'});
set(gca,'YDir','normal');
xlabel({'V_{||}', '(mm/s)'});
xlim([5 25]);
axis square;
colormap(cmap);
cbar = colorbar;
cbar.Location = 'northoutside';
% if ~isempty(clims)
%     caxis(clims);
% end
ylabel(cbar, 'pdf (s * mm^{-1} * cycles^{-1})');
ConfAxis;

% Midlimbs, 6 Feet Down
makeFigure;
hold on;
imagesc(binCenters_vel, binCenters_phase, norm_joint(:,:,7,2));
[xList, yList] = meshgrid(binCenters_vel, binCenters_phase);
contour(xList, yList, norm_joint(:,:,7,2),'edgecolor', 'k');
title({'Midlimbs','6 Feet Down'});
ylabel({'Phase', '(Cycle Fractions)'});
set(gca,'YDir','normal');
xlabel({'V_{||}', '(mm/s)'});
axis square;
colormap(cmap);
cbar = colorbar;
cbar.Location = 'northoutside';
% if ~isempty(clims)
%     caxis(clims);
% end
ylabel(cbar, 'pdf (s * mm^{-1} * cycles^{-1})');
ConfAxis;

%% Plot the data for each limb (Separate Figures)

% for i = 1:3
% makeFigure;
% 
% imagesc(binCenters_phase, 1:6, (N(:,:,i))');
% title(nameList{i});
% ylabel('Feet in Stance');
% set(gca,'YDir','reverse');
% xlabel({'Phase', '(Cycle Fractions)'});
% axis square;
% colormap(cmap);
% cbar = colorbar;
% if ~isempty(clims)
%     caxis(clims);
% end
% ylabel(cbar, 'pdf[ Stance Feet, \Delta\phi ] (feet x cycles)^{-1}');
% ConfAxis;
% 
% end

%% Plot the data for each limb (Single Figure)

% makeFigure;
% for i = 1:3
%     
%     subplot(1,3,i);
%     imagesc(binCenters_phase, 1:6, (N(:,:,i))');
%     title(nameList{i});
%     ylabel('Feet in Stance');
%     set(gca,'YDir','reverse');
%     xlabel({'Phase', '(Cycle Fractions)'});
%     axis square;
%     colormap(cmap);
%     cbar = colorbar;
%     if ~isempty(clims)
%         caxis(clims);
%     end
%     ylabel(cbar, 'pdf[ Stance Feet, \Delta\phi ] (feet x cycles)^{-1}');
%     ConfAxis;
% end

    
%% Plot a line plot partitioned by walking speed (Single Figure)

legendStr = {'1 foot down','2 feet down', '3 feet down', '4 feet down', '5 feet down','6 feet down'};

makeFigure;
i = 1;
for n = 1:length(nameList) % Loop through the limbs
    for m = 1:numSpeeds % Loop through the speeds
    
        % Create the subplot
        subplot(3,3,i);
        
        % Loop through and plot a line for each number of feet down
        for k = 1:size(N_binned,2) 
            hold on;
            plot(binCenters_phase', N_binned(:,k,n,m), 'linewidth', 2, 'color', corder(k+1,:)); 
        end
        if i == 6
            legend(legendStr);    
        end
        xlabel({'Phase', '(Cycle Fractions)'});
        axis square;
        ylabel('pdf (1/cycles)');
        ConfAxis;
        
        % Increment the counter
        i = i+1;
    end
end

%% Plot a line plot partitioned by walking speed (Single Figure)
% NOTE: This version only includes 3,4,5 feet down components

legendStr = {'3 feet down', '4 feet down', '5 feet down'};

makeFigure;
i = 1;
for n = 1:length(nameList) % Loop through the limbs
    for m = 1:numSpeeds % Loop through the speeds
        
        % Create the subplot
        subplot(3,3,i);
        
        % Loop through and plot a line for each number of feet down
        for k = 1:size(N_binned,2)
            hold on;
            if ismember(k,[3,4,5])
                plot(binCenters_phase', N_binned(:,k,n,m), 'linewidth', 2, 'color', corder(k+1,:));
            end
        end
        if i == 6
            legend(legendStr);
        end
        xlabel({'Phase', '(Cycle Fractions)'});
        axis square;
        ylabel('pdf (1/cycles)');
        ConfAxis;
        
        % Increment the counter
        i = i+1;
    end
end

%% Plot a line plot partitioned by walking speed - Just Midlimbs (For IsoD1 Data) (Single Figure)
% NOTE: This version only includes 3,4,5 feet down components

legendStr = {'3 feet down', '4 feet down', '5 feet down'};

makeFigure;
i = 1;
n = 2; % Fix to midlimbs
for m = 1:numSpeeds % Loop through the speeds
    
    % Create the subplot
    subplot(1,3,i);
    
    % Loop through and plot a line for each number of feet down
    for k = 1:size(N_binned,2)
        hold on;
        if ismember(k,[3,4,5])
            plot(binCenters_phase', N_binned(:,k,n,m), 'linewidth', 2, 'color', corder(k+1,:));
        end
    end
    
    if i == 1
        ylim([0 .12]);
    elseif i == 2
        ylim([0 .75]);
    elseif i == 3
        ylim([0 1.25]);
        legend(legendStr);
    end
    xlabel({'Phase', '(Cycle Fractions)'});
    axis square;
    ylabel('pdf (1/cycles)');
    ConfAxis;
    
    % Increment the counter
    i = i+1;
end

end

