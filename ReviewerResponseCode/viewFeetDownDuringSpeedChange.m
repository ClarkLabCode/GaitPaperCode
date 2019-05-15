function viewFeetDownDuringSpeedChange( newData )
% This function generates visualizations of the distributions of number of
% feet down in trajectories with accelerations and deccelerations. The
% intention of this analysis is to show that the distributions don't
% change.

% Define the inputs for the extracting example trajectories
winlen = 100; % Dial this up and down
trigger = newData.forwardSpeed_mmPerSec >= 15 & newData.forwardSpeed_mmPerSec <= 20;  
varList = {'forwardSpeed_mmPerSec',...
           'L1_down_cam', 'L2_down_cam', 'L3_down_cam', ...
           'R1_down_cam', 'R2_down_cam', 'R3_down_cam', ...
           'InstantaneousPhase_L1y', 'InstantaneousPhase_L2y', 'InstantaneousPhase_L3y', ...
           'InstantaneousPhase_R1y', 'InstantaneousPhase_R2y', 'InstantaneousPhase_R3y'};

% Extract the individual trajectories in a cell array 
[D] = getSlicesOneSidedEither(newData, varList, winlen, trigger, 'pre');

% Convert the data cell array into a tensor with dimensions (time, vars, examples) D
D = cat(3,D{:});

% Calculate the mean acceleration over each trajectory
spf = 1/150; % Seconds per Frame
v = squeeze(D(:,1,:));
a = (v(end,:) - v(1,:)) ./ (size(v,1)*spf);

% Calculate the range of speeds for each trajectory
v_range = max(v,[],1) - min(v,[],1);

%% Plot 1 - Distributions of accelerations and decelerations

makeFigure;
plot(v(end,:), a, 'linestyle', 'none', 'marker', '.');
xlabel('v_{||}');
ylabel('dv/dt');
title('Acceleration vs. Forward Speed');
ConfAxis;

% TODO: Convert to a line plot
% Distribution of average accelerations
makeFigure;
histogram(a);
ylabel('Count');
xlabel('dv/dt');
ConfAxis;

% TODO: Convert to a line plot
% Plot of range of velocities
makeFigure;
histogram(v_range);
ylabel('Count');
xlabel('Velocity Range (mm/s)');
title('Max Velocity - Min Velocity');
ConfAxis;

%% Partition the dataset to look at trajectories where the fly speeds up and slows down at the same speeds
% TODO: Convert these to percentiles

% Calculate the instantaneous number of feet down in each condition
D(:,14,:) = sum(D(:,2:7,:),2);

% Isolate trajectories where the fly speeds up
% upIdx = a > 20;
upIdx = a > 10 & a <15;
up = D(:,:,upIdx);

% Isolate trajectories where the fly slows down
% downIdx = a < -15;
downIdx = a > -15 & a < -10;
down = D(:,:,downIdx);

%% Plot 2 - Time-series of the selected trajectories

% Compute the mean and SEM trajectories
% Accelerations
vUp = squeeze(up(:,1,:));
mean_vUp = mean(vUp,2);
sem_vUp = std(vUp, [], 2) ./ size(vUp,2)^(1/2);

% Decelerations
vDown = squeeze(down(:,1,:));
mean_vDown = mean(vDown,2);
sem_vDown = std(vDown, [], 2) ./ size(vDown,2)^(1/2);

% Compute the standard error of the mean for each trajectory
time = [-winlen:0]' ./ 150;
MakeFigure;
PlotXvsY([time,time],[mean_vUp, mean_vDown], 'error', [sem_vUp, sem_vDown]);
ylim([0 28]);
ylabel('Forward Speed (mm/s)');
xlabel('Time (s)');
ConfAxis;


%% Plot 3 - Distribution of footfalls

% Extract out the number of feet down at t=0 from each of 
cntUp = squeeze(up(end,14,:));
cntDown = squeeze(down(end,14,:));

% Calculate the number of feet down distributions for each of these conditions
pdfUp = histcounts(cntUp, [-.5:1:6.5],'normalization','probability');
pdfDown = histcounts(cntDown, [-.5:1:6.5],'normalization','probability');

% Plot the results
makeFigure; 
plot([0:6;0:6]', [pdfUp;pdfDown]', 'linewidth', 2);
xlabel('Number of Feet in Stance');
ylabel('Probability');
legend({'Accelerating', 'Decelerating'});
ConfAxis;

%% Plot 4 - Distribution of limb relative phases

% Extract out the limb phases at time t=0
Phi = squeeze(D(end,8:13,:))';

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

% Labels for each of the pairwise relationships
labelList = {'L1-R1', 'L2-R2', 'L3-R3',...
             'L2-L1', 'L3-L2', 'L3-L1',...
             'R2-R1', 'R3-R2', 'R3-R1'};

% Convert relative phases to cycle fractions
Phi_rel = mod(Phi_rel,(2*pi)) ./(2*pi);

% Get the relative phases for the up and down conditions
phiUp = Phi_rel(upIdx,:);
phiDown = Phi_rel(downIdx,:);

% Generate the historgrams
pEdges = [0:.01:1];
pBinCenters = pEdges(1:end-1) + diff(pEdges)/2;

% Remomve all rows that contain nanPhases from phiUp and phiDown
phiUp = phiUp(~any(isnan(phiUp),2),:);
phiDown = phiDown(~any(isnan(phiDown),2),:);

% Loop through and generate the distributions for each of the relative phases
for prIdx = 1:9
   
   PUp(:,prIdx) = bqksdensity(2*pi*phiUp(:,prIdx), 2*pi*pBinCenters);
   PDown(:,prIdx) = bqksdensity(2*pi*phiDown(:,prIdx), 2*pi*pBinCenters);
    
end

% Plot each of the relative phase distributions
makeFigure;
for prIdx = 1:9
    
    subplot(3,3,prIdx);
    plot([pBinCenters',pBinCenters'],[PUp(:,prIdx),PDown(:,prIdx)], 'linewidth', 2);
    title(labelList{prIdx});
    ylim([0,1]);
    axis('square','tight');
    ylabel('PDF');
    xlabel('\Delta \phi (Cycle Fractions)');
    ConfAxis;
    
    if prIdx == 9
        legend({'Accelerating', 'Decelerating'});
    end
end

end
