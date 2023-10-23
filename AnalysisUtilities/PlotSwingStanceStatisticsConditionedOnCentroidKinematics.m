function PlotSwingStanceStatisticsConditionedOnCentroidKinematics( newData, vfInterval )

fps = 150;
smoothing = true;

%% Prepare data

limbList = {'L1','L2','L3','R1','R2','R3'};

% Define the list of variables to be extracted from the table
tSwingVarList = strcat(limbList, '_SwingDur_mmPerSec');
tStanceVarList = strcat(limbList, '_StanceDur_mmPerSec');
tPeriodVarList = strcat(limbList, '_Period_mmPerSec');

velVarList = {'angVel_radPerSec','forwardSpeed_mmPerSec','translationalSpeed_mmPerSec'};

% Select whether to use smoothed frequency data or not
if smoothing
    velVarList = strcat('smooth_',velVarList);
end

% Extract the needed data
V = newData{:, velVarList};
tSwing = newData{:, tSwingVarList};
tStance = newData{:, tStanceVarList};
tPeriod = newData{:, tPeriodVarList};
videoID = newData.videoID;

% Convert units
V(:,1) = rad2deg(V(:,1));

idx = all((tPeriod < 1e3) & (tPeriod > 0) & (tSwing < 100 | isnan(tSwing)),2);
V = V(idx,:);
tSwing = tSwing(idx,:);
tStance = tStance(idx,:);
tPeriod = tPeriod(idx,:);
videoID = videoID(idx);

% Symmetrize left and right turns
left = V(:,1) < 0;

VSym = V;
tSwingSym = tSwing;
tStanceSym = tStance;
tPeriodSym = tPeriod;

VSym(left, [1,3]) = -VSym(left,[1,3]);

tSwingSym(left,:) = tSwingSym(left, [4,5,6,1,2,3]);
tStanceSym(left,:) = tStanceSym(left, [4,5,6,1,2,3]);
tPeriodSym(left,:) = tPeriodSym(left, [4,5,6,1,2,3]);


%% Conditional averages, with symmetrization

% Swing duration
tic;
seriesLabels = {'O1','O2','O3','I1','I2','I3'};
PlotAverageModulationAsFunctionOfCentroidKinematics(VSym, tSwingSym, videoID,...
'swing duration (ms)',seriesLabels, [0 50], vfInterval );
toc;

% Stance duration
tic;
seriesLabels = {'O1','O2','O3','I1','I2','I3'};
PlotAverageModulationAsFunctionOfCentroidKinematics(VSym, tStanceSym, videoID,...
'stance duration (ms)',seriesLabels, [0 500], vfInterval );
toc;


% Period
tic;
seriesLabels = {'O1','O2','O3','I1','I2','I3'};
PlotAverageModulationAsFunctionOfCentroidKinematics(VSym, tPeriodSym, videoID,...
'period (ms)',seriesLabels, [0 500], vfInterval );
toc;


end