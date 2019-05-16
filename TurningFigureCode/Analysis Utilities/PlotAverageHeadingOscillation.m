function PlotAverageHeadingOscillation(newData)

% Set the forward speed discretization
vfEdges = 0:1:31;
vfCenters = vfEdges(1:end-1);

% Set the number of bootstraps
nboot = 1000;

% Extract the forward speed
% vf = newData.forwardSpeed_mmPerSec;
vf = newData.smooth_forwardSpeed_mmPerSec;

% Compute residual angular velocity from smoothed version
vr = newData.angVel_radPerSec;
vrresid = (vr - newData.smooth_angVel_radPerSec);
vr = rad2deg(vr);
vrresid = rad2deg(vrresid);

% Compute average residual lateral speed from smoothed version
vp = newData.translationalSpeed_mmPerSec;
vpresid = (vp - newData.smooth_translationalSpeed_mmPerSec);

% Extract IDs
id = newData.videoID;

% Discretize the forward speed
d = discretize(vf, vfEdges);
notNan = ~isnan(d);
d = d(notNan); 
vf = vf(notNan);
vr = vr(notNan);
vrresid = vrresid(notNan);
vp = vp(notNan);
vpresid = vpresid(notNan);
id = id(notNan);

%% Plot the mean absolute yaw

% Compute statistics across videos and forward speeds
a = accumarray([id,d], abs(vr), [length(unique(id)), length(vfEdges)-1], @nanmean);

% Compute 95% confidence intervals via bootstrapping
ci = bootci(nboot, {@nanmean, a});
ci = ci';

% Compute the mean over IDs
m = mean(a,1)';

% Plot everything up
MakeFigure;
PlotAsymmetricErrorPatch( vfCenters', m, ci(:,1), ci(:,2), [1 0 0]);
axis('square');
xlabel('v_{||} (mm/s)');
ylabel('|v_{r}| (\circ/s)');
ConfAxis('fontSize', 14);
ylim([0 200]);
yticks(0:50:200);

%% Plot the mean absolute residual yaw

% Compute statistics across videos and forward speeds
a = accumarray([id,d], abs(vrresid), [length(unique(id)), length(vfEdges)-1], @nanmean);

% Compute 95% confidence intervals via bootstrapping
ci = bootci(nboot, {@nanmean, a});
ci = ci';

% Compute the mean over IDs
m = mean(a,1)';

% Plot everything up
MakeFigure;
PlotAsymmetricErrorPatch( vfCenters', m, ci(:,1), ci(:,2), [1 0 0]);
axis('square');
xlabel('v_{||} (mm/s)');
ylabel('|v_{r} - v_{r}^{smoothed}| (\circ/s)');
ConfAxis('fontSize', 14);
ylim([0 200]);
yticks(0:50:200);

%% Plot the mean absolute lateral speed

% Compute statistics across videos and forward speeds
a = accumarray([id,d], abs(vp), [length(unique(id)), length(vfEdges)-1], @nanmean);

% Compute 95% confidence intervals via bootstrapping
ci = bootci(nboot, {@nanmean, a});
ci = ci';

% Compute the mean over IDs
m = mean(a,1)';

% Plot everything up
MakeFigure;
PlotAsymmetricErrorPatch( vfCenters', m, ci(:,1), ci(:,2), [1 0 0]);
axis('square');
xlabel('v_{||} (mm/s)');
ylabel('|v_{\perp}| (mm/s)');
ConfAxis('fontSize', 14);
ylim([0 5]);
yticks(0:1:5);

%% Plot the mean absolute residual yaw

% Compute statistics across videos and forward speeds
a = accumarray([id,d], abs(vpresid), [length(unique(id)), length(vfEdges)-1], @nanmean);

% Compute 95% confidence intervals via bootstrapping
ci = bootci(nboot, {@nanmean, a});
ci = ci';

% Compute the mean over IDs
m = mean(a,1)';

% Plot everything up
MakeFigure;
PlotAsymmetricErrorPatch( vfCenters', m, ci(:,1), ci(:,2), [1 0 0]);
axis('square');
xlabel('v_{||} (mm/s)');
ylabel('|v_{\perp} - v_{\perp}^{smoothed}| (mm/s)');
ConfAxis('fontSize', 14);
ylim([0 5]);
yticks(0:1:5);

end