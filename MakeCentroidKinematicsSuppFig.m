%% Figure 1 - figure supplement 2
% from % "The manifold structure of limb coordination in walking Drosophila"

% Sample centroid trace with and without smoothing
PlotSampleCentroidTraceSmoothing( newData );

% Comparing autocorrelation functions of centroid variables with and without smoothing
corder = [
    0.346666666666667,0.536000000000000,0.690666666666667;
    0.915294117647059,0.281568627450980,0.287843137254902;
    0.441568627450980,0.749019607843137,0.432156862745098
    ];
[ Rraw, t_ms ] = PlotVelocityAutocorrelation(newData, 75, false, 150, corder, false, 1);
[ Rsmoothed, t_ms ] = PlotVelocityAutocorrelation(newData, 75, true, 150, corder, false, 1);

figure('Position',[200,500,1000,1000],'WindowStyle','docked');
hold on;
set(gca, 'colororder', corder);
plot(t_ms, Rsmoothed, '-', 'linewidth', 2);
plot(t_ms, Rraw, '--', 'linewidth', 2);
axis('xy','square');
ConfAxis('fontSize', 14);
xlabel('time (ms)');
ylabel('autocorrelation');
legend({'v_{r}', 'v_{||}','v_{\perp}'});
xticks(0:100:500);

% Average heading oscillation as a function of forward speed
PlotAverageHeadingOscillation(newData);

% Distribution of FWHM of yaw rates around yaw extrema
PlotYawFWHM(newData);

