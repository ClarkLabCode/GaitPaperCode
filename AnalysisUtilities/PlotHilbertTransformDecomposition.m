% PlotHilbertTransformDecomposition.m: Examples for the gait paper

%% Generate a plot of the hilbert tranform decomposition

% Clip the first frame of the trajectory, since it contains NaN values
singleFly = singleFly(2:end,:);

% Define a time Variable
time = (([1:size(singleFly,1)]-1)/150)*1000; % in milliseconds. 150 fps is the frame rate

makeFigure;
subplot(3,1,1);
hold on;
plot(time, -singleFly.L2_yPlot_mm, 'linewidth', 2, 'marker', 'none', 'color', 'b');
xlabel('Time (ms)');
ylabel('Limb Position (mm)');
ConfAxis;
xlim([0, time(end)]);

subplot(3,1,2);
hold on;
plot(time, singleFly.InstantaneousAmplitude_L2y, 'linewidth', 2, 'marker', 'none', 'color', 'b');
ylabel('Instantaneous Amplitude (mm)');
xlabel('Time (ms)');
ylim([0 1]);
ConfAxis;
xlim([0, time(end)]);

subplot(3,1,3);
hold on;
plot(time, mod(singleFly.InstantaneousPhase_L2y,(2*pi))/(2*pi), 'linewidth', 2, 'marker', 'none', 'color', 'b');
ylabel('Phase (Cycle Fractions)');
xlabel('Time (ms)');
ConfAxis;
xlim([0, time(end)]);