function PlotYawFWHM(newData)

winlen = 15;
nboot = 1000;

% Extract turning data
[ ID, V ] = ExtractFoldedTurningData( newData, winlen, true, true);
t = (-winlen:winlen)';

% Get the yaw data
vr = squeeze(V(:,1,:))';

% Get the peak yaws
vrMax = vr(:,t==0);

% Get the IDs
id = squeeze(ID(t==0,4,:));
ni = length(unique(id));

% Compute the deviations from half max
distHalfMax = abs(vr-vrMax/2);

% Find the points closest to half max
[~,halfMaxIdxL] = min( distHalfMax(:,t<0), [], 2);
[~,halfMaxIdxR] = min( distHalfMax(:,t>0), [], 2);

% Compute the half max times
tl = t(t<0);
tr = t(t>=0);
halfMaxL = tl(halfMaxIdxL);
halfMaxR = tr(halfMaxIdxR);

% Compute the full width at half maximum
fullWidthHalfMax = halfMaxR - halfMaxL;

% Compute the distribution of FWHMs
n = accumarray([fullWidthHalfMax,id], 1, [length(t),ni], @sum);

% Compute the mean probability
p = sum(n,2) / length(fullWidthHalfMax);

% Compute confidence intervals
ci = bootci(nboot, {@(x) sum(x,1) / sum(x(:)), n'})';

% Plot the distribution of FWHMs
tau = (0:30)/0.15;
figure('Position',[200,500,500,700],'WindowStyle','docked');
plot(tau, p, '-o','linewidth',2);
axis('square');
xlabel('full width at half maximum yaw (ms)');
ylabel('relative frequency');
ConfAxis('fontSize', 16);

% Plot the distribution of FWHMs with error patches
tau = (0:30)'/0.15;
figure('Position',[200,500,500,700],'WindowStyle','docked');
PlotAsymmetricErrorPatch(tau, p, ci(:,1), ci(:,2), lines(1));
hold on;
plot(tau, p, '-o','linewidth',2);
hold off;
axis('square');
xlabel('full width at half maximum yaw (ms)');
ylabel('relative frequency');
ConfAxis('fontSize', 16);

end