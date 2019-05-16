%% Set all parameters

% Colormap
[ ~, cmp, ~ ] = MakeTurningPaperColormaps();

% Number of trajectories
n = 1000;

% Sampling rate
fps = 150;

% Trajectory half-window length, in seconds
trajLength = 1;

% Minimum and maximum swing durations (seconds)
minT = 0.04;
maxT = 0.1;

% Signal-to-noise ratio for adding white Gaussian noise to limb timeseries
noiseSNR = [25; 5];

% Amplitude of oscillations
mu = 1;

% Strength of self-coupling
alpha = 50;

% Sharpness of swing-stance transition
a = 1e8;

% Pre-sampling relaxation time
t0 = 20;

% Embedding sample half-window, in frames
winlen = 15;

% Limb x-position offsets
xOffset = [-3 -2 -1 1 2 3];

% Right tetrapod coupling matrix
K = [
    0,    -0.5, -0.5, -0.5, -0.5, 1; ...
    -0.5, 0,    -0.5, 1,    -0.5, -0.5; ...
    -0.5, -0.5, 0,    -0.5, 1,    -0.5; ...
    -0.5, 1,    -0.5, 0,    -0.5, -0.5; ...
    -0.5, -0.5, 1,    -0.5, 0,    -0.5; ...
    1,    -0.5, -0.5, -0.5, -0.5, 0
    ];

% Duty factor
duty_factor = 2/3;

% Number of bins for phase distributions
nbins = 50;

% Smoothing parameters
wlen = 15;
forder = 3;
dorder = 0;

%% Configure the integration window

% Trajectory length, in samples
tLength = 2 * trajLength * fps + 1;

% Timestep
dt = 1/(fps);

% Time span
tspan = [-t0-trajLength, trajLength];

% Query times
tq = (-trajLength:dt:trajLength)';

% Randomly sample swing durations from the uniform distribution on [minT,
% maxT] (in seconds)
T_swing = (maxT-minT)*rand(n,1) + minT;

% Final times (in seconds)
t = (-winlen:winlen)' / fps;

% Compute the period in seconds
T = T_swing ./ (1 - duty_factor');

% Compute the frequency in Hz
freq = 1 ./ T;

% Get the bin edges
binEdges = 0:(1/nbins):1;

% Get the bin centers
binCenters = binEdges(1:end-1) + diff(binEdges)/2;

%% Set up a parallel pool

poolobj = gcp('nocreate');
if isempty(poolobj)
    poolobj = parpool('local');
end

%% Integrate the model

% Allocate containers
L = nan(tLength, 6, n);
D = nan(tLength, 6, n);

% Display width
ns = floor(log10(n))+1;

% Iterate in parallel over trajectories
tAll = tic;
parfor ind = 1:n
    tic;
    
    % Set the initial conditions
    X0 = randn(12,1);
    
    % Calculate angular frequencies
    omega_swing = (pi) / T_swing(ind);
    omega_stance = ((1-duty_factor)/duty_factor)*omega_swing;
    
    % Integrate the 12-ODE model
    dX = @(t,x) canonicalGaitModelDynamicalEquations(t, x, a, alpha, mu, K, omega_swing, omega_stance);
    [tout, Xout] = ode45(dX, tspan, X0);
    
    % Interpolate to uniform dt
    Xinterp = interp1(tout, Xout, tq);
    
    % Extract coxa joint angles
    L(:,:,ind) = Xinterp(:,1:6);
    
    % Extract which limbs are in stance phase
    D(:,:,ind) = Xinterp(:,7:end) > 0;
    
    % Print a status update
    fprintf('\tTrajectory %*d of %*d: %f s\n',  ns, ind, ns, n, toc);
    
end
fprintf('\nGenerated %*d trajectories in %f seconds\n', ns, n, toc(tAll));


%% Correct limb ordering for the fact that the model is time-reversal symmetric
% We select hind->front ordering post-hoc
tic;

% Compute the mean relative phase of L3 and L1 using the discrete-time
% analytic signal method
phi = squeeze(nanmean(mod(angle(hilbert(L(:,3,:))) - angle(hilbert(L(:,1,:))),2*pi),1));

% Correct right tetrapod data
idx = phi<pi;
L(:,:,idx) = L(:,[6,5,4,3,2,1],idx);
D(:,:,idx) = D(:,[6,5,4,3,2,1],idx);

% Truncate to the desired window length
w = floor(tLength/2)+1;
L = L(w-winlen:w+winlen,:,:);
D = D(w-winlen:w+winlen,:,:);

fprintf('Corrected limb ordering in %f seconds.\n', toc);

%% Plot the phase distributions for each SNR

% Compute the projection matrix
[ b ] = localsg( wlen, forder, dorder );

for snrInd = 1:length(noiseSNR)
    
    % Add white Gaussian noise
    % (due to the limitations of the awgn function, this must be done on a
    % per-trajectory basis)
    Lnoise = L;
    for ind = 1:size(L,3)
        Lnoise(:,:,ind) = awgn(L(:,:,ind), noiseSNR(snrInd));
    end
    
    % Estimate instantaneous phases
    phi_raw = unwrap(angle(hilbert(Lnoise)));
    phi_smooth = localsgfilt(phi_raw, b, wlen);
    Phi = reshape(permute(phi_smooth, [1,3,2]), [], 6);
    
    % Compute needed pairings
    l2r2 = mod(Phi(:,2) - Phi(:,5), 2*pi)/(2*pi); % L2-R2
    l3l1 = mod(Phi(:,3) - Phi(:,1), 2*pi)/(2*pi); % L3-L1
    
    r2l2 = mod(Phi(:,5) - Phi(:,2), 2*pi)/(2*pi); % R2-L2
    r3r1 = mod(Phi(:,6) - Phi(:,4), 2*pi)/(2*pi); % R3-R1
    
    % Compute the histogram (without symmetrizing)
    N_left  = histcounts2(l2r2, l3l1, binEdges, binEdges, 'normalization','pdf');
    N_right = histcounts2(r2l2, r3r1, binEdges, binEdges, 'normalization','pdf');
    
    MakeFigure;
    plot(t, Lnoise(:,:,1), 'linewidth', 2);
    legend({'L1','L2','L3','R1','R2','R3'});
    xlabel('Time (s)');
    ylabel('y-position (arb. units)');
    title(sprintf('SNR = %d', noiseSNR(snrInd)));
    ConfAxis('fontSize', 16);
    
    MakeFigure;
    subplot(1,2,1);
    imagesc(binCenters, binCenters, (N_left)');
    axis('xy','equal','tight');
    hold on;
    plot(2/3, 2/3, '.', 'Color','k', 'Linewidth', 10, 'MarkerSize', 30 );
    configurePhaseJointDistributionAxisLabels();
    xlabel('\Delta_{L2-R2}\phi (cycles modulo 1)');
    ylabel('\Delta_{L3-L1}\phi (cycles modulo 1)');
    title(sprintf('SNR = %d', noiseSNR(snrInd)));
    cbar = colorbar;
    ylabel(cbar, 'pdf (1/cycles^2)');
    caxis([0 round(max([N_left(:); N_right(:)]))]);
    colormap(cmp);
    
    subplot(1,2,2);
    imagesc(binCenters, binCenters, (N_right)');
    axis('xy','equal','tight');
    hold on;
    plot(1/3, 2/3, '.', 'Color','k', 'Linewidth', 10, 'MarkerSize', 30 );
    configurePhaseJointDistributionAxisLabels();
    xlabel('\Delta_{R2-L2}\phi (cycles modulo 1)');
    ylabel('\Delta_{R3-R1}\phi (cycles modulo 1)');
    title(sprintf('SNR = %d', noiseSNR(snrInd)));
    cbar = colorbar;
    ylabel(cbar, 'pdf (1/cycles^2)');
    caxis([0 round(max([N_left(:); N_right(:)]))]);
    colormap(cmp);
end

%%

function [ dX ] = canonicalGaitModelDynamicalEquations( t, X, a, alpha, mu, K, omega_swing, omega_stance)
% This model is a adaptation of the model developed in Righetti and
% Ijspeert (2008) to hexapod locomotion.
%
% References:
%       [1] Righetti, Ludovic, and Auke Jan Ijspeert. "Pattern generators
%           with sensory feedback for the control of quadruped locomotion."
%           Robotics and Automation, 2008. ICRA 2008. IEEE International
%           Conference on. IEEE, 2008

% Separate state variables
x = X(1:6);
y = X(7:12);

% Compute frequencies
omega = omega_stance + (omega_swing-omega_stance) ./ (1+exp(a*y));

% Compute other needed quantities
r = alpha*(mu - hypot(x,y).^2);

% Evaluate the dynamical equations
dX = [
    r .*x - omega.*y;
    r .*y + omega.*x + K*y;
    ];

end


function plotCanonicalGaitLocations(pointColor, textColor)
if ~((isstring(textColor) || ischar(textColor)) && strcmpi(textColor,'none'))
    pointOpts = {'.', 'Color',pointColor, 'Linewidth', 10, 'MarkerSize', 30};
    textOpts = {'FontSize', 20, 'Color',textColor, 'HorizontalAlignment','Center'};
    
    hold on;
    plot(1/2, 1/3, pointOpts{:});
    text(1/2, 1/3 - 1/48, 'Wave gait', textOpts{:});
    
    plot(1/3, 2/3, pointOpts{:} );
    text(1/3, 2/3 - 1/48, 'Tetrapod gait', textOpts{:} );
    
    plot(2/3, 2/3, pointOpts{:} );
    text(2/3, 2/3 - 1/48, 'Tetrapod gait', textOpts{:} );
    
    plot(1/2, 1, pointOpts{:} );
    text(1/2, 1 - 1/48, 'Tripod gait', textOpts{:} );
end
axis('xy','equal','tight');
end

function configurePhaseJointDistributionAxisLabels()
xlabel('\Delta_{M}\phi (cycles modulo 1)');
ylabel('\Delta_{H-F}\phi (cycles modulo 1)');
ConfAxis;
axis('equal','tight');
xlim([0, 1]);
ylim([0, 1]);
xticks([0, 1/6, 1/3, 1/2, 2/3, 5/6, 1]);
yticks([0, 1/6, 1/3, 1/2, 2/3, 5/6, 1]);
xticklabels({'0', '1/6', '1/3', '1/2', '2/3', '5/6', '1'});
yticklabels({'0', '1/6', '1/3', '1/2', '2/3', '5/6', '1'});
end


function [ b ] = localsg( wlen, forder, dorder )
% Optimized Savitzky-Golay filter projection matrix computation

% Validate inputs
if (round(wlen)~=wlen) || ~mod(wlen,2)
    error('Window length must be an odd integer.');
end
if forder > wlen - 1
    error('Filter order must be less than the window length minus one.');
end
if dorder > forder
    error('Derivative order must be less than or equal to the filter order.');
end

% Compute the half-window
halfwin = fix((wlen-1)/2);
x = (-halfwin:halfwin)';

% Form the Vandermonde matrix
v = x .^ (0:forder);

% Optimize for dorder = 0
if dorder == 0
    b = (v * (v \ eye(wlen)))';
else
    k = (1:(forder-dorder));
    hx = [(zeros(wlen,dorder)), ones(wlen,1)*prod(1:dorder), x.^k .* (k+1)];
    b = (hx * (v \ eye(wlen)))';
end
end

function [ y ] = localsgfilt( x, b, wlen )
% Apply the Savitzky-Golay filter defined by the projection matrix b to
% data stored as columns of the matrix x. This implementation assumes that
% missing (NaN) values occur only at the beginning and end of timeseries

% Allocate output
y = nan(size(x));
if size(x,1) > wlen
    for ind = 1:size(x,3)
        % Find NaN rows
        xnotnan = ~any(any(isnan(x),2),3);
        if nnz(xnotnan) > wlen
            % Remove NaN rows
            x(:,:,ind) = x(xnotnan,:,ind);
            
            % Compute the transient on
            ybegin = fliplr(b(:,(wlen-1)/2+2:end))' * flipud(x(1:wlen,:,ind));
            
            % Compute the steady state output
            ycenter = filter(b(:,(wlen-1)./2+1), 1, x(:,:,ind), [], 1);
            
            % Compute the transient off
            yend = fliplr(b(:,1:(wlen-1)/2))' * flipud(x(end-(wlen-1):end,:,ind));
            
            % Store output
            y(xnotnan,:,ind) = [ybegin; ycenter(wlen:end,:,:); yend];
        end
    end
end
end

