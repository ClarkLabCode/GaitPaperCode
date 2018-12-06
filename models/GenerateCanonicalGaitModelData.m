function [ dataTable, gaitInd, freq, L, D, parpoolTocBytes ] = GenerateCanonicalGaitModelData(n, tripodFraction)
% This model is a adaptation of the model developed in Righetti and
% Ijspeert (2008) to generate canonical hexapod locomotor gaits
%
% References:
%       [1] Righetti, Ludovic, and Auke Jan Ijspeert. "Pattern generators
%           with sensory feedback for the control of quadruped locomotion."
%           Robotics and Automation, 2008. ICRA 2008. IEEE International
%           Conference on. IEEE, 2008

%% Set parameters

% Number of sample trajectories
if nargin < 1
    n = 1e5;
end

% Fraction of tripod gait
if nargin < 2
    tripodFraction = 1/3;
end

% Sampling rate
fps = 150;

% Trajectory half-window length, in seconds
trajLength = 1;

% Minimum and maximum swing durations (seconds)
minT = 0.04;
maxT = 0.1;

% Signal-to-noise ratio for adding white Gaussian noise to limb timeseries
noiseSNR = 25;

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

%% Configure the integration window

% Trajectory length, in samples
tLength = 2 * trajLength * fps + 1;

% Timestep
dt = 1/(fps);

% Time span
tspan = [-t0-trajLength, trajLength];

% Query times
tq = (-trajLength:dt:trajLength)';

% Final times (in ms)
t = (-winlen:winlen)' / fps * 1000;

%% Set coupling matrices and duty factors for canonical gaits

% Define the labels
gaitList = {'tripod','left tetrapod','right tetrapod', 'wave'};

% Tripd coupling matrix
kTripod = toeplitz([0, -1, 1, -1, 1, -1]);

% Left tetrapod coupling matrix
kLeftTetrapod = [
    0,    -0.5, -0.5, -0.5, 1,    -0.5; ...
    -0.5, 0,    -0.5, -0.5, -0.5, 1; ...
    -0.5, -0.5, 0,    1,    -0.5, -0.5; ...
    -0.5, -0.5, 1,    0,    -0.5, -0.5; ...
    1,    -0.5, -0.5, -0.5, 0,    -0.5; ...
    -0.5, 1,    -0.5, -0.5, -0.5, 0
    ];

% Right tetrapod coupling matrix
kRightTetrapod = [
    0,    -0.5, -0.5, -0.5, -0.5, 1; ...
    -0.5, 0,    -0.5, 1,    -0.5, -0.5; ...
    -0.5, -0.5, 0,    -0.5, 1,    -0.5; ...
    -0.5, 1,    -0.5, 0,    -0.5, -0.5; ...
    -0.5, -0.5, 1,    -0.5, 0,    -0.5; ...
    1,    -0.5, -0.5, -0.5, -0.5, 0
    ];

% Wave coupling matrix
kWave = toeplitz([0, 0.5, -0.5, -1, -0.5, 0.5]);

% Coupling matrices for all gaits
Kall = cat(3,kTripod,kLeftTetrapod,kRightTetrapod,kWave);

% Duty factors for all gaits
duty_factor_all = [1/2, 2/3, 2/3, 5/6];

%% Set frequencies and gait types

% Select canonical gait types based on set fraction of tripod, with
% remainder of points divided evenly between tetrapod and wave
nTripod = round(n*tripodFraction);
nLeftTetrapod = round(n*(1-tripodFraction)/4);
nRightTetrapod = round(n*(1-tripodFraction)/4);
nWave = n - (nTripod + nLeftTetrapod + nRightTetrapod);
gaitInd = [
    ones(nTripod,1);
    2*ones(nLeftTetrapod,1);
    3*ones(nRightTetrapod,1);
    4*ones(nWave,1)
    ];

% Randomly sample swing durations from the uniform distribution on [minT,
% maxT] (in seconds)
T_swing = (maxT-minT)*rand(n,1) + minT;

% Compute the period in seconds
T = T_swing ./ (1 - duty_factor_all(gaitInd)');

% Compute the frequency in Hz
freq = 1 ./ T;

%% Set up a parallel pool

poolobj = gcp('nocreate');
if isempty(poolobj)
    poolobj = parpool('local');
end

%% Generate sythetic trajectories

% Begin monitoring data transfer to and from workers
ticBytes(poolobj);

% Allocate containers
L = nan(tLength, 6, n);
D = nan(tLength, 6, n);

% Display width
ns = floor(log10(n))+1;

% Iterate in parallel over trajectories
tAll = tic;
parfor ind = 1:n
    tic;
    
    % Get the coupling matrix and duty factor for the current trjaectory
    K = Kall(:,:,gaitInd(ind));
    duty_factor = duty_factor_all(gaitInd(ind));
    
    % Calculate angular frequencies
    omega_swing = (pi) / T_swing(ind);
    omega_stance = ((1-duty_factor)/duty_factor)*omega_swing;
    
    % Set the initial conditions
    X0 = randn(12,1);
    
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

% Capture parallel pool data transfer information
parpoolTocBytes = tocBytes(poolobj);

%% Correct limb ordering for the fact that the model is time-reversal symmetric
% We select hind->front ordering post-hoc
tic;

% Compute the mean relative phase of L3 and L1 using the discrete-time
% analytic signal method
phi = squeeze(nanmean(mod(angle(hilbert(L(:,3,:))) - angle(hilbert(L(:,1,:))),2*pi),1));

% Correct left tetrapod data
idx = phi<pi & gaitInd==2;
L(:,:,idx) = L(:,[6,5,4,3,2,1],idx);
D(:,:,idx) = D(:,[6,5,4,3,2,1],idx);

% Correct right tetrapod data
idx = phi<pi & gaitInd==3;
L(:,:,idx) = L(:,[6,5,4,3,2,1],idx);
D(:,:,idx) = D(:,[6,5,4,3,2,1],idx);

% Correct wave data
idx = phi>pi & gaitInd==4;
L(:,:,idx) = L(:,[6,5,4,3,2,1],idx);
D(:,:,idx) = D(:,[6,5,4,3,2,1],idx);

fprintf('Corrected limb ordering in %f seconds.\n', toc);

%% Select only a desired temporal window

Lall = L;
Dall = D;
w = floor(tLength/2)+1;
L = Lall(w-winlen:w+winlen,:,:);
D = Dall(w-winlen:w+winlen,:,:);

%% Add noise to synthetic limb data
tic;

% Store the raw limb data
Lraw = L;

% Add white Gaussian noise
% (due to the limitations of the awgn function, this must be done on a
% per-trajectory basis)
for ind = 1:size(L,3)
    L(:,:,ind) = awgn(Lraw(:,:,ind), noiseSNR);
end

fprintf('Added white Gaussian noise with SNR %d in %f seconds.\n', noiseSNR, toc);

%% Compute needed quantities
tic;

% Length of temporal window
nt = 2*winlen+1;

% Labels of gaits
gaitLabel = gaitList(kron(gaitInd, ones(nt, 1)))';

% Define limb x-positions
x = xOffset + zeros(nt, 6, n);

% Define the fly ID
flyID = kron((1:n)', ones(nt, 1));

% Define the frame number
frame = (0:nt*n-1)';

% Define the video id
videoID = ones(nt*n,1);

% Define inverse period (frequency) array in Hz
f = kron(freq, ones(nt,1));

% Define the forward speed via a simplistic mapping
vf = (25-5)/(12.5-4) * (f-2);

% Set the centroid coordinates to zero

fprintf('Completed data preprocessing in %f seconds.\n', toc);

%% Format data into a table for output
tic;

% Reshape data arrays
y       = reshape(permute(L, [1 3 2]), nt*n, 6);
x       = reshape(permute(x, [1 3 2]), nt*n, 6);
down    = reshape(permute(D, [1 3 2]), nt*n, 6);

% Define variable lists
idVarList = {'Frame','flyID', 'uniqueFlyTrajID', 'videoID'};
velVarList = {'forwardSpeed_mmPerSec', 'linVel_mmPerSec'};
limbList = {'L1','L2','L3','R1','R2','R3'};
limbVarListX   = strcat(limbList, '_xPlot');
limbVarListY   = strcat(limbList, '_yPlot');
limbVarListXmm = strcat(limbList, '_xPlot_mm');
limbVarListYmm = strcat(limbList, '_yPlot_mm');
downVarList    = strcat(limbList, '_down_cam');

% Define list of nonsense values required for compatibility
nullVarList = {'angVel_radPerSec','translationalSpeed_mmPerSec', 'xCOM','yCOM', 'Orient_Rad_FullRotation'};

% Generate zero array for nonsense values
nullArray = zeros(nt*n, length(nullVarList));

% Define the variable list
varList = [
    idVarList,velVarList,...
    limbVarListX, limbVarListY, limbVarListXmm, limbVarListYmm, ...
    downVarList, {'Frequency_Hz'}, nullVarList
    ];

% Form the data into a homogeneous numeric array
dataArray = [frame, flyID, flyID, videoID, vf, vf, x, y, x, y, down, f, nullArray];

% Convert the homogeneous array into a table
dataTable = array2table(dataArray, 'VariableNames', varList);

% Append cell array of gait type labels
dataTable.GaitType = gaitLabel;

fprintf('Converted data into output table format in %f seconds.\n', toc);

end

%% Define the dynamical equations of the model

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
