function [ dataTable, parpoolTocBytes ] = GeneratePhaseOscillatorModelData(tSwing, tStance)
%GeneratePhaseOscillatorModelData: A function to generate data using the phase
%oscillator model from Figure 5 in a table format that can be used with our
%experimental data analysis code. 

%% Parse arguments

% Set default swing duration, in milliseconds
if nargin < 1
    tSwing = 40;
end

% Set default mesh of stance durations, in milliseconds
if nargin < 2
    tStance = 1./linspace(1/40, 1/210, 100)';
end

%% Set fixed parameters

fps = 150;      % Simulated frame rate, in Hz
dt = 1;         % Sampling rate, in frames
dtsub = 0.005;  % Integration timestep, in frames

% Trajectory length (in seconds)
trajLength = 20;

% Time interval to allow relaxation (in seconds)
relaxationTime = 20;

% Cross-body coupling strength
coupling = 1/8;

% Initial condition
y0 = 2*pi*[0.93;0.35;0.20;0.25;0.62;0.47];

% Limb x-position offsets
xOffset = [-3 -2 -1 1 2 3];

%% Prepare parameters for integration

% Convert units of swing and stance durations from milliseconds to frames
framesPerMs = fps/1000;
tSwingFrames = tSwing * framesPerMs;
tStanceFrames = tStance * framesPerMs;

% Form the parameter structure
param.Tswing = tSwingFrames;
param.coupling = coupling;
param.Tstance = tStanceFrames(1);

%% Set up a parallel pool

poolobj = gcp('nocreate');
if isempty(poolobj)
    poolobj = parpool('local');
end

%% Integrate the model

% Time interval over which to integrate (in frames)
tlims = [-relaxationTime, trajLength] * fps;

% Define temporal vector
t = (tlims(1):dt:tlims(2))';

% Compute index to use when discarding relaxation window
notDiscardIdx = (t>=0);
nt = nnz(notDiscardIdx);

% Allocate containers
numT = length(tStance);
theta = nan(nt, 6, numT);

% Define display width
dw = floor(log10(numT))+1;

% Begin monitoring data transfer to and from workers
ticBytes(poolobj);

% Iterate over stance durations
tall = tic;
parfor ind = 1:numT
    % Start a timer for this integration
    tic;

    % Get the current value for the stance duration
    p = param;
    p.Tstance = tStanceFrames(ind);

    % Integrate using a fixed-step explicit Runge-Kutta method with
    % subsampling
    temp = ode2subsample(@(t,y)updatemodel(t,y,p),tlims, y0, dt, dtsub);

    % Discard relaxation window
    theta(:,:,ind) = temp(notDiscardIdx,:,:);
    
    % Print a status update
    fprintf('Integrated parameter set %*d of %*d in %f seconds.\n', dw, ind, dw, numT, toc);
end

% Capture parallel pool data transfer information
parpoolTocBytes = tocBytes(poolobj);

% Print a final status update
fprintf('Completed integration in %f seconds\n', toc(tall));

%% Compute needed quantities
tic;

% Compute limb y-positions
y = cos(theta);

% Compute limb swing/stance
down = mod(theta,2*pi)>pi;

% Define limb x-positions
x = xOffset + zeros(nt, 6, numT);

% Define the fly ID
flyID = kron((1:numT)', ones(nt, 1));

% Define the frame number
frame = (0:nt*numT-1)';

% Define the video id
videoID = ones(numT*nt, 1);

% Define inverse period (frequency) array in Hz
f = 1000./(tSwing(:) + tStance(:));
freq = kron(f, ones(nt,1));

% Define the forward speed
vf = (25-5)/(12.5-4) * (freq-2);

fprintf('Completed data preprocessing in %f seconds.\n', toc);

%% Format data into a table for output
tic;

% Reshape data arrays
theta   = reshape(permute(theta,   [1 3 2]), nt*numT, 6);
y       = reshape(permute(y,       [1 3 2]), nt*numT, 6);
x       = reshape(permute(x,       [1 3 2]), nt*numT, 6);
down    = reshape(permute(down,    [1 3 2]), nt*numT, 6);

% Define variable lists
idVarList = {'Frame','flyID', 'uniqueFlyTrajID', 'videoID'};
velVarList = {'forwardSpeed_mmPerSec','linVel_mmPerSec'};
limbList = {'L1','L2','L3','R1','R2','R3'};
phaseVarListX   = strcat('InstantaneousPhase_', limbList, 'x');
phaseVarListY   = strcat('InstantaneousPhase_', limbList, 'y');
rawPhaseVarList   = strcat(limbList,'_Phase');
limbVarListX   = strcat(limbList, '_xPlot');
limbVarListY   = strcat(limbList, '_yPlot');
limbVarListXmm = strcat(limbList, '_xPlot_mm');
limbVarListYmm = strcat(limbList, '_yPlot_mm');
downVarList    = strcat(limbList, '_down_cam');

% Define list of nonsense values required for compatibility 
nullVarList = {'angVel_radPerSec','translationalSpeed_mmPerSec', 'xCOM','yCOM', 'Orient_Rad_FullRotation'};

% Generate zero array for nonsense values
nullArray = zeros(nt*numT, length(nullVarList));

% Define the variable list
varList = [
    idVarList,velVarList,...
    limbVarListX, limbVarListY, limbVarListXmm, limbVarListYmm,...
    rawPhaseVarList, phaseVarListX, phaseVarListY, downVarList, ...
    {'Frequency_Hz'}, nullVarList
    ];

% Form the data into a homogeneous numeric array
dataArray = [frame, flyID, flyID, videoID, vf, vf, x, y, x, y, zeros(numT*nt, 6), theta, theta, down, freq, nullArray];

% Convert the homogeneous array into a table
dataTable = array2table(dataArray, 'VariableNames', varList);

% Append additional nonsense values for compatability
dataTable.xCOM = zeros(height(dataTable),1);
dataTable.yCOM = zeros(height(dataTable),1);
dataTable.Orient_Rad_FullRotation = zeros(height(dataTable),1);
dataTable.linVel_mmPerSec = dataTable.forwardSpeed_mmPerSec;

fprintf('Converted data into output table format in %f seconds.\n', toc);

end

function dydt = updatemodel(t,y,param)
% Define the system of ODEs

% Extract the parameters from the structure
Tswing = param.Tswing;
Tstance = param.Tstance;
coupling = param.coupling;

% L1 to L3
dydt(1,1) = pi/Tswing*fswing(y(1))/(1+fswing(y(2)) + coupling*sin(y(4)-y(1))) + pi/Tstance*fstance(y(1));
dydt(2,1) = pi/Tswing*fswing(y(2))/(1+fswing(y(3)) + coupling*sin(y(5)-y(2))) + pi/Tstance*fstance(y(2));
dydt(3,1) = pi/Tswing*fswing(y(3))/(1              + coupling*sin(y(6)-y(3))) + pi/Tstance*fstance(y(3));

% R1 to R3
dydt(4,1) = pi/Tswing*fswing(y(4))/(1+fswing(y(5)) + coupling*sin(y(1)-y(4))) + pi/Tstance*fstance(y(4));
dydt(5,1) = pi/Tswing*fswing(y(5))/(1+fswing(y(6)) + coupling*sin(y(2)-y(5))) + pi/Tstance*fstance(y(5));
dydt(6,1) = pi/Tswing*fswing(y(6))/(1              + coupling*sin(y(3)-y(6))) + pi/Tstance*fstance(y(6));

end

function out = fswing(th)
% Indicator function for swing
out = mod(th,2*pi)<pi;
end

function out = fstance(th)
% Indicator function for stance
out = mod(th,2*pi)>pi;
end

function [Y,t] = ode2subsample(odefun,tspan, y0, dt, dts, varargin)
%ODE2SUBSAMPLE  Solve differential equations with a non-adaptive method of
%order 2 with a smaller timestep than the timestep at which samples are
%stored (for the sake of memory efficiency).
% Adapted from ODE2 in the original MATLAB suite of non-adaptive ODE
% integrators by JZV
% See https://blogs.mathworks.com/cleve/2014/05/12/ordinary-differential-equation-suite/

%ODE2  Solve differential equations with a non-adaptive method of order 2.
%   Y = ODE2(ODEFUN,TSPAN,Y0) with TSPAN = [T1, T2, T3, ... TN] integrates 
%   the system of differential equations y' = f(t,y) by stepping from T0 to 
%   T1 to TN. Function ODEFUN(T,Y) must return f(t,y) in a column vector.
%   The vector Y0 is the initial conditions at T0. Each row in the solution 
%   array Y corresponds to a time specified in TSPAN.
%
%   Y = ODE2(ODEFUN,TSPAN,Y0,P1,P2...) passes the additional parameters 
%   P1,P2... to the derivative function as ODEFUN(T,Y,P1,P2...). 
%
%   This is a non-adaptive solver. The step sequence is determined by TSPAN
%   but the derivative function ODEFUN is evaluated multiple times per step.
%   The solver implements the improved Euler (Heun's) method of order 2.   
%
%   Example 
%         tspan = 0:0.1:20;
%         y = ode2(@vdp1,tspan,[2 0]);  
%         plot(tspan,y(:,1));
%     solves the system y' = vdp1(t,y) with a constant step size of 0.1, 
%     and plots the first component of the solution.   
%

if ~isnumeric(tspan)
    error('TSPAN should be a vector of integration steps.');
end

if ~isnumeric(y0)
    error('Y0 should be a vector of initial conditions.');
end


try
    f0 = feval(odefun,tspan(1),y0,varargin{:});
catch ME
    msg = ['Unable to evaluate the ODEFUN at t0,y0. ',ME.message];
    error(msg);
end

y0 = y0(:);   % Make a column vector.
if ~isequal(size(y0),size(f0))
    error('Inconsistent sizes of Y0 and f(t0,y0).');
end

if rem(dt, dts) ~= 0
    error('Subsampling timestep must evenly divide sampling timestep');
end

neq = length(y0);

t = (tspan(1):dt:tspan(2))';
N = length(t);
Y = zeros(neq,N);
F = zeros(neq,2);

% Store the initial condition
Y(:,1) = y0;

% Initialize the current state
yi = y0;

% Iterate over sampling steps
for ii = 2:N
    % Get the time
    ti = t(ii-1);
    
    % perform sub-iterations
    while ti < t(ii)
        F(:,1) = feval(odefun,ti,yi,varargin{:});
        F(:,2) = feval(odefun,ti+dts,yi+dts*F(:,1),varargin{:});
        yi = yi + (dts/2)*(F(:,1) + F(:,2));
        ti = ti + dts;
    end
    
    % Store the value
    Y(:,ii) = yi;
end
Y = Y.';
end