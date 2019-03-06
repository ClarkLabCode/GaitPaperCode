function [ newData ] = ComputeHaarWaveletInstantaneousPhase( newData )
%COMPUTEHAARWAVELETINSTANTANEOUSPHASE Instantaneous phase computation using Haar wavelet transforms
%
%
%   References:
%       Goswami, J. C., & Hoefel, A. E. (2004). Algorithms for estimating
%       instantaneous frequency. Signal processing, 84(8), 1423-1427. 

%% Compute instantanoues phases and frequencies

% Define limb list
limbList = {'L1','L2','L3','R1','R2','R3'};
limbXY = [strcat(limbList, 'x'), strcat(limbList, 'y')];

% Define variable list
varList = [strcat(limbList, '_xPlot'), strcat(limbList, '_yPlot')];

% Extract numerical data from table
X = newData{:, varList};

% Find the groups corresponding to individual trajectories
[G, ~] = findgroups(newData.uniqueFlyTrajID);

% Start a timer
tic;

% Compute instantaneous phases and amplitudes
out = cell2mat(splitapply(@(x){haar_wavelet_phase(x)}, X, G));
theta = out(:, 1:12);
amp = out(:, 13:24);

% Approximate instantaneous frequencies
dtheta = [nan(1,12); diff(theta,1,1)];

% Print a timing message to the console
fprintf('Computed instantaneous phases and frequencies using a Haar wavelet transform in %f seconds.\n', toc);

%% Append data to table

if any(contains(newData.Properties.VariableNames, 'InstantaneousPhase'))
    warning('Overwriting existing phase and frequency variables.');
end

% Append phase data to table
phaseVarList = strcat('InstantaneousPhase_', limbXY);
newData{:, phaseVarList} = theta;

% Append frequency data to table
frequencyVarList = strcat('InstantaneousFrequency_', limbXY);
newData{:, frequencyVarList} = dtheta;

% Append amplitude data to table 
amplitudeVarList = strcat('InstantaneousAmplitude_', limbXY);
newData{:, amplitudeVarList} = amp;

end

function [ out ] = haar_wavelet_phase(x)
% A local function to compute the instantaneous phase and amplitude using a
% Haar wavelet transform 

% Center the data
x = x - nanmean(x);

% Shift the data as needed
y = circshift(x,-1);

% Compute the orthogonal terms
c = (x + y)/2;
d = (x - y)/2;

% Set the final elements to NaN
c(end,:) = NaN;
d(end,:) = NaN;

% Compute the phase
theta = unwrap(atan2(d,c));

% Compute the amplitude 
amp = hypot(c,d);

% Combine the data for output
out = [theta, amp];

end