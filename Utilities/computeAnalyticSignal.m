function [ newData ] = computeAnalyticSignal( newData, edgeTrim, centering )
%This function uses the analytic signal method to compute the instantaneous
%amplitudes, phases, and frequencies of the limb data on a per-trajectory
%basis.

% Define number of frames to trim from the ends of each timeseries (to remove edge effects)
if ~exist('edgeTrim', 'var') || isempty(edgeTrim)
    edgeTrim = 3;
end

% Select whether to mean-subtract (centering = false) or center by the
% method suggested by Lamb and Stockl (2014) (centering = true)
if ~exist('centering', 'var') || isempty(centering)
    centering = 0;
end

%% Compute

% Define limb list
limbList = {'L1','L2','L3','R1','R2','R3'};
limbXY = [strcat(limbList, 'x'), strcat(limbList, 'y')];

% Define variable list
varList = [strcat(limbList, '_xPlot_mm'), strcat(limbList, '_yPlot_mm')];

% Extract numerical data from table
data = newData{:, varList};

% Find the groups corresponding to individual trajectories
[G, ~] = findgroups(newData.uniqueFlyTrajID);

% Calculate the analytic signal for each trajectory
[ out ] = splitapply(@(x){trajAnalyticSignal(x, edgeTrim, centering)}, data, G);

% Extract the data
out = cat(1, out{:});
A = cell2mat(out(:,1));
phi = cell2mat(out(:,2));
phiDot = cell2mat(out(:,3));

%% Append data to table

% Append amplitude data to table
amplitudeVarList = strcat('InstantaneousAmplitude_', limbXY);
newData{:, amplitudeVarList} = A;

% Append phase data to table
phaseVarList = strcat('InstantaneousPhase_', limbXY);
newData{:, phaseVarList} = phi;

% Append frequency data to table
frequencyVarList = strcat('InstantaneousFrequency_', limbXY);
newData{:, frequencyVarList} = phiDot;

end

function [ out ] = trajAnalyticSignal( x, edgeTrim, centering )
%This function calculates the instantaneous amplitudes, phases, and
%frequencies for timeseries stored as columns of a data matrix.

% Input checking
if (edgeTrim < 0) || (size(x,1) > edgeTrim)

    % Center the input data
    if centering
        x = x - min(x, [], 1) - (max(x, [], 1) - min(x, [], 1)) / 2;
    else
        x = x - nanmean(x,1);
    end

    % Calculate the analytic signal
    H = hilbert(x);

    % Calculate the instantaneous amplitude
    A = abs(H);

    % Calculate the instantaneous phase
    phi = unwrap(angle(H));

    % Remove ends (to avoid edge effects)
    if edgeTrim > 0
        A(1:edgeTrim,:) = NaN;
        A(end-(edgeTrim-1):end,:)= NaN;

        phi(1:edgeTrim,:) = NaN;
        phi(end-(edgeTrim-1):end,:)= NaN;
    end
else
    phi = nan(size(x));
    A = nan(size(x));
end

% Calculate the instantaneous frequency
phiDot = diff( vertcat(nan(size(phi(1,:))), phi ) );

% Format data for output
out = {A, phi, phiDot};

end
