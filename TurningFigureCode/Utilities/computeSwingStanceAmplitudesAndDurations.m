function [ newData ] = computeSwingStanceAmplitudesAndDurations( newData, logicalType )
% This function computes the Euclidean distance in both the camera frame
% and the fly frame associated with a given swing/stance logical.
% Additionally, this function calculates the durations associated with each
% swing and stance event and paired periods of swing followed by stance. 

% NOTE: In generating this view, we can pull/filter all the first frames of 
% each swing and stance event and get a dataset that contains all the swing
% and stance duration and amplitude information desired where the
% underlying data is organized on a per-step basis.

% Define the conversion from pixels to mm 
mmPerPix = 0.043;

% Set default value for second argument
if nargin < 2
    logicalType = '_down_cam';
end

%% Remove trajectories that are shorter than the desired length from the dataset

% Define the rows that need to be trimmed
minTrajLength = 15;
G = findgroups(newData.uniqueFlyTrajID);
tbl = tabulate(G);
trajToTrim = tbl(tbl(:,2) < minTrajLength, 1);
trimmedRows = ~ismember(G, trajToTrim);

% Trim the dataset
newData = newData(trimmedRows,:);

%% Define variable lists

limbList = {'L1','L2','L3','R1','R2','R3'};

% Define variable list for camera frame limb positions
limbVarListXCam = strcat(limbList, '_xCam');
limbVarListYCam = strcat(limbList, '_yCam');

% Define variable list for egocentric frame limb positions
limbVarListXEgo = strcat(limbList, '_xPlot');
limbVarListYEgo = strcat(limbList, '_yPlot');

%% Extract needed data from the table

% Extract camera frame limb positions
xCam = newData{:, limbVarListXCam};
yCam = newData{:, limbVarListYCam};

% Extract egocentric frame limb positions
xEgo = newData{:, limbVarListXEgo};
yEgo = newData{:, limbVarListYEgo};

% Extract out the forward speed of the fly
vf = newData{:,'forwardSpeed_mmPerSec'};

%% Locate the limb extreme positions

% Define the list of variables
% upDownVarList = strcat(limbList, '_down');
% upDownVarList = strcat(limbList, '_down_cam');
upDownVarList = strcat(limbList, logicalType);

% Get the unique trajectories
G = findgroups(newData.uniqueFlyTrajID);

% Get the up/down logicals
D = newData{:, upDownVarList};

% Get the locations of all paired events in the dataset 
% NOTE: These paired events differ based on the 'pair-type'
% 'ud' = swing --> stance
% 'du' = stance --> swing
% 'uu' = swing --> next swing

tic;
% Swing --> Stance
A_ud = splitapply(@(x) {local_get_indices(x, 'ud', minTrajLength)}, D, G);
A_ud = cell2mat(A_ud);
% Stance --> Swing
A_du = splitapply(@(x) {local_get_indices(x, 'du', minTrajLength)}, D, G);
A_du = cell2mat(A_du);
% Swing --> Next Swing
A_uu = splitapply(@(x) {local_get_indices(x, 'uu', minTrajLength)}, D, G);
A_uu = cell2mat(A_uu);
fprintf('Located paired limb events in %f seconds.\n', toc);

%% Compute the values associated with each swing and stance paired event

% Define the Euclidean distance metric
euc = @(p,a) sqrt(dot(p,p,2)+dot(a,a,2)-2*dot(a,p,2));

% Compute the swing amplitudes in the camera frame 
% NOTE: Put in the location of each swing initiation
tic;
dCam = nan(size(A_ud));
for i = 1:6
    p = [xCam(A_ud(:,i)==1,i), yCam(A_ud(:,i)==1,i)];
    a = [xCam(A_ud(:,i)==2,i), yCam(A_ud(:,i)==2,i)];
    dCam(A_ud(:,i)==1,i) = euc(p,a);
end
fprintf('Computed camera frame displacements in %f seconds.\n', toc);

% Compute the swing amplitudes in the egocentric frame
% NOTE: Put in the location of each swing initiation
tic;
dEgo = nan(size(A_ud));
for i = 1:6
    p = [xEgo(A_ud(:,i)==1,i),yEgo(A_ud(:,i)==1,i)];
    a = [xEgo(A_ud(:,i)==2,i),yEgo(A_ud(:,i)==2,i)];
    dEgo(A_ud(:,i)==1,i) = euc(p,a);
end
fprintf('Computed egocentric frame displacements in %f seconds.\n', toc);

% Compute the swing duration
tic; 
tSwing = nan(size(A_ud));
rowID = [1:size(A_ud,1)]';
for i = 1:6
    srt = rowID(A_ud(:,i)==1); % Swing Starts
    stp = rowID(A_ud(:,i)==2); % Stance Starts
    tSwing(A_ud(:,i)==1,i) = stp-srt;
end
fprintf('Computed the swing durations in %f seconds.\n', toc);

% Compute the average forward speed over the swing event
tic; 
vSwing = nan(size(A_ud));
rowID = [1:size(A_ud,1)]';
for i = 1:6
    srt = rowID(A_ud(:,i)==1); % Swing Starts
    stp = rowID(A_ud(:,i)==2); % Stance Starts
    vSwing(A_ud(:,i)==1,i) = arrayfun(@(a,b) mean(vf(a:b)), srt, stp);
end
fprintf('Computed the average forward speed during swing in %f seconds.\n', toc);


%% Compute the values associated with each stance and swing paired event

% Compute the stance duration
tic; 
tStance = nan(size(A_du));
rowID = [1:size(A_du,1)]';
for i = 1:6
    srt = rowID(A_du(:,i)==1); % Stance Starts
    stp = rowID(A_du(:,i)==2); % Swing Starts
    tStance(A_du(:,i)==1,i) = stp-srt;
end
fprintf('Computed the stance durations in %f seconds.\n', toc);

% Compute the average forward speed over the stance event
tic; 
vStance = nan(size(A_du));
rowID = [1:size(A_du,1)]';
for i = 1:6
    srt = rowID(A_du(:,i)==1); % Stance Starts
    stp = rowID(A_du(:,i)==2); % Swing Starts
    vStance(A_du(:,i)==1,i) = arrayfun(@(a,b) mean(vf(a:b)), srt, stp);
end
fprintf('Computed the average forward speed during stance in %f seconds.\n', toc);

%% Compute the values associated with full swing-stance full period event

% Compute the period duration
tic; 
tPeriod = nan(size(A_uu));
rowID = [1:size(A_uu,1)]';
for i = 1:6
    srt = rowID(A_uu(:,i)==1); % Swing Starts
    stp = rowID(A_uu(:,i)==2); % Next Swing Starts
    tPeriod(A_uu(:,i)==1,i) = stp-srt;
end
fprintf('Computed the periods in %f seconds.\n', toc);

%% Convert units to standard units

% Convert amplitude units from pixels to mm
dCam = mmPerPix * dCam;
dEgo = mmPerPix * dEgo;

% Convert duration units from frames to milliseconds
tSwing = tSwing * 100/15;
tStance = tStance * 100/15;
tPeriod = tPeriod * 100/15;

%% Propagate the new variable values and append to the complete dataset

% Create variables that appropriately propogate the values to the desired rows
[ dCamInterp ] = local_last_neighbor_interp(A_ud, G, 'partial', dCam);
[ dEgoInterp ] = local_last_neighbor_interp(A_ud, G, 'partial', dEgo);
[ tSwingInterp ] = local_last_neighbor_interp(A_ud, G, 'partial', tSwing);
[ tStanceInterp ] = local_last_neighbor_interp(A_du, G, 'partial', tStance);
[ vSwingInterp ] = local_last_neighbor_interp(A_ud, G, 'partial', vSwing);
[ vStanceInterp ] = local_last_neighbor_interp(A_du, G, 'partial', vStance);
[ tPeriodInterp ] = local_last_neighbor_interp(A_uu, G, 'all', tPeriod);

% Define some variable names for the new variables of interest
% Data
dCamVarList = strcat(limbList, '_CamStepAmplitude_mm');
dEgoVarList = strcat(limbList, '_EgoStepAmplitude_mm');
tSwingVarList = strcat(limbList, '_SwingDur_mmPerSec');
tStanceVarList = strcat(limbList, '_StanceDur_mmPerSec');
tPeriodVarList = strcat(limbList, '_Period_mmPerSec');
vSwingVarList = strcat(limbList, '_SwingAvgVf_mmPerSec');
vStanceVarList = strcat(limbList, '_StanceAvgVf_mmPerSec');

% Logicals for swing and stance starts
swingVarList = strcat(limbList, '_SwingStart');
stanceVarList = strcat(limbList, '_StanceStart');

% Append the new measurements to the dataset
% Data
newData{:, dCamVarList} = dCamInterp;
newData{:, dEgoVarList} = dEgoInterp;
newData{:, tSwingVarList} = tSwingInterp;
newData{:, tStanceVarList} = tStanceInterp;
newData{:, tPeriodVarList} = tPeriodInterp;
newData{:, vSwingVarList} = vSwingInterp;
newData{:, vStanceVarList} = vStanceInterp;

% Logicals
newData{:, swingVarList} = (A_ud == 1);
newData{:, stanceVarList} = (A_du == 1);

end

function [XInterp ] = local_last_neighbor_interp(A, G, prop_type, X)
% Last-neighbor interpolation

tic;

% Initialize empty vectors of NaNs to store the final values
XInterp = nan(size(A));

% Define temporary 'queue' variables for storing the most recent value of interest
q = nan(1,6);

% Loop through all the calculated values of interest
for i = 1:size(X,1)
    
    % Loop through each limb
    for j = 1:6
        
        % Update the queued values
        switch prop_type
            
            case 'partial'
                % Propogate a value from A == 1 until A == 2
                if A(i,j) == 1
                    q(j) = X(i,j);
                    % Copy over the current value 
                    XInterp(i,j) = q(j);
                elseif A(i,j) == 2
                    q(j) = NaN;
                end
                
            case 'all'
                % Propogate a value from A == 1 until A == 2
                if A(i,j) == 1
                    q(j) = X(i,j);
                    % Copy over the current value 
                    XInterp(i,j) = q(j);
                elseif A(i,j) == 2
                    XInterp(i,j) = q(j); % This line is the only thing distinguishing the 'all' case from 'partial'
                    q(j) = NaN;                    
                end                
            
            otherwise
                error('Invalid propogation type. Must be "partial" or "all"');
                
        end
        
        % If the current value is NaN but the queue is NOT NaN, propagate the non-NaN value
        % NOTE: When this condition is not met, the value remains NaN
        if isnan(X(i,j)) && ~isnan(q(j))
            
            % Store the values
            XInterp(i,j) = q(j);
            
        end        
        
    end
    
    % When the trajectory changes reset the queues to be NaN values
    if (i<size(XInterp,1)) && (G(i+1)~=G(i))
        q = nan(1,6);
    end
end
fprintf('Interpolated values in %f seconds\n', toc);

end

function [ out ] = local_get_indices(x, pair_type, minTrajLength)
% Local function for identifying swing and stance transitions in each
% trajectory.

% Initialize the output variable to be the same size of the input dataset
out = zeros(size(x,1),6);

% Check that the trajectory is long enough to have multiple steps
if size(x,1) >= minTrajLength
    
    % Loop through each of the limbs
    for n = 1:6
        queue = zeros(1,2);
        i = 1;
        
        switch pair_type
            
            % Swing followed by stance
            case 'ud'
                
                op = false; % This tell us which column to populate                
                % Loop through all the frames
                for m = 1:size(x,1)-1
                    
                    % If cur is swing AND next is stance AND we're looking for swing to stance
                    % Then grab the row position of the swing to stance transition
                    if ~x(m,n) && x(m+1,n) && op
                        queue(i,2) = m+1;
                        i = i + 1;
                        op = false;
                        
                        % If cur is stance AND next is swing AND we're looking for stance to swing
                        % Then grab the row position of the stance to swing transition
                    elseif x(m,n) && ~x(m+1,n) && ~op
                        queue(i,1) = m+1;
                        op = true;
                    end
                end
            
            % Stance followed by swing    
            case 'du'
                
                op = true; % This tell us which column to populate                
                % Loop through all the frames
                for m = 1:size(x,1)-1
                    
                    % If cur is swing AND next is stance AND we're looking for swing to stance
                    % Then grab the row position of the swing to stance transition
                    if ~x(m,n) && x(m+1,n) && op
                        queue(i,1) = m+1;
                        op = false;
                        
                        % If cur is stance AND next is swing AND we're looking for stance to swing
                        % Then grab the row position of the stance to swing transition
                    elseif x(m,n) && ~x(m+1,n) && ~op
                        queue(i,2) = m+1;
                        i = i + 1;
                        op = true;
                    end
                end                
            
            % Swing followed by next swing
            case 'uu'
                
                op = false; % This tell us which column to populate                
                % Loop through all the frames
                for m = 1:size(x,1)-1
                    
                    % Grab all subsequent swing starts
                    if x(m,n) && ~x(m+1,n) && op
                        queue(i,1) = prev+1; % First swing frame
                        queue(i,2) = m; % NOTE: Here we grab the position of the LAST stance frame
                        i = i + 1;
                        % Grab the last stance frame prior to the swing event
                        prev = m;
                        
                    % Grab the first stance to swing event
                    elseif x(m,n) && ~x(m+1,n) && ~op
                        prev = m;
                        op = true;
                    end
                end
                
        end
        
        % Grab all the locations in queue that have both a start and an end
        queue = queue(all(queue,2),:);
        
        % Locate all the locations that are the first event (ex. 'ud' --> start of a swing)
        out(queue(:,1),n) = 1;
        
        % Locate all the positions that are the second event (ex. 'ud' --> start of a stance)
        out(queue(:,2),n) = 2;
        
    end
end
end