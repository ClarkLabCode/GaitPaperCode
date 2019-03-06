function [ f, gof, output ] = calcPowerLaw( newData )
% This function calculates the line of best fit for the powerlaw governing
% stance duration.

%% Get the stance duration and forward velocity data for each of the limbs.

% Define a limb list
limbList = {'L1','L2','L3','R1','R2','R3'};

% Define some variables for gathering the data
velVarList = {'_StanceAvgVf_mmPerSec'};

% Define some logical variables for isolating events
stanceVarList = strcat(limbList, '_StanceStart');
logicals = newData{:,stanceVarList};

data = [];

for i =1:6
    
    % Get the velocity data and associated variable
    steps = newData{logicals(:,i), [strcat(limbList(i), velVarList),strcat(limbList(i),'_StanceDur_mmPerSec')]};
    
    % Filter out NaN values
    steps(isnan(steps(:,2)),:) = [];
    
    % Add the current limb's step to the full dataset
    data = [data; steps];
     
end

% Filter out all the entries that have a NaN value as the forward speed
data(isnan(data(:,1)),:) = [];

% Filter out data where the forward speed is outside the range 5-35 mm/s
data(data(:,1) < 5,:) = [];
data(data(:,1) > 35,:) = [];

vel = data(:,1);
stance = data(:,2);

%% Fit the powerlaw

vel = data(:,1);
stance = data(:,2);

[f, gof, output] = fit(vel,stance,'power1');


end

