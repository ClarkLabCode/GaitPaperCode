function PlotAverageStepLengthAndDirection( newData, vfInterval )
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here

mmPerPix = 0.043;

%% Define variable lists

limbList = {'L1','L2','L3','R1','R2','R3'};

velVarList = {'angVel_radPerSec','forwardSpeed_mmPerSec','translationalSpeed_mmPerSec'};
velVarList = strcat('smooth_',velVarList);

% Define variable list for camera frame limb positions
limbVarListXCam = strcat(limbList, '_xCam');
limbVarListYCam = strcat(limbList, '_yCam');

% Define variable list for egocentric frame limb positions
limbVarListXEgo = strcat(limbList, '_xPlot');
limbVarListYEgo = strcat(limbList, '_yPlot');

% Define the list of swing/stance logicals
upDownVarList = strcat(limbList, '_down_cam');

%% Extract needed data from the table

% Extract velocity data
V = newData{:, velVarList};
V(:,1) = rad2deg(V(:,1));

% Extract camera frame limb positions
xCam = newData{:, limbVarListXCam};
yCam = newData{:, limbVarListYCam};

% Extract egocentric frame limb positions
xEgo = newData{:, limbVarListXEgo};
yEgo = newData{:, limbVarListYEgo};

ID = newData.videoID;

%% Locate the limb extreme positions

% Get the unique trajectories
G = findgroups(newData.uniqueFlyTrajID);

% Get the swing/stance logicals
D = newData{:, upDownVarList};

% Get the locations of all paired AEPs and PEPs in the dataset
tic;
A = splitapply(@(x) {local_get_indices(x)}, D, G);
A = cell2mat(A);
fprintf('Located limb extreme positions in %f seconds.\n', toc);

%% Compute the displacements

% Define the Euclidean distance metric
euc = @(p,a) sqrt(dot(p,p,2)+dot(a,a,2)-2*dot(a,p,2));

% Compute displacements and displacement angles in the camera frame
tic;
dispCam = nan(size(A));
for i = 1:6
    p = [xCam(A(:,i)==1,i), yCam(A(:,i)==1,i)];
    a = [xCam(A(:,i)==2,i), yCam(A(:,i)==2,i)];
    dispCam(A(:,i)==1,i) = euc(p,a);
end
fprintf('Computed camera frame displacements in %f seconds.\n', toc);

% Compute displacements in the egocentric frame
tic;
dispEgo = nan(size(A));
dispAngle = nan(size(A));
for i = 1:6
    p = [xEgo(A(:,i)==1,i),yEgo(A(:,i)==1,i)];
    a = [xEgo(A(:,i)==2,i),yEgo(A(:,i)==2,i)];
    dispEgo(A(:,i)==1,i) = euc(p,a);
    dispAngle(A(:,i)==1,i) = atan2(p(:,2)-a(:,2), p(:,1)-a(:,1));
end
fprintf('Computed egocentric frame displacements in %f seconds.\n', toc);

% Convert units
dispCam = mmPerPix * dispCam;
dispEgo = mmPerPix * dispEgo;

% Adjust displacement angle origin
dispAngle = dispAngle + pi/2;

%% Last-neighbor interpolation

tic;
% Allocate containers
dispCamInterp = nan(size(A));
dispEgoInterp = nan(size(A));
dispAngleInterp = nan(size(A));

% Initialize queues
qCam = nan(1,6);
qEgo = nan(1,6);
qAng = nan(1,6);

for i = 1:size(dispCam,1)
    for j = 1:6
        if isnan(dispCam(i,j)) && ~isnan(qCam(j))
            % Store the values
            dispCamInterp(i,j) = qCam(j);
            dispEgoInterp(i,j) = qEgo(j);
            dispAngleInterp(i,j) = qAng(j);
            
        elseif ~isnan(dispCam(i,j))
            % Update the queues
            qCam(j) = dispCam(i,j);
            qEgo(j) = dispEgo(i,j);
            qAng(j) = dispAngle(i,j);
            
            % Store the values
            dispCamInterp(i,j) = qCam(j);
            dispEgoInterp(i,j) = qEgo(j);
            dispAngleInterp(i,j) = qAng(j);
        end
    end
    
    % Cleare the queues
    if (i<size(dispCamInterp,1)) && (G(i+1)~=G(i))
        qEgo = nan(6,1);
        qCam = nan(6,1);
        qAng = nan(1,6);
    end
end
fprintf('Interpolated values in %f seconds\n', toc);

%% Symmetrize left and right turns

% Copy values
VSym = V;
dispAngleInterpSym = dispAngleInterp;
dispCamInterpSym = dispCamInterp;
dispEgoInterpSym = dispEgoInterp;

% Locate left turns
left = V(:,1) < 0;

% Fold left turns
VSym(left, [1,3]) = -VSym(left,[1,3]);
dispAngleInterpSym(left,:) = 2*pi-dispAngleInterpSym(left, [4,5,6,1,2,3]);
dispCamInterpSym(left,:) = dispCamInterpSym(left, [4,5,6,1,2,3]);
dispEgoInterpSym(left,:) = dispEgoInterpSym(left, [4,5,6,1,2,3]);


%% Plot averages with symmetrization
tic;
seriesLabels = {'O1','O2','O3','I1','I2','I3'};

% Step direction
PlotAverageModulationAsFunctionOfCentroidKinematics(VSym, dispAngleInterpSym,...
    ID, 'step direction (\circ)', seriesLabels, [-40 40], vfInterval, 'circular');

% Step length in camera frame
PlotAverageModulationAsFunctionOfCentroidKinematics(VSym, dispCamInterpSym,...
    ID, 'camera frame step length (mm)', seriesLabels, [0, 3], vfInterval);

fprintf('Computed symmetrized averages in %f seconds\n', toc);


end

%% Local functions

function [ out ] = local_get_indices(x)
out = zeros(size(x,1),6);
if size(x,1) > 15
    for n = 1:6
        queue = zeros(1,2);
        i = 1;
        op = false;
        for m = 1:size(x,1)-1
            if ~x(m,n) && x(m+1,n) && op
                queue(i,2) = m+1;
                i = i + 1;
                op = false;
            elseif x(m,n) && ~x(m+1,n) && ~op
                queue(i,1) = m+1;
                op = true;
            end
        end
        queue = queue(all(queue,2),:);
        out(queue(:,1),n) = 1;
        out(queue(:,2),n) = 2;
    end
end
end
