%% Define variables to extract

t = (-winlen:winlen)';
t_ms = t / 0.15;

limbList = {'L1','L2','L3','R1','R2','R3'};

velVarList = {'angVel_radPerSec','forwardSpeed_mmPerSec','translationalSpeed_mmPerSec'};

downVarList = strcat(limbList, '_down_cam');

limbVarListX = strcat(limbList, '_xPlot_mm');
limbVarListY = strcat(limbList, '_yPlot_mm');

phaseVarList = strcat('InstantaneousPhase_', limbList, 'y');
frequencyVarList = strcat('InstantaneousFrequency_', limbList, 'y');

centroidVarList = {'xCOM','yCOM'};

if smoothing
    velVarList = strcat('smooth_', velVarList);
    phaseVarList = strcat('smooth_', phaseVarList);
    centroidVarList = strcat('smooth_', centroidVarList);
end

idVarList = {'uniqueFlyTrajID','Frame','yawExtremum', 'videoID'};

varList = [idVarList, velVarList, limbVarListX, limbVarListY, phaseVarList, downVarList, frequencyVarList, centroidVarList];

%% Define sample

[ trigger ] = DefineRandomSampleOfTrajectorySegments(newData, ns, winlen, samplingMethod, vfMin, vfMax, vfBin);

%% Extract data

% Extract slices
tic;
[ S ] = getSlicesSimple(newData, varList, winlen, trigger);
S = cat(3, S{:});

% Extract out each desired array
ID = S(:, 1:4, :);
V = S(:, 5:7, :);
L = S(:, 8:19, :);
Phi = S(:, 20:25, :);
D = S(:, 26:31, :);
dPhi = S(:, 32:37,:);
CM = S(:, 38:39,:);
clearvars S;

% Convert units
V(:,1,:) = rad2deg(V(:,1,:));
Phi = mod(Phi, 2*pi) / (2*pi);
dPhi = dPhi / (2*pi) * 150;

% Flip the axes such that +x is rightward and +y is forward
L = - L;

fprintf('Extracted snippets in %f seconds\n', toc);

% Identify turns
isExtremum = squeeze(ID(t==0,3,:));
isTurn = isExtremum ~= 0;

% If turns are forcibly included, fold left and right turns together
if strcmp(samplingMethod, 'turns') && (foldTurns)
    
    % Identify left turns
    left = (isExtremum < 0);
    
    % Fold centroid velocities
    V(:, [1,3], left) = - V(:, [1,3], left);
    
    % Fold limb positions
    L(:, :, left) = L(:,[4,5,6,1,2,3,10,11,12,7,8,9], left);
    L(:, 1:6, left) = - L(:, 1:6, left);
    
    % Fold phases, frequencies, and swing/stance
    Phi(:,:,left) = Phi(:,[4,5,6,1,2,3], left);
    dPhi(:,:,left) = dPhi(:,[4,5,6,1,2,3], left);
    D(:,:,left) = D(:,[4,5,6,1,2,3], left);
    
end

% Mean-subtract the limb positions
Lcentered = L - mean(L,1);

%% Compute the curvature

[ curvature ] = fast_2d_curvature(3,t,squeeze(CM(:,1,:))', squeeze(CM(:,2,:))');

%% Preprocess data

tic;
if meanSubtractTrajectories
    X_V = reshape(permute(V, [3 1 2]), [], size(V,1)*size(V,2));
    X_V = (X_V - nanmean(X_V,1)) ./ nanstd(X_V,0,1);

    X_L = reshape(permute(Lcentered, [3 1 2]), [], size(L,1)*size(L,2));
    X_L = (X_L - nanmean(X_L,1))./nanstd(X_L,0,1);
else
    X_L = reshape(permute(L, [3 1 2]), [], size(L,1)*size(L,2));
    X_V = reshape(permute(V, [3 1 2]), [], size(V,1)*size(V,2));
end


fprintf('Preprocessed centroid and limb data in %f seconds.\n', toc);

%% Separate out turning data

if strcmp(samplingMethod, 'turns')
    X_Lturn = X_L(isTurn,:);
    X_Vturn = X_V(isTurn,:);
    Lturn = L(:,:,isTurn);
    Phiturn = Phi(:,:,isTurn);
    
    Vturn = V(:,:,isTurn);
    IDturn = ID(:,:,isTurn);
    LcenteredTurn = Lcentered(:,:,isTurn);
end