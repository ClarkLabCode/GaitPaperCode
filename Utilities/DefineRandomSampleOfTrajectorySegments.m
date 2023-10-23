function [ trigger ] = DefineRandomSampleOfTrajectorySegments(newData, ns, winlen, samplingMethod, vfMin, vfMax, vfBin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    ns = 1e5;
end

if nargin < 3
    winlen = 15;
end

if nargin < 4
    samplingMethod = 'all';
end

if nargin < 5
    vfMin = 1;
end
if nargin < 6
    vfMax = 40;
end
if nargin < 7
    vfBin = 1;
end

% Only used if sampling method is 'trajectory_sequences'
if strcmp(samplingMethod, 'trajectory_sequences')
    samplesPerTrajectory = 100;
end

%% Define sample
tic;

% Extract the within-trajectory frame number to exclude timepoints that are
% within the half-window-length of either the start or the end of the
% trajectory
g = findgroups(newData.uniqueFlyTrajID);
idx1 = cell2mat(splitapply(@(x) {(1:length(x))'}, g, g));
idx2 = cell2mat(splitapply(@(x) {(length(x):-1:1)'}, g, g));
frameIdx = (idx1 > winlen) & (idx2 > winlen);
clearvars idx1 idx2;

% Extract the forward velocity of the fly
vf = newData.smooth_forwardSpeed_mmPerSec;
vr = newData.smooth_angVel_radPerSec;

% Allocate an indicator
trigger = zeros(height(newData),1,'logical');

% Select which sampling method to use
switch samplingMethod
    case 'all'
        % Sample uniformly from all points satisfying the desired forward
        % velocity bounds, without replacement
        idx = frameIdx & (vf > vfMin) & (vf <= vfMax);
        k = find(idx);
        p = randperm(length(k));
        trigger(k(p(1:ns))) = true;
        clearvars k p idx;
        
    case 'uniform_velocity'
        % Sample such that distributions of forward velocities at time zero
        % is approximately uniform
        
        % Discretize the forward velocity using the desired bins
        vfEdges = (vfMin:vfBin:vfMax)';
        nb = length(vfEdges)-1;
        vfDiscr = discretize(vf, vfEdges);
        
        % Define the number of points to sample per bin
        samplesPerBin = ceil(ns / nb);
        
        % Iterate over bins
        for ind = 1:nb
            % Sample uniformly from each bin, without replacement
            idx = frameIdx & (vfDiscr == ind);
            k = find(idx);
            p = randperm(length(k));
            trigger(k(p(1:samplesPerBin))) = true;
            clearvars k p idx;
        end
        
    case 'trajectory_sequences'
        % Sample some number of frames from each trajectory
        idx = frameIdx & (vf > vfMin) & (vf <= vfMax);
        
        trajList = unique(g(idx));
        n = accumarray(findgroups(g(idx)),1);
        
        % Remove trajectories with less than the desired number of possible
        % samples
        trajList = trajList(n > samplesPerTrajectory);
        
        % Randomly permute the remaining trajectories
        trajList = trajList(randperm(length(trajList)));
        
        ind = 0;
        n = 0;
        while (n < ns) && (ind < length(trajList))
            ind = ind + 1;
            
            % Take samplesPerTrajectory sequential samples at a random
            % offset from the start of the trajectory
            k = find(idx & (g==trajList(ind)));
            
            % Make sure that everything is sequential
            if all(diff(k)==1)
                % If so, do this the fast way
                p = randi([1 length(k)-samplesPerTrajectory],1);
                k = k(p:p+samplesPerTrajectory-1);
                trigger(k) = true;
                n = n + samplesPerTrajectory;
                %             else
                %                 % If not, fall into this slower case
                %                 k = k([(k(2)-k(1)); diff(k)] == 1);
                %                 if length(k) > samplesPerTrajectory
                %                     p = randi([1 length(k)-samplesPerTrajectory],1);
                %                     k = k(p:p+samplesPerTrajectory-1);
                %                     trigger(k) = true;
                %                     n = n + samplesPerTrajectory;
                %                 end
            end
        end
        
    case 'turns'
        
        % Sample uniformly from all points satisfying the desired forward
        % velocity bounds, without replacement
        idx = frameIdx & (vf > vfMin) & (vf <= vfMax);
        k = find(idx);
        p = randperm(length(k));
        trigger(k(p(1:ns))) = true;
        clearvars k p idx;
        
        % Ensure turns are included
        trigger = trigger | ((vf > vfMin) & (vf <= vfMax) & (newData.yawExtremum ~= 0));
    otherwise
        error('invalid sampling method: %s', samplingMethod);
end

clearvars frameIdx;
fprintf('Defined sample in %f seconds\n', toc);

end

