function PlotModelExampleTrajectories( modelData, sample_length )
% This function plots the model example trajectories

% Select a set of trajectories based on forward speed
lut = unique(modelData{:,{'uniqueFlyTrajID','forwardSpeed_mmPerSec'}},'rows');
idx = [1, 30, 62, 72, 82, 99]; % Consider using this set of indices

ids = lut(idx,:);
num_ids = size(ids,1);
speeds = ids(:,2);

%% Plot just the left half of the trajectory
% Define a step order (limbList);
limbList = {'L3','L2','L1'};
limbUpDownList = strcat(limbList,'_down_cam');

makeFigure;
for n = 1:num_ids
    
    subplot(num_ids,1,n);
    
    % Extract the data
    fly = modelData{modelData.uniqueFlyTrajID == ids(n,1),limbUpDownList};
    step = fly(end-sample_length+1:end,:)';
    
    % Configure the plot
    % Set the tick marks on the y-axis
    set(gca, 'ytick', [1:3]);
    set(gca, 'yticklabels', limbList);
    
    % Invert the y-axis
    set(gca,'Ydir','reverse');
    
    % Set the x-axis
    tickLocs = [0:7.5:sample_length];
    set(gca, 'xtick', tickLocs);
    set(gca, 'xticklabels', [0:50:sample_length*1000/150]);

    % Visualize the image
    hold on;
    imagesc(step);
    colormap('gray');
    set(gca, 'Layer','top'); % Makes the tick marks visible over the plot
    
    % Provide an x-axis label
    xlabel('Time (ms)');
    
    % Provide an y-axis label
    ylabel(sprintf('%0.1f mm/s', speeds(n)));
    
    % Configure the axes
    axis('tight');
    ConfAxis;
    
end


%% Plot all limbs

% Define a step order (limbList);
limbList = {'R3','R2','R1','L3','L2','L1'};
limbUpDownList = strcat(limbList,'_down_cam');

makeFigure;
for n = 1:num_ids
    
    subplot(num_ids,1,n);
    
    % Extract the data
    fly = modelData{modelData.uniqueFlyTrajID == ids(n,1),limbUpDownList};
    step = fly(end-sample_length+1:end,:)';
    
    % Configure the plot
    % Set the tick marks on the y-axis
    set(gca, 'ytick', [1:6]);
    set(gca, 'yticklabels', limbList);
    
    % Invert the y-axis
    set(gca,'Ydir','reverse');
    
    % Set the x-axis
    tickLocs = [0:7.5:sample_length];
    set(gca, 'xtick', tickLocs);
    set(gca, 'xticklabels', [0:50:sample_length*1000/150]);

    % Visualize the image
    hold on;
    imagesc(step);
    colormap('gray');
    set(gca, 'Layer','top'); % Makes the tick marks visible over the plot
    
    % Provide an x-axis label
    xlabel('Time (ms)');
    
    % Provide an y-axis label
    ylabel(sprintf('%0.1f mm/s', speeds(n)));
    
    % Configure the axes
    axis('tight');
    ConfAxis;
    
end

end

