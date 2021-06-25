function [ dX ] = savitzkyGolayFilter( newData, varList, dorder, forder, wlen )
%SAVITZKYGOLAYFILTER: A function to smooth data and estimate derivatives
%using Savitzky-Golay filtering on a per-trajectory basis.

%% Form the filter

% NOTE: dt is assumed to be unity
dt = 1;

% Compute the projection matrix
[ b ] = localsg( wlen, forder, dorder );

% Scale the filter appropriately
if dorder > 0
    b = b/(-dt)^dorder;
end

%% Apply the filter to the data

% Extract numerical data from table
data = newData{:, varList};

% Find the groups corresponding to individual trajectories
[G, ~] = findgroups(newData.uniqueFlyTrajID);

% Estimate the derivatives for the timeseries of each trajectory
[ dX ] = splitapply(@(x){localsgfilt( x, b, wlen )}, data, G);

% Combine data for output
dX = cat(1, dX{:});

end

function [ b ] = localsg( wlen, forder, dorder )
% Optimized Savitzky-Golay filter projection matrix computation

% Validate inputs
if (round(wlen)~=wlen) || ~mod(wlen,2)
    error('Window length must be an odd integer.');
end
if forder > wlen - 1
    error('Filter order must be less than the window length minus one.');
end
if dorder > forder
    error('Derivative order must be less than or equal to the filter order.');
end

% Compute the half-window
halfwin = fix((wlen-1)/2);
x = (-halfwin:halfwin)';

% Form the Vandermonde matrix
v = x .^ (0:forder);

% Optimize for dorder = 0
if dorder == 0
    b = (v * (v \ eye(wlen)))';
else
    k = (1:(forder-dorder));
    hx = [(zeros(wlen,dorder)), ones(wlen,1)*prod(1:dorder), x.^k .* (k+1)];
    b = (hx * (v \ eye(wlen)))';
end
end

function [ y ] = localsgfilt( x, b, wlen )
% Apply the Savitzky-Golay filter defined by the projection matrix b to
% data stored as columns of the matrix x. This implementation assumes that
% missing (NaN) values occur only at the beginning and end of timeseries

% Allocate output
y = nan(size(x));
if size(x,1) > wlen
    % Find NaN rows
    xnotnan = ~any(isnan(x),2);
    if nnz(xnotnan) > wlen
        % Remove NaN rows
        x = x(xnotnan,:);
        
        % Compute the transient on
        ybegin = fliplr(b(:,(wlen-1)/2+2:end))' * flipud(x(1:wlen,:));
        
        % Compute the steady state output
        ycenter = filter(b(:,(wlen-1)./2+1), 1, x);
        
        % Compute the transient off
        yend = fliplr(b(:,1:(wlen-1)/2))' * flipud(x(end-(wlen-1):end,:));
        
        % Store output
        y(xnotnan,:) = [ybegin; ycenter(wlen:end,:); yend];
    end
end
end
