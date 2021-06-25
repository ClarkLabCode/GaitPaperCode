function [ Y ]  = ramerDouglasPeucker(X, epsilon)
% A simple MATLAB implementation of the (parametric) Ramer?Douglas?Peucker
% algorithm for simplifying a curve composed of line segments. Input X must
% be an array of two-dimensional row vectors, and epsilon a scalar
% tolerance parameter.

% See https://en.wikipedia.org/wiki/Ramer-Douglas-Peucker_algorithm for
% mathematical details. As this implementation is non-optimized, it is
% likely of approximate complexity O(n^2).

% Find the point with the maximum distance
dmax = 0;
index = 0;
for i = 2:(length(X)-1)
    d = perpendicularDistance(X(i,:), X(1,:), X(end,:));
    if ( d > dmax )
        index = i;
        dmax = d;
    end
end

% If max distance is greater than epsilon, recursively simplify
if ( dmax > epsilon )
    % Recursive call
    recResults1 = ramerDouglasPeucker(X(1:index,:), epsilon);
    recResults2 = ramerDouglasPeucker(X(index:end,:), epsilon);
    
    % Build the result list
    Y = [recResults1(1:length(recResults1)-1,:); recResults2(1:length(recResults2),:)];
else
    Y = [X(1,:); X(end,:)];
end

end

% Utility function to calculate the perpendicular distance from a point
% to a line.
function [ d ] = perpendicularDistance(P, Q1, Q2)

dX = Q2 - Q1;

if vecnorm(dX) == 0
    d = vecnorm(P-Q1);
else
    t = sum((P - Q1).*dX) / vecnorm(dX).^2;
    
    if t < 0
        d = vecnorm(P-Q1);
    elseif t > 1
        d = vecnorm(P-Q2);
    else
        d = vecnorm(P - (Q1 + t*dX));
    end
end
end