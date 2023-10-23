function [ cmpBlueRed, cmpRed, cmpBlue ] = MakeTurningColormaps()

% Define the full colormap
load('blueRedColorMap.mat', 'cmpBlueRed');

if nargout > 1
    [~,ind] = max(vecnorm(cmpBlueRed,2,2));
    cmpRed = cmpBlueRed(ind:end,:);
    cmpBlue = cmpBlueRed(1:ind,:);
end

end
