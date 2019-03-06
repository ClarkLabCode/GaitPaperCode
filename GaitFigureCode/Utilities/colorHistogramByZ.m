function colorHistogramByZ( x, y, z, fn, xedges, yedges )
%This function plots a 2d histogram with bin colors given by the value
%taken by a scalar-valued function with handle fn in each bin.

if exist('sigma', 'var') && ~isempty(sigma) && sigma>0
    filtFlag = true;
else
    filtFlag = false;
end

% Get the histogram counts
if isempty(xedges) && isempty(yedges)
    [N,Xedges,Yedges,binX,binY] = histcounts2(x, y, 50, 'Normalization','Probability');
elseif isscalar(xedges) && iscalar(yedges) && (xedges == yedges)
    [N,Xedges,Yedges,binX,binY] = histcounts2(x, y, xedges, 'Normalization','Probability');
else
    [N,Xedges,Yedges,binX,binY] = histcounts2(x, y, xedges, yedges, 'Normalization','Probability');
end

% Remove null values
idx = (binX ~= 0) & (binY ~= 0);
z = z(idx);
binX = binX(idx);
binY = binY(idx);

% Get the values
[ C ] = accumarray([binX binY], z, [size(N,1), size(N,2)], fn);

% If desired, smoooth the histogram
if filtFlag
    C = imgaussfilt(C,sigma);
    
    imagesc([min(Xedges), max(Xedges)], [min(Yedges), max(Yedges)], C');
else
    % Show the histogram, making the areas where there is no data transparent
    imagesc([min(Xedges), max(Xedges)], [min(Yedges), max(Yedges)], C', 'AlphaData', (N ~= 0)');
end

axis('xy');

end

