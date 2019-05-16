function [ trigger ] = findLocalExtremaKatsov(kind, x, g, w, medge)
%FINDLOCALEXTREMAKATSOV: A function to find local extrema in the vector x for
%each of the groups defined by the vector groups, using the peak-finding
%methods from Katsov et al 2017.

if length(x) ~= length(g)
    error('The length of the input vector must equal the length of the grouping vector');
end

% Find the groups
[G, ~] = findgroups(g);

% Find the extrema
[trigger] = splitapply(@(a){findextrema(a,w,medge,kind)}, x, G);

% Combine data for output
trigger = cat(1, trigger{:});

end

function [ y ] = findextrema(x, w, medge, kind)
switch kind
    case 'maxima'
        y = localmaxm(x,w,medge);
    case 'minima'
        y = localminm(x,w,medge);
    case 'all'
        y = localmaxm(x,w,medge) | localminm(x,w,medge);
end
end

function index = localminm(x, w, medge)

% INPUT: x is column vector; w is half-window length
%
% roger jang 1999 
% ak 2004 generalized
% Adapted from classify_submodes_cmsec_degsec_30fps.m by JZ-V, 2018
% Downloaded from https://datadryad.org/resource/doi:10.5061/dryad.854j2
% [m,n] = size(x);
% if ( m < n ) % row vector

% assume column vector
b1 = zeros(length(x), w);
b2 = zeros(length(x), w);
b3 = zeros(length(x), w);
b4 = zeros(length(x), w);

for i = 1:w
    b1(1+w:end-w,i) = x(1+w:end-w) < x(1+w-i:end-w-i);
	b2(1+w:end-w,i) = x(1+w:end-w) < x(1+w+i:end-w+i);
    b3(1+w:end-w,i) = x(1+w-i:end-w-i) < -medge;
	b4(1+w:end-w,i) = x(1+w+i:end-w+i) < -medge;
end
b5 = x < -medge;

index = ones(size(x));

for j = 1:w
    index = index & b1(:,j) & b2(:,j) & ( b3(:,j) | b4(:,j) );
end

index = index & b5;
index = logical(index);
end

function index = localmaxm(x, w, medge)

% INPUTS: x is column vector; w is half-window length
%
% roger jang 1999
% ak 2004 generalized
% Adapted from classify_submodes_cmsec_degsec_30fps.m by JZ-V, 2018
% Downloaded from https://datadryad.org/resource/doi:10.5061/dryad.854j2

% [m,n] = size(x);
% if ( m < n ) % row vector

% assume column vector
b1 = zeros(length(x), w);
b2 = zeros(length(x), w);
% continuity constraints
b3 = zeros(length(x), w);
b4 = zeros(length(x), w);

% >0 condition ensures search on pos values only;
% revise for valley search.
for i = 1:w
    b1(1+w:end-w,i) = x(1+w:end-w) > x(1+w-i:end-w-i);
	b2(1+w:end-w,i) = x(1+w:end-w) > x(1+w+i:end-w+i);
    b3(1+w:end-w,i) = x(1+w-i:end-w-i) > medge;
	b4(1+w:end-w,i) = x(1+w+i:end-w+i) > medge;
end
b5 = x > medge;

index = ones(size(x));

for j = 1:w
    index = index & b1(:,j) & b2(:,j) & ( b3(:,j) | b4(:,j) );
end

index = index & b5;
index = logical(index);

end
