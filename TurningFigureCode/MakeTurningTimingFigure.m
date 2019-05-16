%% MAKETURNINGTIMINGFIGURE.m
% This script generates the turning timing figures from 
% "The manifold structure of limb coordination in walking Drosophila"

% Set data paths
SetTurningDataPaths;

% Set minimum forward speed for yaw extremum detection
vfMin = 0.5;

%% Load data from file and prepare data

load(wtdatapath, 'newData');
[ newData ] = PrepareTurningData(newData, vfMin);

%% Set preferences

% Set number of permutations for Monte Carlo
nperm = 1e5;

% Set edges for discretization
phiBinEdges = (0:1/180:1)';

% Select whether to plot bootstrapped confidence intervals on PDFs
plotCI = true;

%% Limb position distributions at yaw extrema

% At all points
PlotLimbPositionDistributionsAtAllPoints(newData);

% At yaw extrema
PlotLimbPositionDistributionsAtTimeOfTurn(newData);

%% Phase distributions at yaw extrema

% Compute for all points
[ fPDFnull, fCDFnull, xq, nnull, xnull, idnull ]  = PlotPhaseDistributionsAtAllPoints(newData, phiBinEdges, plotCI);

% Compute for turns
[ fPDFturn, fCDFturn, xq, nturn, xturn, idturn ] = PlotPhaseDistributionsAtTimeOfTurn(newData, phiBinEdges, plotCI);

%% Perform Kuiper tests

p = nan(6,1);
pci = nan(6,2);
vstat = nan(6,1);
for ind = 1:6
    
    countsNull = accumarray([discretize(xnull(:,ind), phiBinEdges), idnull], 1);
    countsTurn = accumarray([discretize(xturn(:,ind), phiBinEdges), idturn], 1);
    
    [p(ind), pci(ind,:), vstat(ind), ~] = kuiper2permtest(countsNull,countsTurn, nperm);
end


