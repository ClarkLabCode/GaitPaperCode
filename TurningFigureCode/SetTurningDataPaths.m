%% SETTURNINGDATAPATHS.m
% This script sets the paths to the .mat files containing the data and adds
% the required code directories to the Matlab path

% Top-level path; will also be used to store cache files for UMAP
basepath = '';

% Path to wild-type data file
wtdatapath = fullfile(basepath, '');

% Path to bristle activation data file
brdatapath = fullfile(basepath, '');

% Add code directories to Matlab path
addpath('Utilities');
addpath('AnalysisUtilities');