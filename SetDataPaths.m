%% SETDATAPATHS.m
% This script sets the paths to the .mat files containing the data and adds
% the required code directories to the Matlab path

% Top-level path; will also be used to store cache files for UMAP
basepath = '';

% Path to wild-type data file
wtdatapath = fullfile(basepath, '20181025_20180530-20180614_IsoD1_Glass_MaskedModel_1000PCs_amplitude_phase_down_downcam_Steps(_down_cam).mat');

% Add code directories to Matlab path
addpath('AnalysisUtilities','Utilities');