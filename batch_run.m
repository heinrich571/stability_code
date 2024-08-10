%% Fresh Start

close all
clear all
clc
path_manager('add')


%% Define Batch Runs

Folders = {'./batch_runs/'};

Filenames(1) = {'Y_Limit_Sensitivity_Test.mat'};
Filenames(2) = {'X_Limit_Sensitivity_Test.mat'};
Filenames(3) = {'Domain_Resolution_Sensitivity_Test.mat'};

Output_Folder = './batch_results/';

N_Workers = 80;


%% Run 'em

BatchX(Folders, Filenames, Output_Folder, N_Workers)


%% Cleanup

path_manager('remove')
