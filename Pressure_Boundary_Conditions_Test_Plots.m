%% Fresh start

close all
clear all
clc

path_manager('add')

%% Load sensitivity test results

Results_Folder = '.\results\tests\';

Case_Names = {'Pressure_Boundary_Condition_Test_WPC_TVANISH_SEXTRAP'
              'Pressure_Boundary_Condition_Test_WLPPE_TVANISH_SEXTRAP'
              'Pressure_Boundary_Condition_Test_WLPPE_TLPPE_SEXTRAP'
              'Pressure_Boundary_Condition_Test_WLPPE_TVANISH_SLPPE'
              'Pressure_Boundary_Condition_Test_WLPPE_TLPPE_SLPPE'};


%% Draw solution for specified eigenvalues

eigenvalue_inds = [3];
X_Stations = [0 5 10 15 20];

N_Cases = length(Case_Names);

for i = 1:N_Cases
    load([Results_Folder Case_Names{i} '.mat']);
    if i == 1
        User_Figures = profileAtX(Solution, eigenvalue_inds, X_Stations, Case_Names{i});
    else
        User_Figures = profileAtX(Solution, eigenvalue_inds, X_Stations, Case_Names{i}, User_Figures);
    end
    clearvars Problem Solution Report
end


