%% Fresh start

close all
clear all
clc

path_manager('add')

%% Load sensitivity test results

Results_Folder = '.\results\tests\';

Case_Names = {'Side_Extrapolation_Test_2nd_derivative'
              'Side_Extrapolation_Test_finite_difference'};


%% Draw solution for specified eigenvalues

eigenvalue_inds = [1 2 3 4 5 6];
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


