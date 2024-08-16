%% Fresh start

close all
clear all
clc
startup

path_manager('add')

%% Load sensitivity test results

% Results_Folder = '.\results\tests\';
Results_Folder = '.\results\';
Results_File   = 'Y_Limit_Sensitivity_Test.mat';

load([Results_Folder Results_File]);


%% Draw solution for specified eigenvalues

eigenvalue_inds = [1 2 3 4 5 6];
X_Stations = [0 5 10 15 20];

omega_ref = Solution(end).Eigenvalues(eigenvalue_inds);

N_Y_Limits = length(Solution);
Y_Limit_vec = zeros([1 N_Y_Limits]);

for i = N_Y_Limits:-1:1
    Y_Limit_vec(i) = max(Solution(i).Domain.vec_Y);
    if i == N_Y_Limits
        User_Figures = profileAtX(Solution(i), eigenvalue_inds, X_Stations, ['$H = ' num2str(Y_Limit_vec(i)) '$']);
    else
        User_Figures = profileAtX(Solution(i), eigenvalue_inds, X_Stations, ['$H = ' num2str(Y_Limit_vec(i)) '$'], User_Figures);
    end

    for j = eigenvalue_inds
        err_omega(j,i) = abs(Solution(i).Eigenvalues(j) - omega_ref(j));
    end
end

for j = eigenvalue_inds
    figure('Name', ['omega_' num2str(j) ' convergence'])
    plot(Y_Limit_vec, err_omega(j,:), '-o')
    title(['$\omega_' num2str(j) ' = ' num2str(omega_ref(j)) '$'])
    xlabel('$H$')
    ylabel(['$\varepsilon_{\omega_' num2str(j) '}$'])
end


