%% Fresh start

close all
clear all
clc

path_manager('add')

%% Load sensitivity test results

Results_Folder = '.\results\tests\';
Results_File   = 'X_Limit_Sensitivity_Test.mat';

load([Results_Folder Results_File]);


%% Draw solution for specified eigenvalues

eigenvalue_inds = [1 2 3 4 5 6];
X_Stations = [0 2 5 10];

omega_ref = Solution(end).Eigenvalues(eigenvalue_inds);

N_X_Limits = length(Solution);
X_Limit_vec = zeros([1 N_X_Limits]);
err_omega = zeros([1 N_X_Limits]);

for i = 1:N_X_Limits
    X_Limit_vec(i) = max(Solution(i).Domain.vec_X);
    if i == 1
        User_Figures = profileAtX(Solution(i), eigenvalue_inds, X_Stations, ['$L = ' num2str(2*X_Limit_vec(i)) '$']);
    else
        User_Figures = profileAtX(Solution(i), eigenvalue_inds, X_Stations, ['$L = ' num2str(2*X_Limit_vec(i)) '$'], User_Figures);
    end
    for j = eigenvalue_inds
        err_omega(j,i) = abs(Solution(i).Eigenvalues(j) - omega_ref(j));
    end
end

for j = eigenvalue_inds
    figure('Name', ['omega_' num2str(j) ' convergence'])
    plot(2*X_Limit_vec, err_omega(j,:))
    title(['$\omega_' num2str(j) ' = ' num2str(omega_ref(j)) '$'])
    xlabel('$L$')
    ylabel(['$\varepsilon_{\omega_' num2str(j) '}$'])
end


