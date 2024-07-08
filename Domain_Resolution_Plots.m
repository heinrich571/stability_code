%% Fresh start

close all
clear all
clc
startup

path_manager('add')


%% Load results

Results_Folder = '.\results\tests\';
Results_File   = 'Domain_Resolution_Test.mat';

load([Results_Folder Results_File]);


%% Calculate graph variables

N_eigenvalues = length(Solution(1).Eigenvalues);
Nx_ref = Problem(end).Domain.Nx + 1;
Ny_ref = Problem(end).Domain.Ny + 1;
eigenvalue_ref = zeros([N_eigenvalues 1]);
for j = 1:N_eigenvalues
    eigenvalue_ref(j) = Solution(end).Eigenvalues(j);
    mat_X_ref  = Solution(end).Domain.mat_X;
    mat_Y_ref  = Solution(end).Domain.mat_Y;
    u_ref(:,:,j) = reshape(Solution(end).Eigenfunctions.u(:,j), Ny_ref, Nx_ref);
    v_ref(:,:,j) = reshape(Solution(end).Eigenfunctions.v(:,j), Ny_ref, Nx_ref);
    w_ref(:,:,j) = reshape(Solution(end).Eigenfunctions.w(:,j), Ny_ref, Nx_ref);
    p_ref(:,:,j) = reshape(Solution(end).Eigenfunctions.p(:,j), Ny_ref, Nx_ref);
end

N_cases = length(Solution);
initializer = zeros([N_cases 1]);
log_N = initializer;
log_T = initializer;
Nx    = initializer;
t     = initializer;
initializer = zeros([N_cases N_eigenvalues]);
err_u = initializer;
err_v = initializer;
err_w = initializer;
err_p = initializer;
for i = 1:N_cases
    Nx(i) = Problem(i).Domain.Nx + 1;
    Ny(i) = Problem(i).Domain.Ny + 1;
    log_N(i) = log(Nx(i));
    if exist('Monitor', 'var')
        t(i)  = Monitor(i).Time;
        log_T(i) = log(t(i));
    end

    for j = 1:N_eigenvalues
        eigenvalue(i,j) = Solution(i).Eigenvalues(j);
        err_eigenvalue(i,j) = abs(eigenvalue(i,j) - eigenvalue_ref(j));
        mat_X = Solution(i).Domain.mat_X;
        mat_Y = Solution(i).Domain.mat_Y;
        u = reshape(Solution(i).Eigenfunctions.u(:,j), Ny(i), Nx(i));
        v = reshape(Solution(i).Eigenfunctions.v(:,j), Ny(i), Nx(i));
        w = reshape(Solution(i).Eigenfunctions.w(:,j), Ny(i), Nx(i));
        p = reshape(Solution(i).Eigenfunctions.p(:,j), Ny(i), Nx(i));
        err_u(i,j) = max(abs(interp2(mat_X, mat_Y, u, mat_X_ref, mat_Y_ref)-u_ref(:,:,j)), [], 'all');
        err_v(i,j) = max(abs(interp2(mat_X, mat_Y, v, mat_X_ref, mat_Y_ref)-v_ref(:,:,j)), [], 'all');
        err_w(i,j) = max(abs(interp2(mat_X, mat_Y, w, mat_X_ref, mat_Y_ref)-w_ref(:,:,j)), [], 'all');
        err_p(i,j) = max(abs(interp2(mat_X, mat_Y, p, mat_X_ref, mat_Y_ref)-p_ref(:,:,j)), [], 'all');
    end
end


%% Plots

% time vs. resolution
if exist('Monitor', 'var')
    figure('Name', 'Solution time vs. domain resolution', 'NumberTitle', 'off')
    plot(Nx, t, '-o')
    xlabel('$N$')
    ylabel('$t$')
    title('Solution time vs. domain resolution')
    grid minor
end


% Convergence of eigenvalues
eigenvalues_to_plot = [1 2 3 4 5];
for j = eigenvalues_to_plot
    figure('Name', ['Eigenvalue #' num2str(j) ' convergence'], 'NumberTitle', 'off')
    subplot(1,2,1)
    plot(Nx, real(eigenvalue(:,j)'))
    xlabel('$N$')
    ylabel(['$\mathrm{Re}\{\omega_' num2str(j) '\}$'])
    subplot(1,2,2)
    plot(Nx, imag(eigenvalue(:,j)'))
    xlabel('$N$')
    ylabel(['$\mathrm{Im}\{\omega_' num2str(j) '\}$'])

    figure('Name', ['Eigenvalue #' num2str(j) ' convergence'], 'NumberTitle', 'off')
    plot(Nx, err_eigenvalue(:,j)', '-^')
    xlabel('$N$')
    ylabel(['$\varepsilon_{\omega_' num2str(j) '}$'])
end


% Eigenvalue maps
for j = eigenvalues_to_plot
    figure('Name', ['Eigenvalue #' num2str(j) ' map'], 'NumberTitle', 'off')
    for i = 1:N_cases
        plot(real(eigenvalue(i,j)), imag(eigenvalue(i,j)), 'x')
    end
    xlabel(['$\mathrm{Re}\{\omega_ ' num2str(j) '\}$'])
    ylabel(['$\mathrm{Im}\{\omega_ ' num2str(j) '\}$'])
end