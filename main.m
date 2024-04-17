%% Fresh start

close all
clearvars
clc
path_manager('add')
startup


%% Define problem parameters

% Grid parameters

X_normalized_Limit  = 20;
Y_normalized_Limit  = 20;
Y_normalized_Median = 2.4;

Nx = 40;
Ny = Nx/2;

beta = 0.255;


% Base flow parameters

Definitions.initguess            = [1.22 1.24];
Definitions.maxIterations        = 1e2;
Definitions.convergenceTolerance = 1e-6;


% Flags and other settings

Flags.Display_Domain    = 1;
Flags.Display_Base_Flow = 1;


%% Get the computational domain

Domain = generate_domain(X_normalized_Limit, Y_normalized_Limit, Y_normalized_Median, Nx, Ny);
if Flags.Display_Domain
    show_domain(Domain);
end


%% Get the base flow

Definitions.interval = Domain.vec_Y;
Base_Flow = get_base_flow(Definitions);
if Flags.Display_Base_Flow
    show_baseflow(Domain, Base_Flow);
end


%% Create Matrices A & B for the BiGlobal eigenvalue problem

[mat_A, mat_B] = create_eigenvalue_matrices(Domain, Base_Flow, beta);
[eigenfunctions_matrix,eigenvalues_matrix] = eig(mat_A, mat_B);
% evalue = diag(eigenvalues_matrix);
% smat_A = sparse(mat_A);
% smat_B = sparse(mat_B);
% 
% 
% %% Solve the eigenvalue problem
% 
% Solution = solve_BiGlobal_eigenvalue_problem();
% Tests = check_solution(Solution, mat_A, mat_B);
% 
% %% Plot the results
% 
% plot_results
% 
% 
% %% Cleanup and remove added paths
% 
% path_manager('remove')

save('output_test.mat', 'eigenfunctions_matrix', 'eigenvalues_matrix')


