%% Fresh start

close all
clearvars
clc
path_manager('add')


%% Define problem parameters

% Grid parameters

X_normalized_Limit  = 5;
Y_normalized_Limit  = 10;
Y_normalized_Median = 2.4;

Nx = 50;
Ny = Nx/2;


% Base flow parameters

Definitions.initguess            = [1.22 1.24];
Definitions.maxIterations        = 1e2;
Definitions.convergenceTolerance = 1e-6;


%% Get the computational domain

Domain = generate_domain(X_normalized_Limit, Y_normalized_Limit, Y_normalized_Median, Nx, Ny);


%% Get the base flow

Definitions.interval = Domain.vec_Y;
Base_Flow = get_base_flow(Definitions);


%% Create Matrices A & B for the BiGlobal eigenvalue problem

[mat_A, mat_B] = create_eigenvalue_matrices(Domain, Base_Flow, 0);
smat_A = sparse(mat_A);
smat_B = sparse(mat_B);


%% Solve the eigenvalue problem

Solution = solve_BiGlobal_eigenvalue_problem();
Tests = check_solution(Solution, mat_A, mat_B);

%% Plot the results

plot_results


%% Cleanup and remove added paths

path_manager('remove')


