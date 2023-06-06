%% Fresh start

close all
clearvars
clc

Path_List = {'./grid_generation'
             './baseflow_definition'
             './'};
for i = 1:numel(Path_List)
    addpath(genpath(Path_List{i}));
end


%% Define problem parameters

% Grid parameters
X_Limit  = 5;
Y_Limit  = 10;
Y_Median = 0.3; 
Nx       = 100;
Ny       = Nx/2;

% Base flow parameters
Definitions.initguess            = [1.22 1.24];
Definitions.maxIterations        = 1e2;
Definitions.convergenceTolerance = 1e-6;
Definitions.interval             = [0 Y_Limit];
Definitions.a                    = 1;
Definitions.nu                   = 1.46e-5;


%% Create computational domain & derivative matrices

Domain = generate_domain(X_Limit, Y_Limit, Y_Median, Nx, Ny);


%% Get the base flow

Definitions.interval = Domain.vec_Y;
Base_Flow = get_base_flow(Domain, Definitions);


%% Create Matrices A & B for the BiGlobal eigenvalue problem

[mat_A, mat_B] = create_eigenvalue_matrices(Domain, Base_Flow);


%% Solve the eigenvalue problem

Solution = solve_BiGlobal_eigenvalue_problem();
% Add an output check function to the solution:
%   a) the function should test wether the obtained solution actually
%   solves the eigenvalue problem


%% Plot the results

plot_results


%% Cleanup and remove added paths

for i = 1:numel(Path_List)
    rmpath(genpath(Path_List{i}));
end


