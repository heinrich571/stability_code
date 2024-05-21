%% Fresh start

% close all
clearvars
clc
path_manager('add')
startup


%% Define problem parameters

% Parallel runs
% N_Workers = 1;
N_Workers = 6;
% N_Workers = 160;

p = gcp('nocreate');

if N_Workers ~= 1 && isempty(p)
    parpool(N_Workers);
end

% Grid parameters

X_normalized_Limit  = 40;
Y_normalized_Limit  = 25;
Y_normalized_Median = 2.4;

Nx = 50;
Ny = 50;

beta = 2;

% Case ID
% Case_ID = ['FFTest_Y' num2str(Y_normalized_Limit)];
% Case_ID = ['DomainARTest_X' num2str(X_normalized_Limit)];
% Case_ID = ['DomainResTest_Nx' num2str(Nx) '_Ny' num2str(Ny)];
Case_ID = ['DomainWallTest_' 'BCWPC'];


% Base flow parameters

Definitions.initguess            = [1.22 1.24];
Definitions.maxIterations        = 1e2;
Definitions.convergenceTolerance = 1e-6;


% Flags and other settings

Flags.Display_Domain    = 0;
Flags.Display_Base_Flow = 0;


%% Get the computational domain

dispstatus('DOMAIN GENERATION')
dispstatus('DOMAIN GENERATION', 0)

Domain = generate_domain(X_normalized_Limit, Y_normalized_Limit, Y_normalized_Median, Nx, Ny);
if Flags.Display_Domain
    show_domain(Domain);
end

dispstatus('Domain generation', 1)
dispstatus()


%% Get the base flow

dispstatus('BASE FLOW CALCULATION')
dispstatus('BASE FLOW CALCULATION', 0)

Definitions.interval = Domain.vec_Y;
Base_Flow = get_base_flow(Definitions);
if Flags.Display_Base_Flow
    show_baseflow(Domain, Base_Flow);
end

dispstatus('BASE FLOW CALCULATION', 1)
dispstatus()


%% Create Matrices A & B for the BiGlobal eigenvalue problem

dispstatus('GENERALIZED EIGENVALUE MATRICES FORMULATION')
dispstatus('GENERALIZED EIGENVALUE MATRICES FORMULATION', 0)

[mat_A, mat_B] = create_eigenvalue_matrices(Domain, Base_Flow, beta);

dispstatus('GENERALIZED EIGENVALUE MATRICES FORMULATION', 1)
dispstatus()


%% Solve the problem

dispstatus('EIGENVALUES CALCULATION')
dispstatus('EIGENVALUES CALCULATION', 0)

tic;
[Solution_Filtered, Solution_Raw] = biglobalsolver(mat_A, mat_B, Domain);
toc;

dispstatus('EIGENVALUES CALCULATION', 1)
dispstatus()


%% Save results

dispstatus('SAVE RESULTS')
dispstatus('SAVE RESULTS', 0)

restults_dir = './results/';
if ~isfolder(restults_dir)
    mkdir(restults_dir);
end
save([restults_dir Case_ID '.mat'], 'Case_ID', 'Definitions', 'Flags', 'mat_A', 'mat_B', 'Base_Flow', 'Domain', 'Solution_Filtered', 'Solution_Raw');

dispstatus('SAVE RESULTS', 1)
dispstatus()


%% Draw eigenfunctions

% omega_ind = 1;
% 
% ploteigfun(Domain, Solution_Filtered, 'u', omega_ind)
% ploteigfun(Domain, Solution_Filtered, 'v', omega_ind)
% ploteigfun(Domain, Solution_Filtered, 'w', omega_ind)
% ploteigfun(Domain, Solution_Filtered, 'p', omega_ind)

% ploteigenvalues(Solution_Filtered)

check_results
