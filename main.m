%% Fresh start

% close all
clearvars
clc
path_manager('add')
startup


%% Define problem parameters

Problem.Computation.N_Workers = 1;

Problem.Domain.Nx       = 40;
Problem.Domain.Ny       = 40;
Problem.Domain.X_Limit  = 200;
Problem.Domain.Y_Limit  = 300;
Problem.Domain.Y_Median = 1*2.4;

Problem.Physics.Beta                  = 0.2;
Problem.Physics.Number_Of_Eigenvalues = 20;

Problem.Base_Flow_Settings.initguess            = 1.23258765682022 + [-1 1]*1e-5;
Problem.Base_Flow_Settings.maxIterations        = 1e2;
Problem.Base_Flow_Settings.convergenceTolerance = 1e-6;

% Sides_Boundary_Condition = 'zero_2nd_derivative';
Sides_Boundary_Condition = 'Linear_Extrapolation';

Problem.Boundary_Conditions.Top.u   = 'Dirichlet';
Problem.Boundary_Conditions.Top.v   = 'Dirichlet';
Problem.Boundary_Conditions.Top.w   = 'Dirichlet';
Problem.Boundary_Conditions.Top.p   = 'Dirichlet';

Problem.Boundary_Conditions.Right.u = Sides_Boundary_Condition;
Problem.Boundary_Conditions.Right.v = Sides_Boundary_Condition;
Problem.Boundary_Conditions.Right.w = Sides_Boundary_Condition;
Problem.Boundary_Conditions.Right.p = Sides_Boundary_Condition;

Problem.Boundary_Conditions.Left.u  = Sides_Boundary_Condition;
Problem.Boundary_Conditions.Left.v  = Sides_Boundary_Condition;
Problem.Boundary_Conditions.Left.w  = Sides_Boundary_Condition;
Problem.Boundary_Conditions.Left.p  = Sides_Boundary_Condition;

Problem.Boundary_Conditions.Wall.u  = 'Dirichlet';
Problem.Boundary_Conditions.Wall.v  = 'Dirichlet';
Problem.Boundary_Conditions.Wall.w  = 'Dirichlet';
Problem.Boundary_Conditions.Wall.p  = 'LPPE';


Problem.Flags.Display_Domain    = 0;
Problem.Flags.Display_Base_Flow = 0;
Problem.Flags.Generate_Report   = 0;

Case_ID        = 'test_lppe_linextrap';
Results_Folder = '.\results\';


%% Send to solver

[Domain, Base_Flow, Solution] = BiGlobalTemporalSolver(Problem);


%% Save results

dispstatus('SAVE RESULTS')
dispstatus('SAVE RESULTS', 0)

if ~isfolder(Results_Folder)
    mkdir(Results_Folder);
end
save([Results_Folder Case_ID '.mat'], 'Problem', 'Solution');

dispstatus('SAVE RESULTS', 1)
dispstatus()


%% Draw eigenfunctions

% Generate a report on the validity of the results against the
% Navier-Stokes equations and satisfying the original problem
Report = struct();
if Problem.Flags.Generate_Report
    Report.Navier_Stokes_Check = validateNS(Domain, Base_Flow, Problem, Solution);
    Report.EVP_Check = verifyEVP(Domain, Base_Flow, Problem, Solution);
end

Eigenvalue_Indices = 1;
Options.Solution_Index = 1;
Options.X_Limit = Problem.Domain.X_Limit;
Options.Y_Limit = Problem.Domain.Y_Limit;
view_results(Case_ID, Results_Folder, Eigenvalue_Indices, Options)
if ~isempty(fieldnames(Report))
    errors_plots(Solution, Report, Eigenvalue_Indices)
end
