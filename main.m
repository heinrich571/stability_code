%% Fresh start

% close all
clearvars
clc
path_manager('add')
startup


%% Define problem parameters

Problem.Computation.N_Workers = 1;

Problem.Domain.Nx       = 50;
Problem.Domain.Ny       = 50;
Problem.Domain.X_Limit  = 20;
Problem.Domain.Y_Limit  = 20;
Problem.Domain.Y_Median = 2.4;

Problem.Physics.Beta                  = 2;
Problem.Physics.Number_Of_Eigenvalues = 40;

Problem.Base_Flow_Settings.initguess            = [1.22 1.24];
Problem.Base_Flow_Settings.maxIterations        = 1e2;
Problem.Base_Flow_Settings.convergenceTolerance = 1e-6;

Problem.Boundary_Conditions.Wall.Pressure = 'LPPE';
Problem.Boundary_Conditions.Sides         = '2nd_derivative_extrapolation';

Problem.Flags.Display_Domain    = 0;
Problem.Flags.Display_Base_Flow = 0;

Case_ID        = 'test_new_solver';
Results_Folder = '.\results\';


%% Send to solver

Solution = BiGlobalTemporalSolver(Problem);


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

check_results
