%% Fresh start

close all
clear all
clc
startup

path_manager('add')


%% Define base problem parameters

Basic_Problem.Computation.N_Workers = 1;

Basic_Problem.Domain.Nx       = -1;
Basic_Problem.Domain.Ny       = -1;
Basic_Problem.Domain.X_Limit  = 20;
Basic_Problem.Domain.Y_Limit  = 20;
Basic_Problem.Domain.Y_Median = 2.4;

Basic_Problem.Physics.Beta                  = 2;
Basic_Problem.Physics.Number_Of_Eigenvalues = 40;

Basic_Problem.Base_Flow_Settings.initguess            = [1.22 1.24];
Basic_Problem.Base_Flow_Settings.maxIterations        = 1e2;
Basic_Problem.Base_Flow_Settings.convergenceTolerance = 1e-6;

Basic_Problem.Boundary_Conditions.Wall.Pressure = 'LPPE';
Basic_Problem.Boundary_Conditions.Sides         = '2nd_derivative_extrapolation';

Basic_Problem.Flags.Display_Domain    = 0;
Basic_Problem.Flags.Display_Base_Flow = 0;
Basic_Problem.Flags.Generate_Report   = 0;

Case_ID        = 'Domain_Resolution_Test';
Results_Folder = '.\results\';


%% Define domain resolutions for convergence test

Nx_vec = [12:2:60];
Ny_vec = [12:2:60];


%% Solve problems

N = length(Nx_vec);
Problem = repmat(Basic_Problem, [N 1]);
Problem(1).Domain.Nx = Nx_vec(1);
Problem(1).Domain.Ny = Ny_vec(1);
Current_Case = ['DomainResTest_Nx_' num2str(Problem(1).Domain.Nx) '_Ny_' num2str(Problem(1).Domain.Ny)];
disp(['Now solving for: ' Current_Case])
disp('')

[Solution, Report] = BiGlobalTemporalSolver(Problem(1));

Solution = repmat(Solution, [N 1]);
Report   = repmat(Report, [N 1]);

parfor i = 2:N
    Problem(i).Domain.Nx = Nx_vec(i);
    Problem(i).Domain.Ny = Ny_vec(i);
    
    Current_Case = ['DomainResTest_Nx_' num2str(Problem(i).Domain.Nx) '_Ny_' num2str(Problem(i).Domain.Ny)];
    disp(['Now solving for: ' Current_Case])
    disp('')
    
    [Solution(i), Report(i)] = BiGlobalTemporalSolver(Problem(i));
end


%% Save results

dispstatus('SAVE RESULTS')
dispstatus('SAVE RESULTS', 0)

if ~isfolder(Results_Folder)
    mkdir(Results_Folder);
end
save([Results_Folder Case_ID '.mat'], 'Problem', 'Solution', 'Report');

dispstatus('SAVE RESULTS', 1)
dispstatus()


%% Cleanup

path_manager('remove')
