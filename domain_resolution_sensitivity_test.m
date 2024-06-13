%% Fresh start

close all
clear all
clc
startup

path_manager('add')


%% Define base problem parameters

Basic_Problem.Computation.N_Workers = 1;

Basic_Problem.Domain.Nx       = -999;
Basic_Problem.Domain.Ny       = -999;
Basic_Problem.Domain.X_Limit  = 10;
Basic_Problem.Domain.Y_Limit  = 20;
Basic_Problem.Domain.Y_Median = 2.4;

Basic_Problem.Physics.Beta                  = 2;
Basic_Problem.Physics.Number_Of_Eigenvalues = 40;

Basic_Problem.Base_Flow_Settings.initguess            = [1.22 1.24];
Basic_Problem.Base_Flow_Settings.maxIterations        = 1e2;
Basic_Problem.Base_Flow_Settings.convergenceTolerance = 1e-6;

Problem.Boundary_Conditions.Wall.Pressure = 'LPPE';

Basic_Problem.Flags.Display_Domain    = 1;
Basic_Problem.Flags.Display_Base_Flow = 0;

Case_ID        = 'Domain_Resolution_Test';
Results_Folder = '.\results\tests\';


%% Define domain resolutions for convergence test

Nx_vec = [20:2:60];
Ny_vec = [20:2:60];


%% Solve problems

for i = 1:length(Nx_vec)
    Problem(i) = Basic_Problem;
    Problem(i).Domain.Nx = Nx_vec(i);
    Problem(i).Domain.Ny = Ny_vec(i);
    
    Current_Case = ['DomainResTest_Nx_' num2str(Problem(i).Domain.Nx) '_Ny_' num2str(Problem(i).Domain.Ny)];
    disp(['Now solving for: ' Current_Case])
    disp('')
    
    Solution(i) = BiGlobalTemporalSolver(Problem(i));
end


%% Save results

dispstatus('SAVE RESULTS')
dispstatus('SAVE RESULTS', 0)

if ~isfolder(Results_Folder)
    mkdir(Results_Folder);
end
save([Results_Folder Case_ID '.mat'], 'Problem', 'Solution');

dispstatus('SAVE RESULTS', 1)
dispstatus()


%% Cleanup

path_manager('remove')
