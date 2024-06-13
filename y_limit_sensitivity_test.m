%% Fresh start

close all
clear all
clc
startup


%% Define base problem parameters

Basic_Problem.Computation.N_Workers = 1;

Basic_Problem.Domain.Nx       = 40;
Basic_Problem.Domain.Ny       = 40;
Basic_Problem.Domain.X_Limit  = 20;
Basic_Problem.Domain.Y_Limit  = 20;
Basic_Problem.Domain.Y_Median = 2.4;

Basic_Problem.Physics.Beta                  = 2;
Basic_Problem.Physics.Number_Of_Eigenvalues = 40;

Basic_Problem.Base_Flow_Settings.initguess            = [1.22 1.24];
Basic_Problem.Base_Flow_Settings.maxIterations        = 1e2;
Basic_Problem.Base_Flow_Settings.convergenceTolerance = 1e-6;

Basic_Problem.Flags.Display_Domain    = 0;
Basic_Problem.Flags.Display_Base_Flow = 0;

Case_ID        = 'Y_Limit_Sensitivity_Test';
Results_Folder = '.\results\tests\';


%% Define domain resolutions for convergence test

Y_Limit_vec = [10 15 20 25];


%% Solve problems

for i = 1:length(Y_Limit_vec)
    Problem(i) = Basic_Problem;
    Problem(i).Domain.Y_Limit = Y_Limit_vec(i);
    
    Current_Case = ['Y_Limit_Sensitivity_Y_Limit_' num2str(Problem(i).Domain.Y_Limit)];
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
