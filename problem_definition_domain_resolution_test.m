%% Fresh Start

close all
clear all
clc


%% Common definitions

Problem.Computation.N_Workers = 1;

Problem.Domain.Nx       = [];
Problem.Domain.Ny       = [];
Problem.Domain.X_Limit  = 200;
Problem.Domain.Y_Limit  = 300;
Problem.Domain.Y_Median = 2.4;

Problem.Physics.Beta                  = 0.25;
Problem.Physics.Number_Of_Eigenvalues = 20;

Problem.Base_Flow_Settings.initguess            = 1.23258765682022 + [-1 1]*1e-5;
Problem.Base_Flow_Settings.maxIterations        = 1e2;
Problem.Base_Flow_Settings.convergenceTolerance = 1e-6;

Sides_Boundary_Condition = 'zero_2nd_derivative';

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

Case_Name = 'Domain_Resolution_Sensitivity_Test';


%% Generate Problems (huehue)

Nx_vec = 20:10:40;
Ny_vec = 20:20:200;

Problem = repmat(Problem, [length(Nx_vec)*length(Ny_vec) 1]);

k = 1;
for i = 1:length(Nx_vec)
    Nx = Nx_vec(i);
    for j = 1:length(Ny_vec)
        Ny = Ny_vec(j);
        Problem(k).Domain.Nx = Nx;
        Problem(k).Domain.Ny = Ny;
        k = k + 1;
    end
end


%% Save

Save_Folder = '.\batch_runs\';

if ~isfolder(Save_Folder)
    mkdir(Save_Folder)
end
save([Save_Folder '\' Case_Name '.mat'], 'Case_Name', 'Problem')
