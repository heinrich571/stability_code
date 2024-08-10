%% Fresh Start

close all
clear all
clc


%% Common definitions

Problem.Computation.N_Workers = 1;

Problem.Domain.Nx       = [];
Problem.Domain.Ny       = [];
Problem.Domain.X_Limit  = [];
Problem.Domain.Y_Limit  = 200;
Problem.Domain.Y_Median = 2.4;

Problem.Physics.Beta                  = 2;
Problem.Physics.Number_Of_Eigenvalues = 20;

Problem.Base_Flow_Settings.initguess            = 1.232587656820289 + [-1 1]*1e-4;
Problem.Base_Flow_Settings.maxIterations        = 1e2;
Problem.Base_Flow_Settings.convergenceTolerance = 1e-6;

Sides_Boundary_Condition = 'zero_2nd_derivative_extrapolation';

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

Case_Name = 'X_Limit_Sensitivity_Test';


%% Generate Problems (huehue)

Nx_max = 80;
X_Limit_max = 150;

Average_Cell_Size = X_Limit_max/Nx_max;

X_Limit_vec = 10:10:X_Limit_max;
Nx_vec = round(X_Limit_vec/Average_Cell_Size);
Ny_vec = 80*ones(size(X_Limit_vec));
Nx_vec(logical(mod(Nx_vec, 2))) = Nx_vec(logical(mod(Nx_vec, 2)))+1;

Problem = repmat(Problem, [length(Nx_vec) 1]);

for i = 1:length(X_Limit_vec)
    Problem(i).Domain.X_Limit = X_Limit_vec(i);
    Problem(i).Domain.Nx = Nx_vec(i);
    Problem(i).Domain.Ny = Ny_vec(i);
end


%% Save

Save_Folder = '.\batch_runs\';

if ~isfolder(Save_Folder)
    mkdir(Save_Folder)
end
save([Save_Folder '\' Case_Name '.mat'], 'Case_Name', 'Problem')
