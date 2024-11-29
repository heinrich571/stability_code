%% Fresh Start

close all
clear all
clc

path_manager('add')


%% Common definitions

Problem.Computation.N_Workers = 1;

Problem.Domain.Nx       = [];
Problem.Domain.Ny       = [];
Problem.Domain.X_Limit  = [];
Problem.Domain.Y_Limit  = 100;
Problem.Domain.Y_Median = 2.4;

Problem.Physics.Beta                  = 0.25;
Problem.Physics.Number_Of_Eigenvalues = 20;

Problem.Base_Flow_Settings.initguess            = 1.232587656820289 + [-1 1]*1e-4;
Problem.Base_Flow_Settings.maxIterations        = 1e2;
Problem.Base_Flow_Settings.convergenceTolerance = 1e-6;

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

Case_Name = 'X_Limit_Sensitivity_Test';


%% Generate Problems (huehue)

Nx_min = 40;
Ny_min = 40;
X_Limit_max = 200;
xcd = 30;
ycd = Problem.Domain.Y_Limit;

rhon_req = Nx_min*Ny_min/(2*xcd*ycd);

X_Limit_vec = 50:10:X_Limit_max;
Problem = repmat(Problem, [length(X_Limit_vec) 1]);

Nx = Nx_min;
Ny = Nx;

for i = 1:length(X_Limit_vec)
    % Calculate required Nx & Ny
    Problem(i).Domain.X_Limit = X_Limit_vec(i);
    Problem(i).Domain.Nx = Nx;
    Problem(i).Domain.Ny = Ny;
    Domain = generate_domain(Problem(i));
    Nxcd = size(Domain.vec_X(Domain.vec_X >= -xcd & Domain.vec_X <= xcd), 1);
    Nycd = size(Domain.vec_X(Domain.vec_Y >= -ycd & Domain.vec_Y <= ycd), 1);
    rhon = Nxcd*Nycd/(2*xcd*ycd);
    while rhon < rhon_req
        Problem(i).Domain.Nx = Nx;
        Problem(i).Domain.Ny = Ny;
        Domain = generate_domain(Problem(i));
        Nxcd = size(Domain.vec_X(Domain.vec_X >= -xcd & Domain.vec_X <= xcd), 1);
        Nycd = size(Domain.vec_X(Domain.vec_Y >= -ycd & Domain.vec_Y <= ycd), 1);
        rhon = Nxcd*Nycd/(2*xcd*ycd);
        Nx = Nx + 6;
        Ny = Nx;
    end
    Problem(i).Domain.Nx = Nx;
    Problem(i).Domain.Ny = Problem(i).Domain.Nx;
    disp(i)
end

for i = 1:length(X_Limit_vec)
    Domain = generate_domain(Problem(i));
    show_domain(Domain)
    axis equal
end


%% Save

Save_Folder = '.\batch_runs\';

if ~isfolder(Save_Folder)
    mkdir(Save_Folder)
end
save([Save_Folder '\' Case_Name '.mat'], 'Case_Name', 'Problem')


%% Cleanup

path_manager('rm')
