function [Domain, Base_Flow, Solution] = BiGlobalTemporalSolver(Problem)

% Expand problem parameters for convenient code
Nx = Problem.Domain.Nx + 1;
Ny = Problem.Domain.Ny + 1;
N_Workers = Problem.Computation.N_Workers;

% Get computation parameters
p = gcp('nocreate');
if N_Workers ~= 1 && isempty(p)
    parpool(N_Workers);
end

% Generate computational domain
dispstatus('DOMAIN GENERATION')
dispstatus('DOMAIN GENERATION', 0)

Domain = generate_domain(Problem);
if Problem.Flags.Display_Domain
    show_domain(Domain);
end

dispstatus('Domain generation', 1)
dispstatus()

% Get the base flow
dispstatus('BASE FLOW CALCULATION')
dispstatus('BASE FLOW CALCULATION', 0)

Base_Flow_Settings = Problem.Base_Flow_Settings;
Base_Flow_Settings.interval = Domain.vec_Y;
Base_Flow = get_base_flow(Base_Flow_Settings);

if Problem.Flags.Display_Base_Flow
    show_baseflow(Domain, Base_Flow);
end

dispstatus('BASE FLOW CALCULATION', 1)
dispstatus()

% Create eigenvalue problem matrices

dispstatus('GENERALIZED EIGENVALUE MATRICES FORMULATION')
dispstatus('GENERALIZED EIGENVALUE MATRICES FORMULATION', 0)

[A, B] = evpmatrices(Domain, Base_Flow, Problem);

dispstatus('GENERALIZED EIGENVALUE MATRICES FORMULATION', 1)
dispstatus()

% Find the eigenvalues and eigenfunctions
dispstatus('EIGENVALUES CALCULATION')
dispstatus('EIGENVALUES CALCULATION', 0)

if Problem.Flags.Display_Operators
    plot_operators(A, B, Domain);
end

tic;
[efmat, evmat, convergence_flag] = eigs(A, B, ...
                                        Problem.Physics.Number_Of_Eigenvalues,...
                                        'smallestabs', 'MaxIterations', 400, 'Display', true, ...
                                        'Tolerance', 1e-8);
% Insert u,v,w Dirichlet values at the top and botttom boundary, and
% Dirichlet p values at the top boundary
efmat = place_trivial_values_at_boundaries(efmat, Nx, Ny);

toc;

dispstatus('EIGENVALUES CALCULATION', 1)
dispstatus()

% Organize raw output
Solution_Raw.Eigenvalues      = diag(evmat);
Solution_Raw.Eigenfunctions.u = get_eigenfunction_of('u', efmat, Nx, Ny);
Solution_Raw.Eigenfunctions.v = get_eigenfunction_of('v', efmat, Nx, Ny);
Solution_Raw.Eigenfunctions.w = get_eigenfunction_of('w', efmat, Nx, Ny);
Solution_Raw.Eigenfunctions.p = get_eigenfunction_of('p', efmat, Nx, Ny);

% Normalize the solution for consistency, and build output variable
% nrm = Solution_Raw.Eigenfunctions.p(Ny,:);
nrm = zeros([1 length(Solution_Raw.Eigenvalues)]);
for i = 1 : length(Solution_Raw.Eigenvalues)
    abs_p = abs(Solution_Raw.Eigenfunctions.p(:,i));
    ind_max = find(abs_p == max(abs_p), 1, 'first');
    nrm(i) = Solution_Raw.Eigenfunctions.p(ind_max,i);
end

Solution.Domain           = Domain;
Solution.Physics          = Problem.Physics;
Solution.Eigenvalues      = Solution_Raw.Eigenvalues;
Solution.Eigenfunctions.u = Solution_Raw.Eigenfunctions.u./nrm;
Solution.Eigenfunctions.v = Solution_Raw.Eigenfunctions.v./nrm;
Solution.Eigenfunctions.w = Solution_Raw.Eigenfunctions.w./nrm;
Solution.Eigenfunctions.p = Solution_Raw.Eigenfunctions.p./nrm;
Solution.Normalizers      = nrm;
Solution.Convergence_Flag = convergence_flag;

end


% Supporting functions
function var_eigfun = get_eigenfunction_of(Variable_Name, Eigenfunctions_Matrix, Nx, Ny)

varnames = {'u' , 'v' , 'w' , 'p'};
varorder = find(strcmp(varnames, Variable_Name), 1, 'first');
var_inds = (1:(Nx*Ny)) + (varorder-1)*Nx*Ny;

var_eigfun = Eigenfunctions_Matrix(var_inds,:);

end

function efmat = place_trivial_values_at_boundaries(efmat, Nx, Ny)

sz = size(efmat, 2);

u_top_inds    = get_variable_top_indices   ('u', Nx, Ny);
v_top_inds    = get_variable_top_indices   ('v', Nx, Ny);
w_top_inds    = get_variable_top_indices   ('w', Nx, Ny);
p_top_inds    = get_variable_top_indices   ('p', Nx, Ny);
u_bottom_inds = get_variable_bottom_indices('u', Nx, Ny);
v_bottom_inds = get_variable_bottom_indices('v', Nx, Ny);
w_bottom_inds = get_variable_bottom_indices('w', Nx, Ny);

inds = sort([u_top_inds v_top_inds w_top_inds p_top_inds u_bottom_inds v_bottom_inds w_bottom_inds]);

for i = inds
    efmat = [efmat(1:i-1,:) ; zeros([1 sz]) ; efmat(i:end,:)];
end

end
