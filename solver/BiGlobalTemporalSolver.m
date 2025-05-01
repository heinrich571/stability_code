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

tic;
[eigenfunctions_matrix, eigenvalues_matrix, convergence_flag] = eigs(sparse(A),...
                                                                     sparse(B),...
                                                                     Problem.Physics.Number_Of_Eigenvalues,...
                                                                     'smallestabs', 'MaxIterations', 400, 'Display', true, ...
                                                                     'Tolerance', 1e-8);
% Insert u,v,w Dirichlet values at the top and botttom boundary, and
% Dirichlet p values at the top boundary
u_top_inds = get_var_top_inds('u', Nx, Ny);
v_top_inds = get_var_top_inds('v', Nx, Ny);
w_top_inds = get_var_top_inds('w', Nx, Ny);
p_top_inds = get_var_top_inds('p', Nx, Ny);
u_bottom_inds = get_var_bottom_inds('u', Nx, Ny);
v_bottom_inds = get_var_bottom_inds('v', Nx, Ny);
w_bottom_inds = get_var_bottom_inds('w', Nx, Ny);
inds = sort([u_top_inds v_top_inds w_top_inds p_top_inds u_bottom_inds v_bottom_inds w_bottom_inds]);
sz = size(eigenfunctions_matrix, 2);
for i = inds
    eigenfunctions_matrix = [eigenfunctions_matrix(1:i-1,:) ; zeros([1 sz]) ; eigenfunctions_matrix(i:end,:)];
end

toc;

dispstatus('EIGENVALUES CALCULATION', 1)
dispstatus()

% Organize raw output
Solution_Raw.Eigenvalues      = diag(eigenvalues_matrix);
Solution_Raw.Eigenfunctions.u = get_eigenfunction_of(eigenfunctions_matrix, 'u', Nx, Ny);
Solution_Raw.Eigenfunctions.v = get_eigenfunction_of(eigenfunctions_matrix, 'v', Nx, Ny);
Solution_Raw.Eigenfunctions.w = get_eigenfunction_of(eigenfunctions_matrix, 'w', Nx, Ny);
Solution_Raw.Eigenfunctions.p = get_eigenfunction_of(eigenfunctions_matrix, 'p', Nx, Ny);

% Normalize the solution for consistency, and build output variable
% nrm = Solution_Raw.Eigenfunctions.p(Ny,:);
nrm = zeros([1 length(Solution_Raw.Eigenvalues)]);
for i = 1:length(Solution_Raw.Eigenvalues)
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
function var_eigfun = get_eigenfunction_of(eigenfunctions_matrix, varname, Nx, Ny)

varnames = {'u' , 'v' , 'w' , 'p'};
varorder = find(strcmp(varnames, varname), 1, 'first');
var_inds = (1:(Nx*Ny)) + (varorder-1)*Nx*Ny;

var_eigfun = eigenfunctions_matrix(var_inds,:);

end
