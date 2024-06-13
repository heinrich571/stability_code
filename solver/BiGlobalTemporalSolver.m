function Solution = BiGlobalTemporalSolver(Problem)

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

[mat_A, mat_B] = create_eigenvalue_matrices(Domain, Base_Flow, Problem);

dispstatus('GENERALIZED EIGENVALUE MATRICES FORMULATION', 1)
dispstatus()

% Find the eigenvalues and eigenfunctions
dispstatus('EIGENVALUES CALCULATION')
dispstatus('EIGENVALUES CALCULATION', 0)

tic;
[eigenfunctions_matrix, eigenvalues_matrix] = eigs(sparse(mat_A), sparse(mat_B), Problem.Physics.Number_Of_Eigenvalues, 'SM');
toc;

dispstatus('EIGENVALUES CALCULATION', 1)
dispstatus()

% Organize raw output
Solution_Raw.Eigenvalues      = diag(eigenvalues_matrix);
Solution_Raw.Eigenfunctions.u = get_eigenfunction_of(eigenfunctions_matrix, 'u', Nx, Ny);
Solution_Raw.Eigenfunctions.v = get_eigenfunction_of(eigenfunctions_matrix, 'v', Nx, Ny);
Solution_Raw.Eigenfunctions.w = get_eigenfunction_of(eigenfunctions_matrix, 'w', Nx, Ny);
Solution_Raw.Eigenfunctions.p = get_eigenfunction_of(eigenfunctions_matrix, 'p', Nx, Ny);

% Get eigenvalues within the specified limitations
eigenvalue_max_magnitude = 100;
inds = find(sqrt(real(Solution_Raw.Eigenvalues).^2 + imag(Solution_Raw.Eigenvalues).^2) <= eigenvalue_max_magnitude);

% Normalize the solution for consistency, and build output variable
nrm = max(Solution_Raw.Eigenfunctions.p(:,inds), [], 1);
% nrm = 1;

Solution.Domain           = Domain;
Solution.Eigenvalues      = Solution_Raw.Eigenvalues(inds);
Solution.Eigenfunctions.u = Solution_Raw.Eigenfunctions.u(:,inds)./nrm;
Solution.Eigenfunctions.v = Solution_Raw.Eigenfunctions.v(:,inds)./nrm;
Solution.Eigenfunctions.w = Solution_Raw.Eigenfunctions.w(:,inds)./nrm;
Solution.Eigenfunctions.p = Solution_Raw.Eigenfunctions.p(:,inds)./nrm;

% Generate a report on the validity of the results against the Navier-Stokes equations
Report = validateNS(Domain, Base_Flow, Problem, Solution);
Solution.Report = Report;

end


% Supporting functions
function var_eigfun = get_eigenfunction_of(eigenfunctions_matrix, varname, Nx, Ny)

varnames = {'u' , 'v' , 'w' , 'p'};
varorder = find(strcmp(varnames, varname), 1, 'first');
var_inds = (1:(Nx*Ny)) + (varorder-1)*Nx*Ny;

var_eigfun = eigenfunctions_matrix(var_inds,:);

end
