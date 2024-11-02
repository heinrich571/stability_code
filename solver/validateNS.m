function Navier_Stokes_Check = validateNS(Domain, Base_Flow, Problem, Solution)

%% Build eigenvalue problem A and B matrices

[A, B] = bglnsematrices(Domain, Base_Flow, Problem);


%% Check the validity of the solution

Ny  = size(Domain.mat_X, 1);
Nx  = size(Domain.mat_X, 2);

N_eigenvalues  = length(Solution.Eigenvalues);

% Initializations
Absolute_Error = zeros([size(A, 1) N_eigenvalues]);

Continuity_Absolute_Error = zeros([Ny Nx N_eigenvalues]);
X_Momentum_Absolute_Error = zeros([Ny Nx N_eigenvalues]);
Y_Momentum_Absolute_Error = zeros([Ny Nx N_eigenvalues]);
Z_Momentum_Absolute_Error = zeros([Ny Nx N_eigenvalues]);

Maximum_Continuity_Absolute_Error = zeros([1 N_eigenvalues]);
Maximum_X_Momentum_Absolute_Error = zeros([1 N_eigenvalues]);
Maximum_Y_Momentum_Absolute_Error = zeros([1 N_eigenvalues]);
Maximum_Z_Momentum_Absolute_Error = zeros([1 N_eigenvalues]);

ind_Maximum_Continuity_Absolute_Error = zeros([1 N_eigenvalues]);
ind_Maximum_X_Momentum_Absolute_Error = zeros([1 N_eigenvalues]);
ind_Maximum_Y_Momentum_Absolute_Error = zeros([1 N_eigenvalues]);
ind_Maximum_Z_Momentum_Absolute_Error = zeros([1 N_eigenvalues]);

Maximum_Continuity_Absolute_Error_X = zeros([1 N_eigenvalues]);
Maximum_Continuity_Absolute_Error_Y = zeros([1 N_eigenvalues]);
Maximum_X_Momentum_Absolute_Error_X = zeros([1 N_eigenvalues]);
Maximum_X_Momentum_Absolute_Error_Y = zeros([1 N_eigenvalues]);
Maximum_Y_Momentum_Absolute_Error_X = zeros([1 N_eigenvalues]);
Maximum_Y_Momentum_Absolute_Error_Y = zeros([1 N_eigenvalues]);
Maximum_Z_Momentum_Absolute_Error_X = zeros([1 N_eigenvalues]);
Maximum_Z_Momentum_Absolute_Error_Y = zeros([1 N_eigenvalues]);

for i = 1:N_eigenvalues
    omega = Solution.Eigenvalues(i);
    q = [Solution.Eigenfunctions.u(:,i) ; Solution.Eigenfunctions.v(:,i) ; Solution.Eigenfunctions.w(:,i) ; Solution.Eigenfunctions.p(:,i)]*Solution.Normalizers(i);
    LHS_matrix = A*q;
    RHS_matrix = omega*B*q;
    Absolute_Error(:,i) = abs((LHS_matrix - RHS_matrix)./LHS_matrix);
    
    X_Momentum_Absolute_Error(:,:,i) = reshape(Absolute_Error(0*Nx*Ny+1:1*Nx*Ny,i), Ny, Nx);
    Y_Momentum_Absolute_Error(:,:,i) = reshape(Absolute_Error(1*Nx*Ny+1:2*Nx*Ny,i), Ny, Nx);
    Z_Momentum_Absolute_Error(:,:,i) = reshape(Absolute_Error(2*Nx*Ny+1:3*Nx*Ny,i), Ny, Nx);
    Continuity_Absolute_Error(:,:,i) = reshape(Absolute_Error(3*Nx*Ny+1:4*Nx*Ny,i), Ny, Nx);
    
    [Maximum_Continuity_Absolute_Error(i), ind_Maximum_Continuity_Absolute_Error(i)] = max(Continuity_Absolute_Error(:,:,i), [], 'all');
    [Maximum_X_Momentum_Absolute_Error(i), ind_Maximum_X_Momentum_Absolute_Error(i)] = max(X_Momentum_Absolute_Error(:,:,i), [], 'all');
    [Maximum_Y_Momentum_Absolute_Error(i), ind_Maximum_Y_Momentum_Absolute_Error(i)] = max(Y_Momentum_Absolute_Error(:,:,i), [], 'all');
    [Maximum_Z_Momentum_Absolute_Error(i), ind_Maximum_Z_Momentum_Absolute_Error(i)] = max(Z_Momentum_Absolute_Error(:,:,i), [], 'all');
    
    Maximum_Continuity_Absolute_Error_X(i) = Domain.mat_X(ind_Maximum_Continuity_Absolute_Error(i));
    Maximum_Continuity_Absolute_Error_Y(i) = Domain.mat_Y(ind_Maximum_Continuity_Absolute_Error(i));
    Maximum_X_Momentum_Absolute_Error_X(i) = Domain.mat_X(ind_Maximum_X_Momentum_Absolute_Error(i));
    Maximum_X_Momentum_Absolute_Error_Y(i) = Domain.mat_Y(ind_Maximum_X_Momentum_Absolute_Error(i));
    Maximum_Y_Momentum_Absolute_Error_X(i) = Domain.mat_X(ind_Maximum_Y_Momentum_Absolute_Error(i));
    Maximum_Y_Momentum_Absolute_Error_Y(i) = Domain.mat_Y(ind_Maximum_Y_Momentum_Absolute_Error(i));
    Maximum_Z_Momentum_Absolute_Error_X(i) = Domain.mat_X(ind_Maximum_Z_Momentum_Absolute_Error(i));
    Maximum_Z_Momentum_Absolute_Error_Y(i) = Domain.mat_Y(ind_Maximum_Z_Momentum_Absolute_Error(i));
end

Navier_Stokes_Check.A = sparse(A);
Navier_Stokes_Check.B = sparse(B);

Navier_Stokes_Check.Errors.Continuity.Absolute.Map = Continuity_Absolute_Error;
Navier_Stokes_Check.Errors.X_Momentum.Absolute.Map = X_Momentum_Absolute_Error;
Navier_Stokes_Check.Errors.Y_Momentum.Absolute.Map = Y_Momentum_Absolute_Error;
Navier_Stokes_Check.Errors.Z_Momentum.Absolute.Map = Z_Momentum_Absolute_Error;

Navier_Stokes_Check.Errors.Continuity.Absolute.Maximum = Maximum_Continuity_Absolute_Error;
Navier_Stokes_Check.Errors.X_Momentum.Absolute.Maximum = Maximum_X_Momentum_Absolute_Error;
Navier_Stokes_Check.Errors.Y_Momentum.Absolute.Maximum = Maximum_Y_Momentum_Absolute_Error;
Navier_Stokes_Check.Errors.Z_Momentum.Absolute.Maximum = Maximum_Z_Momentum_Absolute_Error;

Navier_Stokes_Check.Errors.Continuity.Absolute.Maximum_X = Maximum_Continuity_Absolute_Error_X;
Navier_Stokes_Check.Errors.X_Momentum.Absolute.Maximum_X = Maximum_X_Momentum_Absolute_Error_X;
Navier_Stokes_Check.Errors.Y_Momentum.Absolute.Maximum_X = Maximum_Y_Momentum_Absolute_Error_X;
Navier_Stokes_Check.Errors.Z_Momentum.Absolute.Maximum_X = Maximum_Z_Momentum_Absolute_Error_X;

Navier_Stokes_Check.Errors.Continuity.Absolute.Maximum_Y = Maximum_Continuity_Absolute_Error_Y;
Navier_Stokes_Check.Errors.X_Momentum.Absolute.Maximum_Y = Maximum_X_Momentum_Absolute_Error_Y;
Navier_Stokes_Check.Errors.Y_Momentum.Absolute.Maximum_Y = Maximum_Y_Momentum_Absolute_Error_Y;
Navier_Stokes_Check.Errors.Z_Momentum.Absolute.Maximum_Y = Maximum_Z_Momentum_Absolute_Error_Y;


end
