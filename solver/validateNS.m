function Navier_Stokes_Check = validateNS(Domain, Base_Flow, Problem, Solution)

%% Initialize matrix elements

initializer = zeros(size(Domain.Dx));
Ny  = size(Domain.mat_X, 1);
Nx  = size(Domain.mat_X, 2);
Z   = initializer;
I   = eye(Ny*Nx);


% Compute matrix terms
phi   = flip(Base_Flow.phi);
dphi  = flip(Base_Flow.dphi);
ddphi = flip(Base_Flow.ddphi);

mat_X = Domain.mat_X;
Dx    = Domain.Dx;
Dy    = Domain.Dy;
D2x   = Domain.D2x;
D2y   = Domain.D2y;

beta = Problem.Physics.Beta;

mat_phi   = repmat(phi  , [1,Nx]);
mat_dphi  = repmat(dphi , [1,Nx]);
mat_ddphi = repmat(ddphi, [1,Nx]);
U  = mat_X.*mat_dphi;
U  = diag(U(:), 0);
V  = diag(-mat_phi(:), 0);
Ux = diag(mat_dphi(:), 0);
Uy = mat_X.*mat_ddphi;
Uy = diag(Uy(:), 0);
Vx = Z;
Vy = diag(-mat_dphi(:), 0);
L  = U*Dx + V*Dy - (D2x + D2y - beta^2*I);


%% LHS matrix entries
% x-momentum
a11 = L + Ux;
a12 = Uy;
a13 = Z;
a14 = Dx;

% y-momentum
a21 = Vx;
a22 = L + Vy;
a23 = Z;
a24 = Dy;

% z-momentum
a31 = Z;
a32 = Z;
a33 = L;
a34 = 1i*beta*I;

% continuity
a41 = Dx;
a42 = Dy;
a43 = 1i*beta*I;
a44 = Z;


%% RHS matrix entries
% x-momentum
b11 = 1i*I;
b12 = Z;
b13 = Z;
b14 = Z;

% y-momentum

b21 = Z;
b22 = 1i*I;
b23 = Z;
b24 = Z;

% z-momentum
b31 = Z;
b32 = Z;
b33 = 1i*I;
b34 = Z;

% continuity
b41 = Z;
b42 = Z;
b43 = Z;
b44 = Z;


%% Build eigenvalue problem A and B matrices

mat_A = [a11 a12 a13 a14
         a21 a22 a23 a24
         a31 a32 a33 a34
         a41 a42 a43 a44];

mat_B = [b11 b12 b13 b14
         b21 b22 b23 b24
         b31 b32 b33 b34
         b41 b42 b43 b44];


%% Check the validity of the solution

N_eigenvalues  = length(Solution.Eigenvalues);

% Initializations
Absolute_Error = zeros([size(mat_A, 1) N_eigenvalues]);

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
    LHS_matrix = mat_A*q;
    RHS_matrix = omega*mat_B*q;
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

Navier_Stokes_Check.mat_A = mat_A;
Navier_Stokes_Check.mat_B = mat_B;

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
