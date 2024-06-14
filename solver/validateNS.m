function Report = validateNS(Domain, Base_Flow, Problem, Solution)

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

% continuity
a11 = Dx;
a12 = Dy;
a13 = 1i*beta*I;
a14 = Z;

% x-momentum
a21 = L + Ux;
a22 = Uy;
a23 = Z;
a24 = Dx;

% y-momentum
a31 = Vx;
a32 = L + Vy;
a33 = Z;
a34 = Dy;

% z-momentum
a41 = Z;
a42 = Z;
a43 = L;
a44 = 1i*beta*I;


%% RHS matrix entries

% continuity
b11 = Z;
b12 = Z;
b13 = Z;
b14 = Z;

% x-momentum
b21 = 1i*I;
b22 = Z;
b23 = Z;
b24 = Z;

% y-momentum

b31 = Z;
b32 = 1i*I;
b33 = Z;
b34 = Z;

% z-momentum
b41 = Z;
b42 = Z;
b43 = 1i*I;
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
    q = [Solution.Eigenfunctions.u(:,i) ; Solution.Eigenfunctions.v(:,i) ; Solution.Eigenfunctions.w(:,i) ; Solution.Eigenfunctions.p(:,i)];
    LHS_matrix = mat_A*q;
    RHS_matrix = omega*mat_B*q;
    Absolute_Error(:,i) = abs(LHS_matrix - RHS_matrix);
    
    Continuity_Absolute_Error(:,:,i) = reshape(Absolute_Error(0*Nx*Ny+1:1*Nx*Ny,i), Ny, Nx);
    X_Momentum_Absolute_Error(:,:,i) = reshape(Absolute_Error(1*Nx*Ny+1:2*Nx*Ny,i), Ny, Nx);
    Y_Momentum_Absolute_Error(:,:,i) = reshape(Absolute_Error(2*Nx*Ny+1:3*Nx*Ny,i), Ny, Nx);
    Z_Momentum_Absolute_Error(:,:,i) = reshape(Absolute_Error(3*Nx*Ny+1:4*Nx*Ny,i), Ny, Nx);
    
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

Report.Errors.Continuity.Absolute.Map = Continuity_Absolute_Error;
Report.Errors.X_Momentum.Absolute.Map = X_Momentum_Absolute_Error;
Report.Errors.Y_Momentum.Absolute.Map = Y_Momentum_Absolute_Error;
Report.Errors.Z_Momentum.Absolute.Map = Z_Momentum_Absolute_Error;

Report.Errors.Continuity.Absolute.Maximum = Maximum_Continuity_Absolute_Error;
Report.Errors.X_Momentum.Absolute.Maximum = Maximum_X_Momentum_Absolute_Error;
Report.Errors.Y_Momentum.Absolute.Maximum = Maximum_Y_Momentum_Absolute_Error;
Report.Errors.Z_Momentum.Absolute.Maximum = Maximum_Z_Momentum_Absolute_Error;

Report.Errors.Continuity.Absolute.Maximum_X = Maximum_Continuity_Absolute_Error_X;
Report.Errors.X_Momentum.Absolute.Maximum_X = Maximum_X_Momentum_Absolute_Error_X;
Report.Errors.Y_Momentum.Absolute.Maximum_X = Maximum_Y_Momentum_Absolute_Error_X;
Report.Errors.Z_Momentum.Absolute.Maximum_X = Maximum_Z_Momentum_Absolute_Error_X;
Report.Errors.Continuity.Absolute.Maximum_Y = Maximum_Continuity_Absolute_Error_Y;
Report.Errors.X_Momentum.Absolute.Maximum_Y = Maximum_X_Momentum_Absolute_Error_Y;
Report.Errors.Y_Momentum.Absolute.Maximum_Y = Maximum_Y_Momentum_Absolute_Error_Y;
Report.Errors.Z_Momentum.Absolute.Maximum_Y = Maximum_Z_Momentum_Absolute_Error_Y;


end
