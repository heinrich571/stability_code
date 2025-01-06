function [A, B] = bglnsematrices(Domain, Base_Flow, Problem)

% Initialize matrix elements
initializer = zeros(size(Domain.Dx));
Ny  = size(Domain.mat_X, 1);
Nx  = size(Domain.mat_X, 2);
Z   = initializer;
I   = eye(Ny*Nx);

% Compute matrix elements
phi   = flip(Base_Flow.phi);
dphi  = flip(Base_Flow.dphi);
ddphi = flip(Base_Flow.ddphi);

mat_X = Domain.mat_X;
Dx    = Domain.Dx;
Dy    = Domain.Dy;
D2x   = Domain.D2x;
D2y   = Domain.D2y;

beta = Problem.Physics.Beta;

mat_phi   = repmat(phi  , [1 Nx]);
mat_dphi  = repmat(dphi , [1 Nx]);
mat_ddphi = repmat(ddphi, [1 Nx]);
U  = mat_X.*mat_dphi;
U  = diag(U(:), 0);
V  = diag(-mat_phi(:), 0);
Ux = diag(mat_dphi(:), 0);
Uy = mat_X.*mat_ddphi;
Uy = diag(Uy(:), 0);
Vx = Z;
Vy = diag(-mat_dphi(:), 0);
L  = (D2x + D2y - beta^2*I) - U*Dx - V*Dy;


% LHS matrix entries
% x-momentum
a11 = L - Ux; a12 = -Uy; a13 = Z; a14 = -Dx;

% y-momentum
a21 = -Vx; a22 = L - Vy; a23 = Z; a24 = -Dy;

% z-momentum
a31 = Z; a32 = Z; a33 = L; a34 = -1i*beta*I;

% continuity
a41 = Dx; a42 = Dy; a43 = 1i*beta*I; a44 = Z;

% a44 = a44 + 1e-8*eye(Nx*Ny);

% RHS matrix entries
% x-momentum
b11 = -1i*I; b12 = Z; b13 = Z; b14 = Z;

% y-momentum
b21 = Z; b22 = -1i*I; b23 = Z; b24 = Z;

% z-momentum
b31 = Z; b32 = Z; b33 = -1i*I; b34 = Z;

% continuity
b41 = Z; b42 = Z;  b43 = Z; b44 = Z;


% Build eigenvalue problem A and B matrices
A = [a11 a12 a13 a14
     a21 a22 a23 a24
     a31 a32 a33 a34
     a41 a42 a43 a44];

B = [b11 b12 b13 b14
     b21 b22 b23 b24
     b31 b32 b33 b34
     b41 b42 b43 b44];

end
