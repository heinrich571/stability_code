function [A, B] = bglnsematrices(Domain, Base_Flow, Problem)

% Get spanwise wave number
beta = Problem.Physics.Beta;

% Build zeros and identity matrices
Ny = size(Domain.mat_X, 1);
Nx = size(Domain.mat_X, 2);
Z  = sparse(zeros(Nx * Ny));
I  = sparse(eye(Nx * Ny));

% Get base flow solution data
phi   = repmat(flip(Base_Flow.phi), [1 Nx]);
dphi  = repmat(flip(Base_Flow.dphi), [1 Nx]);
ddphi = repmat(flip(Base_Flow.ddphi), [1 Nx]);

% Get differentiation matrices
Dx    = Domain.Dx;
Dy    = Domain.Dy;
D2x   = Domain.D2x;
D2y   = Domain.D2y;
xmat  = Domain.mat_X;

% Various matrix elements
U  = sparse(xmat .* dphi);
U  = sparse(diag(U(:), 0));
V  = sparse(diag(-phi(:), 0));
Ux = sparse(diag(dphi(:), 0));
Uy = sparse(xmat .* ddphi);
Uy = sparse(diag(Uy(:), 0));
Vx = sparse(Z);
Vy = sparse(diag(-dphi(:), 0));
L  = sparse((D2x + D2y - beta^2*I) - U*Dx - V*Dy);

% Matrix entries
% x-momentum
a11 = L - Ux; a12 = -Uy; a13 = Z; a14 = -Dx;    b11 = -1i*I; b12 = Z; b13 = Z; b14 = Z;

% y-momentum
a21 = -Vx; a22 = L - Vy; a23 = Z; a24 = -Dy;    b21 = Z; b22 = -1i*I; b23 = Z; b24 = Z;

% z-momentum
a31 = Z; a32 = Z; a33 = L; a34 = -1i*beta*I;    b31 = Z; b32 = Z; b33 = -1i*I; b34 = Z;

% Continuity
a41 = Dx; a42 = Dy; a43 = 1i*beta*I; a44 = Z;   b41 = Z; b42 = Z;  b43 = Z; b44 = Z;

% Build EVP A matrix
A = [a11 a12 a13 a14
     a21 a22 a23 a24
     a31 a32 a33 a34
     a41 a42 a43 a44];

% Build EVP B matrix
B = [b11 b12 b13 b14
     b21 b22 b23 b24
     b31 b32 b33 b34
     b41 b42 b43 b44];

end
