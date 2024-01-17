function [mat_A, mat_B] = create_eigenvalue_matrices(Domain, Base_Flow, beta)

% Initialize matrix elements

initializer = zeros(size(Domain.Dx));
Ny  = size(Domain.mat_X, 1)-1;
Nx  = size(Domain.mat_X, 2)-1;
Z   = initializer;
I   = eye((Ny+1)*(Nx+1));


% Compute matrix terms

phi   = flip(Base_Flow.phi);
dphi  = flip(Base_Flow.dphi);
ddphi = flip(Base_Flow.ddphi);

mat_X = Domain.mat_X;
Dx    = Domain.Dx;
Dy    = Domain.Dy;

mat_phi    = repmat(phi  , [1,Nx+1]);
mat_dphi   = repmat(dphi , [1,Nx+1]);
mat_ddphi  = repmat(ddphi, [1,Nx+1]);
U   = mat_X.*mat_dphi;
U   = diag(U(:), 0);
V   = diag(-mat_phi(:), 0);
Ux  = diag(mat_dphi(:), 0);
Uy  = mat_X.*mat_ddphi;
Uy  = diag(Uy(:), 0);
Vx  = Z;
Vy  = diag(-mat_dphi(:), 0);
L   = (U - Dx)*Dx - (-V - Dy)*Dy + beta^2;


% LHS matrix entries

a11 = Dx;
a12 = Dy;
a13 = 1i*beta*I;
a14 = Z;

a21 = L + Ux;
a22 = Uy;
a23 = Z;
a24 = Dx;

a31 = Vx;
a32 = L + Vy;
a33 = Z;
a34 = Dy;

a41 = Z;
a42 = Z;
a43 = L;
a44 = 1i*beta*I;


% RHS matrix entries

b11 = Z;
b12 = Z;
b13 = Z;
b14 = Z;

b21 = 1i*I;
b22 = Z;
b23 = Z;
b24 = Z;

b31 = Z;
b32 = 1i*I;
b33 = Z;
b34 = Z;

b41 = Z;
b42 = Z;
b43 = 1i*I;
b44 = Z;


% Build output matrices

mat_A = [a11 a12 a13 a14
         a21 a22 a23 a24
         a31 a32 a33 a34
         a41 a42 a43 a44];

mat_B = [b11 b12 b13 b14
         b21 b22 b23 b24
         b31 b32 b33 b34
         b41 b42 b43 b44];


end
