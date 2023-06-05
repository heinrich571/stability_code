function [mat_A, mat_B] = create_eigenvalue_matrices(Domain, Base_Flow)

% Initialize matrix elements
initializer = zeros(size(Domain.Dx));
Z   = initializer;
U   = initializer;
V   = initializer;
Ux  = initializer;
Uy  = initializer;
Vx  = initializer;
Vy  = initializer;
L2D = initializer;
I   = eye(Ny+1, Nx+1);

% Compute matrix terms
Ux  = ;
Uy  = ;
Vx  = ;
Vy  = ;
L2D = U*Dx + V*Dy - 1/Re*(D2x + D2y - beta^2*I);

% LHS matrix entries
a11 = L2D + Ux;
a12 = Uy;
a13 = Z;
a14 = Dx;

a21 = Vx;
a22 = L2D + Vy;
a23 = Z;
a24 = Dy;

a31 = Z;
a32 = Z;
a33 = L2D;
a34 = 1i*beta*I;

a41 = Dx;
a42 = Dy;
a43 = 1i*beta*I;
a44 = Z;

% RHS matrix entries
b11 = 1i*I;
b12 = Z;
b13 = Z;
b14 = Z;

b21 = Z;
b22 = 1i*I;
b23 = Z;
b24 = Z;

b31 = Z;
b32 = Z;
b33 = 1i*I;
b34 = Z;

b41 = Z;
b42 = Z;
b43 = Z;
b44 = Z;

% Build matrices
mat_A = [a11 a12 a13 a14
         a21 a22 a23 a24
         a31 a32 a33 a34
         a41 a42 a43 a44];

mat_B = [b11 b12 b13 b14
         b21 b22 b23 b24
         b31 b32 b33 b34
         b41 b42 b43 b44];

end