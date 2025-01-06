%% Fresh Start

% close all
clear all
clc

addpath('../grid_generation/')


%% Domain Generation
% Define domain
Nx = 40;
Ny = 120;
Problem.Domain.Nx       = Nx;
Problem.Domain.Ny       = Ny;
Problem.Domain.X_Limit  = 200;
Problem.Domain.Y_Limit  = 300;
Problem.Domain.Y_Median = 6*2.4;

% Generate domain
Domain = generate_domain(Problem);


%% Build Test Function

mat_X = Domain.mat_X;
mat_Y = Domain.mat_Y;

% f = sin(mat_X).*cos(mat_Y);
% 
% analytical.dfdx =  cos(mat_X).*cos(mat_Y);
% analytical.dfdy = -sin(mat_X).*sin(mat_Y);
% 
% analytical.d2fdx2 = -f;
% analytical.d2fdy2 = -f;
% 
% analytical.d2fdxdy = -cos(mat_X).*sin(mat_Y);
% analytical.d2fdydx = analytical.d2fdxdy;

a = 1;
b = -0.5;
c = 0.02;
f = a*exp(b*mat_Y.^2).*cos(c*mat_X);

analytical.dfdx = -a*c*exp(b*mat_Y.^2).*sin(c*mat_X);
analytical.dfdy = 2*a*b*mat_Y.*exp(b*mat_Y.^2).*cos(c*mat_X);

analytical.d2fdx2 = -a*c^2*exp(b*mat_Y.^2).*cos(c*mat_X);
analytical.d2fdy2 = 2*a*b*exp(b*mat_Y.^2).*(2*b*mat_Y.^2 + 1).*cos(c*mat_X);

analytical.d2fdxdy = -2*a*b*c*mat_Y.*exp(b*mat_Y.^2).*sin(c*mat_X);
analytical.d2fdydx = analytical.d2fdxdy;

surf(mat_X, mat_Y, f)


%% Test Derivatives
% Numerical partial derivatives
Dx   = Domain.Dx;
Dy   = Domain.Dy;
D2x  = Domain.D2x;
D2y  = Domain.D2y;
DxDy = Dx*Dy;
DyDx = Dy*Dx;
Nx = length(Domain.vec_X);
Ny = length(Domain.vec_Y);


numerical.dfdx = reshape(Dx*f(:), [Ny Nx]);
numerical.dfdy = reshape(Dy*f(:), [Ny Nx]);

numerical.d2fdx2 = reshape(D2x*f(:), [Ny Nx]);
numerical.d2fdy2 = reshape(D2y*f(:), [Ny Nx]);

numerical.d2fdxdy = reshape(DxDy*f(:), [Ny Nx]);
numerical.d2fdydx = reshape(DyDx*f(:), [Ny Nx]);


%% Compute Errors

error.dfdx = numerical.dfdx - analytical.dfdx;
error.dfdy = numerical.dfdy - analytical.dfdy;

error.d2fdx2 = numerical.d2fdx2 - analytical.d2fdx2;
error.d2fdy2 = numerical.d2fdy2 - analytical.d2fdy2;

error.d2fdxdy = numerical.d2fdxdy - analytical.d2fdxdy;
error.d2fdydx = numerical.d2fdydx - analytical.d2fdydx;


%% Plot results

% dfdx
figure('Name', 'dfdx', 'NumberTitle', 'off')
subplot(1,2,1)
surf(mat_X, mat_Y, analytical.dfdx)
title('Analytical')
xlabel('x')
ylabel('y')
zlabel('(df/dx)_{analytical}')
view(-45,15)
subplot(1,2,2)
surf(mat_X, mat_Y, numerical.dfdx)
title('Numerical')
xlabel('x')
ylabel('y')
zlabel('(df/dx)_{numerical}')
view(-45,15)

% dfdx - error
figure('Name', 'err dfdx', 'NumberTitle', 'off')
surf(mat_X, mat_Y, error.dfdx)
xlabel('x')
ylabel('y')
zlabel('$\varepsilon_{df/dx}$')
view(-45,15)

% dfdy
figure('Name', 'dfdy', 'NumberTitle', 'off')
subplot(1,2,1)
surf(mat_X, mat_Y, analytical.dfdy)
title('Analytical')
xlabel('x')
ylabel('y')
zlabel('(df/dy)_{analytical}')
view(-45,15)
subplot(1,2,2)
surf(mat_X, mat_Y, numerical.dfdy)
title('Numerical')
xlabel('x')
ylabel('y')
zlabel('(df/dy)_{numerical}')
view(-45,15)

% dfdy - error
figure('Name', 'err dfdy', 'NumberTitle', 'off')
surf(mat_X, mat_Y, error.dfdy)
xlabel('x')
ylabel('y')
zlabel('$\varepsilon_{df/dy}$')
view(-45,15)

% d2fdx2
figure('Name', 'd2fdx2', 'NumberTitle', 'off')
subplot(1,2,1)
surf(mat_X, mat_Y, analytical.d2fdx2)
title('Analytical')
xlabel('x')
ylabel('y')
zlabel('(d2f/dx2)_{analytical}')
view(-45,15)
subplot(1,2,2)
surf(mat_X, mat_Y, numerical.d2fdx2)
title('Numerical')
xlabel('x')
ylabel('y')
zlabel('(d2f/dx2)_{numerical}')
view(-45,15)

% d2fdx2 - error
figure('Name', 'err d2fdx2', 'NumberTitle', 'off')
surf(mat_X, mat_Y, error.d2fdx2)
xlabel('x')
ylabel('y')
zlabel('$\varepsilon_{d2f/dx2}$')
view(-45,15)

% d2fdy2
figure('Name', 'd2fdy2', 'NumberTitle', 'off')
subplot(1,2,1)
surf(mat_X, mat_Y, analytical.d2fdy2)
title('Analytical')
xlabel('x')
ylabel('y')
zlabel('(d2f/dy2)_{analytical}')
view(-45,15)
subplot(1,2,2)
surf(mat_X, mat_Y, numerical.d2fdy2)
title('Numerical')
xlabel('x')
ylabel('y')
zlabel('(d2f/dy2)_{numerical}')
view(-45,15)

% d2fdy2 - error
figure('Name', 'err d2fdy2', 'NumberTitle', 'off')
surf(mat_X, mat_Y, error.d2fdy2)
xlabel('x')
ylabel('y')
zlabel('$\varepsilon_{d2f/dy2}$')
view(-45,15)

% d2fdxdy
figure('Name', 'd2fdxdy', 'NumberTitle', 'off')
subplot(1,2,1)
surf(mat_X, mat_Y, analytical.d2fdxdy)
title('Analytical')
xlabel('x')
ylabel('y')
zlabel('(d2f/dxdy)_{analytical}')
view(-45,15)
subplot(1,2,2)
surf(mat_X, mat_Y, numerical.d2fdxdy)
title('Numerical')
xlabel('x')
ylabel('y')
zlabel('(d2f/dxdy)_{numerical}')
view(-45,15)

% d2fdxdy - error
figure('Name', 'err d2fdxdy', 'NumberTitle', 'off')
surf(mat_X, mat_Y, error.d2fdxdy)
xlabel('x')
ylabel('y')
zlabel('$\varepsilon_{d2f/dxdy}$')
view(-45,15)

% d2fdydx
figure('Name', 'd2fdydx', 'NumberTitle', 'off')
subplot(1,2,1)
surf(mat_X, mat_Y, analytical.d2fdydx)
title('Analytical')
xlabel('x')
ylabel('y')
zlabel('(d2f/dydx)_{analytical}')
view(-45,15)
subplot(1,2,2)
surf(mat_X, mat_Y, numerical.d2fdydx)
title('Numerical')
xlabel('x')
ylabel('y')
zlabel('(d2f/dydx)_{numerical}')
view(-45,15)

% d2fdydx - error
figure('Name', 'err d2fdydx', 'NumberTitle', 'off')
surf(mat_X, mat_Y, error.d2fdydx)
xlabel('x')
ylabel('y')
zlabel('$\varepsilon_{d2f/dydx}$')
view(-45,15)


%% Cleanup

rmpath('../grid_generation/')