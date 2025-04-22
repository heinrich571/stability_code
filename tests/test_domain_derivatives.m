%% Fresh Start

% close all
clear all
clc

addpath('../grid_generation/')


%% Domain Generation
% Define domain
Problem.Domain.Nx       = 60;
Problem.Domain.Ny       = 60;
Problem.Domain.X_Limit  = 200;
Problem.Domain.Y_Limit  = 300;
Problem.Domain.Y_Median = 3 * 2.4;

% Generate domain
Domain = generate_domain(Problem);


%% Build Test Function

casenum = 2;

xmat = Domain.mat_X;
ymat = Domain.mat_Y;

analytical = get_analytical_functions(xmat, ymat, casenum);
f = analytical.f;


%% Test Derivatives
% Numerical partial derivatives
nx   = length(Domain.vec_X);
ny   = length(Domain.vec_Y);
Dx   = Domain.Dx;
Dy   = Domain.Dy;
D2x  = Domain.D2x;
D2y  = Domain.D2y;
DxDy = Dx*Dy;
DyDx = Dy*Dx;

numerical.dfdx = reshape(Dx * f(:), [ny nx]);
numerical.dfdy = reshape(Dy * f(:), [ny nx]);

numerical.d2fdx2 = reshape(D2x * f(:), [ny nx]);
numerical.d2fdy2 = reshape(D2y * f(:), [ny nx]);

numerical.d2fdxdy = reshape(DxDy * f(:), [ny nx]);
numerical.d2fdydx = reshape(DyDx * f(:), [ny nx]);


%% Compute Errors

error.dfdx = numerical.dfdx - analytical.dfdx;
error.dfdy = numerical.dfdy - analytical.dfdy;

error.d2fdx2 = numerical.d2fdx2 - analytical.d2fdx2;
error.d2fdy2 = numerical.d2fdy2 - analytical.d2fdy2;

error.d2fdxdy = numerical.d2fdxdy - analytical.d2fdxdy;
error.d2fdydx = numerical.d2fdydx - analytical.d2fdydx;


%% Plot results

% f
figure('Name', 'f(x,y)', 'NumberTitle', 'off')
surf(xmat, ymat, f)
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$f\left(x,y\right)$', 'Interpreter', 'Latex')

% dfdx
figure('Name', 'dfdx', 'NumberTitle', 'off')
subplot(1,2,1)
surf(xmat, ymat, analytical.dfdx)
title('Analytical')
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\left(\frac{df}{dx}\right)_{\mathrm{analytical}}$', 'Interpreter', 'Latex')
view(-45,15)
ax1 = gca;
subplot(1,2,2)
surf(xmat, ymat, numerical.dfdx)
title('Numerical')
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\left(\frac{df}{dx}\right)_{\mathrm{numerical}}$', 'Interpreter', 'Latex')
view(-45,15)
ax2 = gca;
linkaxes([ax1 ax2])

% dfdx - error
figure('Name', 'err dfdx', 'NumberTitle', 'off')
surf(xmat, ymat, error.dfdx)
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\varepsilon_{df/dx}$', 'Interpreter', 'Latex')
view(-45,15)

% dfdy
figure('Name', 'dfdy', 'NumberTitle', 'off')
subplot(1,2,1)
surf(xmat, ymat, analytical.dfdy)
title('Analytical')
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\left(\frac{df}{dy}\right)_{\mathrm{analytical}}$', 'Interpreter', 'Latex')
view(-45,15)
ax1 = gca;
subplot(1,2,2)
surf(xmat, ymat, numerical.dfdy)
title('Numerical')
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\left(\frac{df}{dy}\right)_{\mathrm{numerical}}$', 'Interpreter', 'Latex')
view(-45,15)
ax2 = gca;
linkaxes([ax1 ax2])

% dfdy - error
figure('Name', 'err dfdy', 'NumberTitle', 'off')
surf(xmat, ymat, error.dfdy)
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\varepsilon_{df/dy}$', 'Interpreter', 'Latex')
view(-45,15)

% d2fdx2
figure('Name', 'd2fdx2', 'NumberTitle', 'off')
subplot(1,2,1)
surf(xmat, ymat, analytical.d2fdx2)
title('Analytical')
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\left((\frac{d^2f}{dx^2}\right)_{\mathrm{analytical}}$', 'Interpreter', 'Latex')
view(-45,15)
ax1 = gca;
subplot(1,2,2)
surf(xmat, ymat, numerical.d2fdx2)
title('Numerical')
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\left(\frac{d^2f}{dx^2}\right)_{\mathrm{numerical}}$', 'Interpreter', 'Latex')
view(-45,15)
ax2 = gca;
linkaxes([ax1 ax2])

% d2fdx2 - error
figure('Name', 'err d2fdx2', 'NumberTitle', 'off')
surf(xmat, ymat, error.d2fdx2)
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\varepsilon_{d^2f/dx^2}$', 'Interpreter', 'Latex')
view(-45,15)

% d2fdy2
figure('Name', 'd2fdy2', 'NumberTitle', 'off')
subplot(1,2,1)
surf(xmat, ymat, analytical.d2fdy2)
title('Analytical')
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\left(\frac{d^2f}{dy^2}\right)_{\mathrm{analytical}}$', 'Interpreter', 'Latex')
view(-45,15)
ax1 = gca;
subplot(1,2,2)
surf(xmat, ymat, numerical.d2fdy2)
title('Numerical')
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\left(\frac{d^2f}{dy^2}\right)_{\mathrm{numerical}}$', 'Interpreter', 'Latex')
view(-45,15)
ax2 = gca;
linkaxes([ax1 ax2])

% d2fdy2 - error
figure('Name', 'err d2fdy2', 'NumberTitle', 'off')
surf(xmat, ymat, error.d2fdy2)
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\varepsilon_{d2f/dy2}$', 'Interpreter', 'Latex')
view(-45,15)

% d2fdxdy
figure('Name', 'd2fdxdy', 'NumberTitle', 'off')
subplot(1,2,1)
surf(xmat, ymat, analytical.d2fdxdy)
title('Analytical')
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\left(\frac{d^2f}{dxdy}\right)_{\mathrm{analytical}}$', 'Interpreter', 'Latex')
view(-45,15)
ax1 = gca;
subplot(1,2,2)
surf(xmat, ymat, numerical.d2fdxdy)
title('Numerical')
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\left(\frac{d^2f}{dxdy}\right)_{\mathrm{numerical}}$', 'Interpreter', 'Latex')
view(-45,15)
ax2 = gca;
linkaxes([ax1 ax2])

% d2fdxdy - error
figure('Name', 'err d2fdxdy', 'NumberTitle', 'off')
surf(xmat, ymat, error.d2fdxdy)
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\varepsilon_{d2f/dxdy}$', 'Interpreter', 'Latex')
view(-45,15)

% d2fdydx
figure('Name', 'd2fdydx', 'NumberTitle', 'off')
subplot(1,2,1)
surf(xmat, ymat, analytical.d2fdydx)
title('Analytical')
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\left(\frac{d^2f}{dydx}\right)_{\mathrm{analytical}}$', 'Interpreter', 'Latex')
view(-45,15)
ax1 = gca;
subplot(1,2,2)
surf(xmat, ymat, numerical.d2fdydx)
title('Numerical')
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\left(\frac{d^2f}{dxdy}\right)_{\mathrm{numerical}}$', 'Interpreter', 'Latex')
view(-45,15)
ax2 = gca;
linkaxes([ax1 ax2])

% d2fdydx - error
figure('Name', 'err d2fdydx', 'NumberTitle', 'off')
surf(xmat, ymat, error.d2fdydx)
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\varepsilon_{d2f/dydx}$', 'Interpreter', 'Latex')
view(-45,15)


%% Cleanup

rmpath('../grid_generation/')


%% Supprting Functions

function analytical = get_analytical_functions(xmat, ymat, casenum)

switch casenum
    case 1
        f = sin(xmat).*cos(ymat);

        analytical.f = f;

        analytical.dfdx =  cos(xmat).*cos(ymat);
        analytical.dfdy = -sin(xmat).*sin(ymat);

        analytical.d2fdx2 = -f;
        analytical.d2fdy2 = -f;

        analytical.d2fdxdy = -cos(xmat).*sin(ymat);
        analytical.d2fdydx = analytical.d2fdxdy;
    case 2
        a = 1;
        b = -0.5;
        c = 0.02;
        f = a*exp(b*ymat.^2).*cos(c*xmat);

        analytical.f = f;

        analytical.dfdx = -a*c*exp(b*ymat.^2).*sin(c*xmat);
        analytical.dfdy =  2*a*b*ymat.*exp(b*ymat.^2).*cos(c*xmat);

        analytical.d2fdx2 = -a*c^2*exp(b*ymat.^2).*cos(c*xmat);
        analytical.d2fdy2 =  2*a*b*exp(b*ymat.^2).*(2*b*ymat.^2 + 1).*cos(c*xmat);

        analytical.d2fdxdy = -2*a*b*c*ymat.*exp(b*ymat.^2).*sin(c*xmat);
        analytical.d2fdydx =  analytical.d2fdxdy;
    otherwise
        error("get_analytical_functions - invalid 'casenum' value");
end

end
