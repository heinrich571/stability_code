%% Fresh Start

close all
clear all
clc

path_manager('add')


%% Domain Generation
% Define domain
Problem.Domain.Nx       = 40;
Problem.Domain.Ny       = 40;
Problem.Domain.X_Limit  = 100;
Problem.Domain.Y_Limit  = 50;
Problem.Domain.Y_Median = 3 * 2.4;

% Generate domain
Domain = generate_domain(Problem);


%% Build Test Function

casenum = 3;

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

% ymax_for_plot = 5;
% ymin_for_plot = 0;
ymax_for_plot = inf;
ymin_for_plot = -inf;
% ymax_for_plot =  Problem.Domain.Y_Limit;
% ymin_for_plot = -Problem.Domain.Y_Limit;

% f
figure('Name', 'f(x,y)', 'NumberTitle', 'off')
surf(xmat, ymat, f)
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$f\left(x,y\right)$', 'Interpreter', 'Latex')
ylim([ymin_for_plot ymax_for_plot])
view(-140,45)

% dfdx
figure('Name', 'dfdx', 'NumberTitle', 'off')
subplot(1,2,1)
surf(xmat, ymat, analytical.dfdx)
title('Analytical')
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\left(\frac{df}{dx}\right)_{\mathrm{analytical}}$', 'Interpreter', 'Latex')
ylim([ymin_for_plot ymax_for_plot])
view(-45,15)
ax1 = gca;
subplot(1,2,2)
surf(xmat, ymat, numerical.dfdx)
title('Numerical')
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\left(\frac{df}{dx}\right)_{\mathrm{numerical}}$', 'Interpreter', 'Latex')
ylim([ymin_for_plot ymax_for_plot])
view(-45,15)
ax2 = gca;
linkaxes([ax1 ax2])
linkprop([ax1 ax2], {'CameraPosition','CameraUpVector'});

% dfdx - error
figure('Name', 'err dfdx', 'NumberTitle', 'off')
surf(xmat, ymat, error.dfdx)
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\varepsilon_{df/dx}$', 'Interpreter', 'Latex')
ylim([ymin_for_plot ymax_for_plot])
view(-45,15)

% dfdy
figure('Name', 'dfdy', 'NumberTitle', 'off')
subplot(1,2,1)
surf(xmat, ymat, analytical.dfdy)
title('Analytical')
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\left(\frac{df}{dy}\right)_{\mathrm{analytical}}$', 'Interpreter', 'Latex')
ylim([ymin_for_plot ymax_for_plot])
view(-45,15)
ax1 = gca;
subplot(1,2,2)
surf(xmat, ymat, numerical.dfdy)
title('Numerical')
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\left(\frac{df}{dy}\right)_{\mathrm{numerical}}$', 'Interpreter', 'Latex')
ylim([ymin_for_plot ymax_for_plot])
view(-45,15)
ax2 = gca;
linkaxes([ax1 ax2])
linkprop([ax1 ax2], {'CameraPosition','CameraUpVector'});

% dfdy - error
figure('Name', 'err dfdy', 'NumberTitle', 'off')
surf(xmat, ymat, error.dfdy)
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\varepsilon_{df/dy}$', 'Interpreter', 'Latex')
ylim([ymin_for_plot ymax_for_plot])
view(-45,15)

% d2fdx2
figure('Name', 'd2fdx2', 'NumberTitle', 'off')
subplot(1,2,1)
surf(xmat, ymat, analytical.d2fdx2)
title('Analytical')
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\left((\frac{d^2f}{dx^2}\right)_{\mathrm{analytical}}$', 'Interpreter', 'Latex')
ylim([ymin_for_plot ymax_for_plot])
view(-45,15)
ax1 = gca;
subplot(1,2,2)
surf(xmat, ymat, numerical.d2fdx2)
title('Numerical')
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\left(\frac{d^2f}{dx^2}\right)_{\mathrm{numerical}}$', 'Interpreter', 'Latex')
ylim([ymin_for_plot ymax_for_plot])
view(-45,15)
ax2 = gca;
linkaxes([ax1 ax2])
linkprop([ax1 ax2], {'CameraPosition','CameraUpVector'});

% d2fdx2 - error
figure('Name', 'err d2fdx2', 'NumberTitle', 'off')
surf(xmat, ymat, error.d2fdx2)
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\varepsilon_{d^2f/dx^2}$', 'Interpreter', 'Latex')
ylim([ymin_for_plot ymax_for_plot])
view(-45,15)

% d2fdy2
figure('Name', 'd2fdy2', 'NumberTitle', 'off')
subplot(1,2,1)
surf(xmat, ymat, analytical.d2fdy2)
title('Analytical')
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\left(\frac{d^2f}{dy^2}\right)_{\mathrm{analytical}}$', 'Interpreter', 'Latex')
ylim([ymin_for_plot ymax_for_plot])
view(-45,15)
ax1 = gca;
subplot(1,2,2)
surf(xmat, ymat, numerical.d2fdy2)
title('Numerical')
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\left(\frac{d^2f}{dy^2}\right)_{\mathrm{numerical}}$', 'Interpreter', 'Latex')
ylim([ymin_for_plot ymax_for_plot])
view(-45,15)
ax2 = gca;
linkaxes([ax1 ax2])
linkprop([ax1 ax2], {'CameraPosition','CameraUpVector'});

% d2fdy2 - error
figure('Name', 'err d2fdy2', 'NumberTitle', 'off')
surf(xmat, ymat, error.d2fdy2)
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\varepsilon_{d2f/dy2}$', 'Interpreter', 'Latex')
ylim([ymin_for_plot ymax_for_plot])
view(-45,15)

% d2fdxdy
figure('Name', 'd2fdxdy', 'NumberTitle', 'off')
subplot(1,2,1)
surf(xmat, ymat, analytical.d2fdxdy)
title('Analytical')
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\left(\frac{d^2f}{dxdy}\right)_{\mathrm{analytical}}$', 'Interpreter', 'Latex')
ylim([ymin_for_plot ymax_for_plot])
view(-45,15)
ax1 = gca;
subplot(1,2,2)
surf(xmat, ymat, numerical.d2fdxdy)
title('Numerical')
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\left(\frac{d^2f}{dxdy}\right)_{\mathrm{numerical}}$', 'Interpreter', 'Latex')
ylim([ymin_for_plot ymax_for_plot])
view(-45,15)
ax2 = gca;
linkaxes([ax1 ax2])
linkprop([ax1 ax2], {'CameraPosition','CameraUpVector'});

% d2fdxdy - error
figure('Name', 'err d2fdxdy', 'NumberTitle', 'off')
surf(xmat, ymat, error.d2fdxdy)
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\varepsilon_{d2f/dxdy}$', 'Interpreter', 'Latex')
ylim([ymin_for_plot ymax_for_plot])
view(-45,15)

% d2fdydx
figure('Name', 'd2fdydx', 'NumberTitle', 'off')
subplot(1,2,1)
surf(xmat, ymat, analytical.d2fdydx)
title('Analytical')
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\left(\frac{d^2f}{dydx}\right)_{\mathrm{analytical}}$', 'Interpreter', 'Latex')
ylim([ymin_for_plot ymax_for_plot])
view(-45,15)
ax1 = gca;
subplot(1,2,2)
surf(xmat, ymat, numerical.d2fdydx)
title('Numerical')
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\left(\frac{d^2f}{dxdy}\right)_{\mathrm{numerical}}$', 'Interpreter', 'Latex')
ylim([ymin_for_plot ymax_for_plot])
view(-45,15)
ax2 = gca;
linkaxes([ax1 ax2])
linkprop([ax1 ax2], {'CameraPosition','CameraUpVector'});

% d2fdydx - error
figure('Name', 'err d2fdydx', 'NumberTitle', 'off')
surf(xmat, ymat, error.d2fdydx)
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
zlabel('$\varepsilon_{d2f/dydx}$', 'Interpreter', 'Latex')
ylim([ymin_for_plot ymax_for_plot])
view(-45,15)


%% Cleanup

path_manager('remove')


%% Supprting Functions

function analytical = get_analytical_functions(xmat, ymat, casenum)

switch casenum
    case 1
        wx = 0.01;
        wy = 0.01;
        f = sin(wx * xmat).*cos(wy * ymat);

        analytical.f = f;

        analytical.dfdx =  wx *  cos(wx * xmat).*cos(wy * ymat);
        analytical.dfdy =  wy * -sin(wx * xmat).*sin(wy * ymat);

        analytical.d2fdx2 = wx^2 * -f;
        analytical.d2fdy2 = wy^2 * -f;

        analytical.d2fdxdy = wx * wy * -cos(wx * xmat).*sin(wy * ymat);
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
    case 3
        Base_Flow_Definitions.initguess            = 1.23258765682022 + [-1 1]*1e-5;
        Base_Flow_Definitions.maxIterations        = 1e2;
        Base_Flow_Definitions.convergenceTolerance = 1e-6;
        
        Base_Flow_Definitions.interval = ymat(:,1);
        Base_Flow = get_base_flow(Base_Flow_Definitions);

        nx     = length(xmat(1,:));
        zmat   = zeros(size(xmat));
        phi    = flip(Base_Flow.phi);
        dphi   = flip(Base_Flow.dphi);
        ddphi  = flip(Base_Flow.ddphi);
        mat_phi   = repmat(phi  , [1 nx]);
        mat_dphi  = repmat(dphi , [1 nx]);
        mat_ddphi = repmat(ddphi, [1 nx]);
        
        U =  xmat .* mat_dphi;
        V = -mat_phi;

        Ux =  mat_dphi;
        Uy =  xmat .* mat_ddphi;
        Vx =  zmat;
        Vy = -mat_dphi;

        Uxx =  zmat;
        Uyy =  zmat; % THIS IS NOT TRUE!!!
        Vxx =  zmat;
        Vyy = -mat_ddphi;

        Uxy = mat_ddphi;
        Vxy = zmat;

        % U
        analytical.f = U;

        analytical.dfdx = Ux;
        analytical.dfdy = Uy;

        analytical.d2fdx2 = Uxx;
        analytical.d2fdy2 = Uyy;

        analytical.d2fdxdy = Uxy;
        analytical.d2fdydx = Uxy;

        % % V
        % analytical.f = V;
        % 
        % analytical.dfdx = Vx;
        % analytical.dfdy = Vy;
        % 
        % analytical.d2fdx2 = Vxx;
        % analytical.d2fdy2 = Vyy;
        % 
        % analytical.d2fdxdy = Vxy;
        % analytical.d2fdydx = Vxy;
    otherwise
        error("get_analytical_functions - invalid 'casenum' value");
end

end
