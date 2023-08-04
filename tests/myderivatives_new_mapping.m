%% Fresh start

% close all
clearvars
clc


%% Define the computational domain

x_limit  = 2;
y_limit = 2;
y_median = 0.5;

Nx = 100;
Ny = Nx/2;


%% Generate the Chebyshev interval

[Dx_cheb, xhat] = cheb(Nx);
[Dy_cheb, yhat] = cheb(Ny);


%% Convert the Chebyshev interval to the physical domain

% x1 = xhat;

a_y = y_median*y_limit/(y_limit-2*y_median);
b_y = 1+2*a_y/y_limit;
y = (a_y*(1+yhat)./(b_y-yhat));

f = 0.01;
e = 1-f;
x1 = e*xhat.^3 + f*xhat;
ix_mid = Nx/2 + 1;
m = 1;
n = 0;
c = pi/2*(m-n);
d = pi/2*(m+n+1);

x = x1;
% x = x_limit*[cos(c*x1(1:ix_mid)+d)+1 ; -cos(c*x1(ix_mid+1:end)+d)-1];

[X, Y] = meshgrid(x, y);

% Dx_physical_domain = diag(-x_limit*sign(xhat).*pi/2.*(3*e*xhat.^2+f).*sin(pi/2*(e*xhat.^3+f*xhat)+pi))*Dx_cheb;
Dx_physical_domain = diag(1./(3*e*xhat.^2+f))*Dx_cheb;
Dy_physical_domain = diag(a_y*(1+b_y)./((y+a_y).^2))*Dy_cheb;
Ix = eye(Nx+1);
Iy = eye(Ny+1);

Dx = kron(Dx_physical_domain, Iy);
Dy = kron(Ix, Dy_physical_domain);


%% Plot results

figure('Name', 'New mapping', 'NumberTitle', 'off')
plot(X, Y, 'k', 'LineWidth', 0.5)
hold on
plot(X', Y', 'k', 'LineWidth', 0.5)
title('Computational Domain')
xlabel('x')
ylabel('y')
grid off

figure('Name', 'x mapping', 'NumberTitle', 'off')
plot(x, xhat, '-o')
ylim([-2 1])
YLimits = get(gca, 'YLim');
for i = 1:numel(xhat)
    plot(x(i)*ones([1 2]), [YLimits(1), xhat(i)], 'r')
end
xlabel('$x$', 'Interpreter', 'LaTeX', 'FontSize', 20)
ylabel('$\hat{x}$', 'Interpreter', 'LaTeX', 'FontSize', 20)

figure('Name', 'Intermediate step x mapping', 'NumberTitle', 'off')
plot(x1, xhat, '-o')
ylim([-2 1])
YLimits = get(gca, 'YLim');
for i = 1:numel(xhat)
    plot(x1(i)*ones([1 2]), [YLimits(1), xhat(i)], 'r')
end
xlabel('$x$', 'Interpreter', 'LaTeX', 'FontSize', 20)
ylabel('$\hat{x}$', 'Interpreter', 'LaTeX', 'FontSize', 20)

figure('Name', 'y mapping', 'NumberTitle', 'off')
plot(yhat, y, '-o')
xlim([-2 1])
XLimits = get(gca, 'XLim');
for i = 1:numel(yhat)
    plot([XLimits(1), yhat(i)], y(i)*ones([1 2]), 'r')
end
xlabel('$\hat{y}$', 'Interpreter', 'LaTeX', 'FontSize', 20)
ylabel('$y$', 'Interpreter', 'LaTeX', 'FontSize', 20)

% f = X.*Y;
% fx = Y;
% fy = X;
% f = f(:);

f   =  sin(X).*cos(Y);
fx  =  cos(X).*cos(Y);
fy  = -sin(X).*sin(Y);
fxx = -sin(X).*cos(Y);
fyy = -sin(X).*cos(Y);
f   =  f(:);

figure('Name', 'df/dx derivative test', 'NumberTitle', 'off')
subplot(1,2,1)
contour(X, Y, reshape(fx,Ny+1,Nx+1))
title('df/dx exact derivative')
xlabel('x')
ylabel('y')
grid off
subplot(1,2,2)
contour(X, Y, reshape(Dx*f,Ny+1,Nx+1))
title('df/dx numerical derivative')
xlabel('x')
ylabel('y')
grid off


figure('Name', 'df/dx derivative differences', 'NumberTitle', 'off')
surf(X, Y, reshape(fx,Ny+1,Nx+1)-reshape(Dx*f,Ny+1,Nx+1))
title('df/dx numerical derivative error')
xlabel('x')
ylabel('y')


figure('Name', 'df/dy derivative test', 'NumberTitle', 'off')
subplot(1,2,1)
contour(X, Y, reshape(fy,Ny+1,Nx+1))
title('df/dx exact derivative')
xlabel('x')
ylabel('y')
grid off
subplot(1,2,2)
contour(X, Y, reshape(Dy*f,Ny+1,Nx+1))
title('df/dx numerical derivative')
xlabel('x')
ylabel('y')
grid off


figure('Name', 'df/dy derivative differences', 'NumberTitle', 'off')
surf(X, Y, reshape(fy,Ny+1,Nx+1)-reshape(Dy*f,Ny+1,Nx+1))
title('df/dy numerical derivative error')
xlabel('x')
ylabel('y')


figure('Name', 'd2f/dx2 derivative test', 'NumberTitle', 'off')
subplot(1,2,1)
contour(X, Y, reshape(fxx,Ny+1,Nx+1))
title('d^2f/dx^2 exact derivative')
xlabel('x')
ylabel('y')
grid off
subplot(1,2,2)
contour(X, Y, reshape(Dx*Dx*f,Ny+1,Nx+1))
title('df/dx numerical derivative')
xlabel('x')
ylabel('y')
grid off


figure('Name', 'd^2f/dx^2 derivative differences', 'NumberTitle', 'off')
surf(X, Y, reshape(fxx,Ny+1,Nx+1)-reshape(Dx*Dx*f,Ny+1,Nx+1))
title('df/dx numerical derivative error')
xlabel('x')
ylabel('y')


figure('Name', 'd^2f/dy^2 derivative test', 'NumberTitle', 'off')
subplot(1,2,1)
contour(X, Y, reshape(fyy,Ny+1,Nx+1))
title('df/dx exact derivative')
xlabel('x')
ylabel('y')
grid off
subplot(1,2,2)
contour(X, Y, reshape(Dy*Dy*f,Ny+1,Nx+1))
title('df/dx numerical derivative')
xlabel('x')
ylabel('y')
grid off


figure('Name', 'd^2f/dy^2 derivative differences', 'NumberTitle', 'off')
surf(X, Y, reshape(fyy,Ny+1,Nx+1)-reshape(Dy*Dy*f,Ny+1,Nx+1))
title('df/dy numerical derivative error')
xlabel('x')
ylabel('y')
