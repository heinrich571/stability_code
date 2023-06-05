%% Fresh start

close all
clearvars
clc


%% Define the computational domain

x_limit  = 20;
y_limit  = 20;
x_meidan = 5;
y_meidan = 5;

Nx = 50;
Ny = Nx;


%% Generate the Chebyshev interval

[Dx, xhat] = cheb(Nx);
[Dy, yhat] = cheb(Ny);


%% Convert the Chebyshev interval to the physical domain

a_x_pos = x_meidan*x_limit/(x_limit-2*x_meidan);
b_x_pos = 1+2*a_x_pos/x_limit;
x_pos = a_x_pos*(1+xhat)./(b_x_pos-xhat);
Dx_pos = diag(a_x_pos*(1+b_x_pos)./((x_pos+a_x_pos).^2))*Dx;

x_meidan = -x_meidan;
x_limit  = -x_limit;

a_x_neg = x_meidan*x_limit/(x_limit-2*x_meidan);
b_x_neg = 1+2*a_x_neg/x_limit;
x_neg = a_x_neg*(1+xhat)./(b_x_neg-xhat);
Dx_neg = diag(a_x_neg*(1+b_x_neg)./((x_neg+a_x_neg).^2))*Dx;

x_meidan = -x_meidan;
x_limit  = -x_limit;

a_y = y_meidan*y_limit/(y_limit-2*y_meidan);
b_y = 1+2*a_y/y_limit;
y = a_y*(1+yhat)./(b_y-yhat);
Dy = diag(a_y*(1+b_y)./((y+a_y).^2))*Dy;

[x_pos_mesh, ~] = meshgrid(x_pos, y);
[x_neg_mesh, y_mesh] = meshgrid(x_neg, y);

Ix = eye(Nx+1);
Iy = eye(Ny+1);

DDx_pos = kron(Dx_pos, Iy);
DDx_neg = kron(Dx_neg, Iy);
DDy = kron(Ix, Dy);


%% Draw the computational domain

figure('Name', 'Computational Domain', 'NumberTitle', 'off')
plot(x_pos_mesh, y_mesh, 'k', 'LineWidth', 0.5)
hold on
plot(x_pos_mesh', y_mesh', 'k', 'LineWidth', 0.5)
plot(x_neg_mesh, y_mesh, 'k', 'LineWidth', 0.5)
plot(x_neg_mesh', y_mesh', 'k', 'LineWidth', 0.5)
title('Computational Domain')
xlabel('x')
ylabel('y')
grid off


%% Check the derivatives

% f_pos  =  sin(x_pos_mesh).*cos(y_mesh);
% fx_pos =  cos(x_pos_mesh).*cos(y_mesh);
% fy_pos = -sin(x_pos_mesh).*sin(y_mesh);
% f_pos  = f_pos(:);
% 
% f_neg  =  sin(x_neg_mesh).*cos(y_mesh);
% fx_neg =  cos(x_neg_mesh).*cos(y_mesh);
% fy_neg = -sin(x_neg_mesh).*sin(y_mesh);
% f_neg  = f_neg(:);
% 
% % fx_pos2plot = 
% 
% % d^2f/dx^2
% figure('Name', 'df/dx derivative test', 'NumberTitle', 'off')
% subplot(1,2,1)
% contour(x_pos_mesh,y_mesh,reshape(fx_pos,Ny+1,Nx+1))
% contour(x_neg_mesh,y_mesh,reshape(fx_neg,Ny+1,Nx+1))
% title('df/dx exact derivative')
% xlabel('x')
% ylabel('y')
% grid off
% subplot(1,2,2)
% contour(x_pos_mesh,y_mesh,reshape(DDx_pos*f_pos,Ny+1,Nx+1))
% contour(x_neg_mesh,y_mesh,reshape(DDx_neg*f_neg,Ny+1,Nx+1))
% title('df/dx numerical derivative')
% xlabel('x')
% ylabel('y')
% grid off
% 
% figure('Name', 'df/dx error', 'NumberTitle', 'off')
% contour(x_pos_mesh,y_mesh,reshape(DDx_pos*f_pos,Ny+1,Nx+1)-reshape(fx_pos,Ny+1,Nx+1))
% contour(x_neg_mesh,y_mesh,reshape(DDx_neg*f_neg,Ny+1,Nx+1)-reshape(fx_neg,Ny+1,Nx+1))
% title('df/dx error (numerical-exact)')
% xlabel('x')
% ylabel('y')
% grid off
% colorbar
% 
% 
% % d^2f/dy^2
% figure('Name', 'df/dy derivative test', 'NumberTitle', 'off')
% subplot(1,2,1)
% contour(x_pos_mesh,y_mesh,reshape(fy_pos,Ny+1,Nx+1))
% contour(x_neg_mesh,y_mesh,reshape(fy_neg,Ny+1,Nx+1))
% title('df/dy exact derivative')
% xlabel('x')
% ylabel('y')
% grid off
% subplot(1,2,2)
% contour(x_pos_mesh,y_mesh,reshape(DDy*f_pos,Ny+1,Nx+1))
% contour(x_neg_mesh,y_mesh,reshape(DDy*f_neg,Ny+1,Nx+1))
% title('df/dy numerical derivative')
% xlabel('x')
% ylabel('y')
% grid off
% 
% figure('Name', 'df/dy error', 'NumberTitle', 'off')
% contour(x_pos_mesh,y_mesh,reshape(DDy*f_pos,Ny+1,Nx+1)-reshape(fy_pos,Ny+1,Nx+1))
% contour(x_neg_mesh,y_mesh,reshape(DDy*f_neg,Ny+1,Nx+1)-reshape(fy_neg,Ny+1,Nx+1))
% title('df/dy error (numerical-exact)')
% xlabel('x')
% ylabel('y')
% grid off
% colorbar

