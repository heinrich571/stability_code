clc
clear all
close all


Ny = 11;
Ly = 5;
y = linspace(-Ly,Ly,Ny+1)';
hy = 2*Ly/Ny;

% 4th order
Dy = 1/hy*(1/12*diag(ones(Ny-1,1),-2) -2/3*diag(ones(Ny,1),-1) + ...
           2/3*diag(ones(Ny,1),1)   -1/12*diag(ones(Ny-1,1),2));
Dy(1,1:5) = 1/hy*[-25/12 4 -3 4/3 -1/4]; % at the wall
Dy(2,1:5) = 1/hy*[-3/12 -10/12 18/12 -6/12 1/12]; % one point from the wall
Dy(Ny+1,Ny+1-(0:4)) = -1/hy*[-25/12 4 -3 4/3 -1/4]; % at the wall
Dy(Ny,Ny+1-(0:4)) = -1/hy*[-3/12 -10/12 18/12 -6/12 1/12]; % one point from the wall

Nx = 10;
Lx = 4;
x = linspace(-Lx,Lx,Nx+1)';
hx = 2*Lx/Nx;

% 4th order
Dx = 1/hx*(1/12*diag(ones(Nx-1,1),-2) -2/3*diag(ones(Nx,1),-1) + ...
           2/3*diag(ones(Nx,1),1)   -1/12*diag(ones(Nx-1,1),2));
Dx(1,1:5) = 1/hx*[-25/12 4 -3 4/3 -1/4]; % at the wall
Dx(2,1:5) = 1/hx*[-3/12 -10/12 18/12 -6/12 1/12]; % one point from the wall
Dx(Nx+1,Nx+1-(0:4)) = -1/hx*[-25/12 4 -3 4/3 -1/4]; % at the wall
Dx(Nx,Nx+1-(0:4)) = -1/hx*[-3/12 -10/12 18/12 -6/12 1/12]; % one point from the wall

% make in 2d:
[x,y]=meshgrid(x,y);


f = sin(x).*cos(y);
fx = cos(x).*cos(y);
fy = sin(x).*-sin(y);

f = f(:);

Iy = eye(Ny+1);
Ix = eye(Nx+1);
DDx = kron(Dx,Iy);
DDy = kron(Ix,Dy);

figure
subplot(1,2,1)
contour(x,y,reshape(fx,Ny+1,Nx+1))
subplot(1,2,2)
contour(x,y,reshape(DDx*f,Ny+1,Nx+1))


figure
subplot(1,2,1)
contour(x,y,reshape(fy,Ny+1,Nx+1))
subplot(1,2,2)
contour(x,y,reshape(DDy*f,Ny+1,Nx+1))



% % for f as a matrix:
% figure
% subplot(1,2,1)
% contour(x,y,fy)
% subplot(1,2,2)
% contour(x,y,Dy*f)
% 
% 
% figure
% subplot(1,2,1)
% contour(x,y,fx)
% subplot(1,2,2)
% contour(x,y,(Dx*f')')


