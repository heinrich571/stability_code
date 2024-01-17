%% Fresh start

close all
clear all
clc

addpath(genpath('../../grid_generation/'))


%% Define flow parameters

a  = 1;
nu = 14.61e-6;
rho = 1.225;

% Grid parameters
eta_Limit    = 20;
X_norm_Limit = eta_Limit;
eta_Median   = 2.4;
Nx           = 100;
Ny           = Nx/2;

% Hiemenz flow solver parameters
Definitions.initguess            = [1.22 1.24];
Definitions.maxIterations        = 1e2;
Definitions.convergenceTolerance = 1e-6;


%% Generate the computational domain

Domain = generate_domain(X_norm_Limit, eta_Limit, eta_Median, Nx, Ny);


%% Solve the Hiemenz problem

Definitions.interval = flip(Domain.vec_Y);
% Definitions.interval = [0 eta_Limit];
Base_Flow            = get_base_flow(Definitions);
Base_Flow.y          = sqrt(nu/a)*Base_Flow.eta;


%% Create the dimensional velocity field

Domain.vec_X = sqrt(nu/a)*Domain.vec_X;
Domain.vec_Y = sqrt(nu/a)*Domain.vec_Y;
Domain.mat_X = sqrt(nu/a)*Domain.mat_X;
Domain.mat_Y = sqrt(nu/a)*Domain.mat_Y;
Domain.Dx    = sqrt(a/nu)*Domain.Dx;
Domain.Dy    = sqrt(a/nu)*Domain.Dy;
Domain.D2x   = (a/nu)*Domain.D2x;
Domain.D2y   = (a/nu)*Domain.D2y;

Base_Flow.y     = flip(Base_Flow.y);
Base_Flow.phi   = flip(Base_Flow.phi);
Base_Flow.dphi  = flip(Base_Flow.dphi);
Base_Flow.ddphi = flip(Base_Flow.ddphi);

u = Domain.mat_X.*a.*repmat(Base_Flow.dphi, [1 length(Domain.vec_X)]);
v = -sqrt(a*nu)*repmat(Base_Flow.phi, [1 length(Domain.vec_X)]);

du_dx   = Domain.Dx* u(:);
du_dy   = Domain.Dy* u(:);
d2u_dx2 = Domain.D2x*u(:);
d2u_dy2 = Domain.D2y*u(:);

dv_dx   = Domain.Dx* v(:);
dv_dy   = Domain.Dy* v(:);
d2v_dx2 = Domain.D2x*v(:);
d2v_dy2 = Domain.D2y*v(:);


%% Check continuity

continuity_err = du_dx + dv_dy;
continuity_err = reshape(abs(continuity_err), [length(Domain.vec_Y) length(Domain.vec_X)]);

figure('Name', 'Continuity equation check', 'NumberTitle', 'off')
contourf(Domain.mat_X, Domain.mat_Y, continuity_err)
title('continuity')
xlabel('x')
ylabel('y')


%% Check x-momentum

dp_dx = -rho*a^2*Domain.mat_X;
x_momentum_err = u(:).*du_dx + v(:).*du_dy + 1/rho*dp_dx(:) - nu*(d2u_dx2 + d2u_dy2); % RHS - LHS
x_momentum_err = reshape(abs(x_momentum_err), [Ny+1 Nx+1]);

figure('Name', 'x-momentum equation check', 'NumberTitle', 'off')
contourf(Domain.mat_X, Domain.mat_Y, x_momentum_err)
title('x-momentum')
xlabel('x')
ylabel('y')


%% Check y-momentum

dp_dy = -rho*sqrt(a^3*nu)*(Base_Flow.phi.*Base_Flow.dphi + Base_Flow.ddphi);
dp_dy = repmat(dp_dy, [1 Nx+1]);
y_momentum_err = u(:).*dv_dx + v(:).*dv_dy + 1/rho*dp_dy(:) - nu*(d2v_dx2 + d2v_dy2); % RHS - LHS
% y_momentum_err = - nu*(d2v_dx2 + d2v_dy2); % RHS - LHS
y_momentum_err = reshape(abs(y_momentum_err), [Ny+1 Nx+1]);

figure('Name', 'y-momentum equation check', 'NumberTitle', 'off')
contourf(Domain.mat_X, Domain.mat_Y, y_momentum_err)
title('y-momentum')
xlabel('x')
ylabel('y')














