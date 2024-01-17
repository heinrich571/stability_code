%% Fresh start

close all
clear all
clc


%% Define flow parameters for the analysis

% Dimensional quantities
nu = 1;
a_vec = [1e0 1e1 1e2 1e3];

% Grid parameters
eta_Limit    = 10;
X_norm_Limit = eta_Limit;
eta_Median   = 0.2;
Nx           = 100;
Ny           = Nx/2;

% Hiemenz flow solver parameters
Definitions.initguess            = [1.22 1.24];
Definitions.maxIterations        = 1e2;
Definitions.convergenceTolerance = 1e-6;


%% Generate computational domain

Domain = generate_domain(X_norm_Limit, eta_Limit, eta_Median, Nx, Ny);


%% Compute flow fields

% Hiemenz
Definitions.interval = flip(Domain.vec_Y);
Base_Flow            = get_base_flow(Definitions);

% "Dimensionalize"
for i = 1:length(a_vec)
    a = a_vec(i);
    y{i} = flip(sqrt(nu/a)*Base_Flow.eta);
    vec_X{i} = sqrt(nu/a)*Domain.vec_X;
    vec_Y{i} = sqrt(nu/a)*Domain.vec_Y;
    mat_X{i} = sqrt(nu/a)*Domain.mat_X;
    mat_Y{i} = sqrt(nu/a)*Domain.mat_Y;
    Dx{i}    = sqrt(a/nu)*Domain.Dx;
    Dy{i}    = sqrt(a/nu)*Domain.Dy;
    D2x{i}   = (a/nu)*Domain.D2x;
    D2y{i}   = (a/nu)*Domain.D2y;
    phi{i}   = flip(Base_Flow.phi);
    dphi{i}  = flip(Base_Flow.dphi);
    ddphi{i} = flip(Base_Flow.ddphi);

    u{i} = mat_X{i}.*a.*repmat(dphi{i}, [1 length(vec_X{i})]);
    v{i} = -sqrt(a*nu)*repmat(phi{i}, [1 length(vec_X{i})]);
    
    delta_99{i} = interp1(u{i}(:,Nx/2), y{i}, 0.99*u{i}(1,Nx/2));
    Re{i} = delta_99{i}*sqrt(a/nu);

    mat_X_norm{i} = mat_X{i}/delta_99{i};
    mat_Y_norm{i} = mat_Y{i}/delta_99{i};
    u_norm{i} = u{i}/sqrt(a*nu);
    v_norm{i} = v{i}/sqrt(a*nu);

    figure('Name', ['Re = ' num2str(Re{i})])
    contourf(mat_X_norm{i}, mat_Y_norm{i}, u_norm{i}, 30)
    title('u (normalized)')
    xlabel('x/\delta_{99%}')
    ylabel('y/\delta_{99%}')

    figure('Name', ['Re = ' num2str(Re{i})])
    contourf(mat_X_norm{i}, mat_Y_norm{i}, v_norm{i}, 30)
    title('v (normalized)')
    xlabel('x/\delta_{99%}')
    ylabel('y/\delta_{99%}')
end

