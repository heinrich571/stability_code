%% Fresh Start

close all
clear all
clc

addpath('..\..\baseflow_definition\')
addpath('..\..\grid_generation\')


%% Define Parameters

ymedian  = 3; % value of median y value

ny_vec   = 10:10:150; % number of wall normal points
ymax_vec = 20:20:300; % values of maximum practical farfield distance from wall


%% Base Flow Calculation

Problem.Domain.Nx       = 20;
Problem.Domain.X_Limit  = 100;
Problem.Domain.Y_Median = ymedian;

Base_Flow_Definitions.initguess            = 1.23258765682022 + [-1 1]*1e-5;
Base_Flow_Definitions.maxIterations        = 1e2;
Base_Flow_Definitions.convergenceTolerance = 1e-6;

for i = 1:length(ny_vec)
    ny = ny_vec(i);

    for j = 1:length(ymax_vec)
        ymax = ymax_vec(j);

        Problem.Domain.Ny      = ny;
        Problem.Domain.Y_Limit = ymax;
        Domain = generate_domain(Problem);

        Base_Flow_Definitions.interval = Domain.vec_Y;
        Base_Flow(i,j) = get_base_flow(Base_Flow_Definitions);
    end
end


%% Process Errors

err_max = zeros(size(Base_Flow));
% dphi_ref = Base_Flow(end,end).dphi(abs(Base_Flow(end,end).eta - ymedian) < 1e-6);
for i = 1:length(ny_vec)
    for j = 1:length(ymax_vec)
        dphi_ref = interp1(Base_Flow(end,end).eta, Base_Flow(end,end).dphi, ymedian, 'pchip');
        err_max(i,j) = max(abs(Base_Flow(i,j).dphi(abs(Base_Flow(i,j).eta - ymedian) < 1e-6) - dphi_ref));
    end
end


%% Plots

figure('Name', 'BL Profile')
for i = 1:length(ny_vec)
    ny = ny_vec(i);
    plot(Base_Flow(i,end).dphi, Base_Flow(i,end).eta, '-o', 'DisplayName', ['$N_y = ' num2str(ny) '$'])
end
xlabel('$\phi''$')
ylabel('$y$')
legend('interpreter', 'latex', 'location', 'northwest')
ylim([0 10])

[NY, YMAX] = meshgrid(ny_vec, ymax_vec);

figure('Name', 'log10(error) surface')
surf(NY, YMAX, err_max')
xlabel('$N_y$')
ylabel('$y_{max}$')
zlabel('$\varepsilon(y_{median})$')


%% Cleanup Paths

rmpath('..\..\baseflow_definition\')
rmpath('..\..\grid_generation\')
