%% Fresh Start

close all
clearvars -except Problem Solution
clc
startup

path_manager('add')


%% Input

% Results file settings
Results_Folder = '.\batch_results\Domain_Resolution_Sensitivity_Test\';
Results_File   = 'Domain_Resolution_Sensitivity_Test_2nd_zero_deriv.mat';
% Results_File   = 'Domain_Resolution_Sensitivity_Test_linextrap.mat';

% Eigenvalues to plot
% Eigenvalues_To_Plot = [1 8];
Eigenvalues_To_Plot = [1 9];
% Eigenvalues_To_Plot = [1 8];


%% Load Results

if ~exist('Solution', 'var')
    load([Results_Folder '\' Results_File]);
end


%% Calculate Graph Variables

% Find the number of eigenvalues
N_eigenvalues = length(Solution(1).Eigenvalues);

% Find the number of nodes in the x & y directions
Nx_vec = 0;
Ny_vec = 0;
for i = 1:length(Solution)
    Nx = length(Solution(i).Domain.vec_X);
    Ny = length(Solution(i).Domain.vec_Y);
    if ~sum(Nx_vec == Nx)
        Nx_vec(end+1) = Nx;
    end
    if ~sum(Ny_vec == Ny)
        Ny_vec(end+1) = Ny;
    end
end
Nx_vec = Nx_vec(2:end) - 1;
Ny_vec = Ny_vec(2:end) - 1;


%% Plots

% Eigenspectra
symbol_list = {'o' 's' 'd' '^' 'x' '*'};
figure('Name', 'Eigenspectra', 'NumberTitle', 'off')
k = 1;
for i = 1:length(Nx_vec)
    figure('Name', ['Eigenspectra for Nx = ' num2str(Nx_vec(i))], 'NumberTitle', 'off')
    for j = 1:length(Ny_vec)
        current_symbol_ind = mod(j, length(symbol_list))+1;
        plot(real(Solution(k).Eigenvalues), imag(Solution(k).Eigenvalues), ...
            symbol_list{current_symbol_ind}, 'DisplayName', ['N_y = ' num2str(Ny_vec(j))], ...
            'LineWidth', 0.5, 'MarkerSize', 10)
        title(['Eigenspectra for $N_x = ' num2str(Nx_vec(i)) '$'])
        xlabel('$\omega_r$')
        ylabel('$\omega_i$')
        xlim([-5 5])
        ylim([-4 0])
        grid off

        k = k + 1;
    end
end

% time vs. resolution
if exist('Monitor', 'var')
    figure('Name', 'Solution time vs. domain resolution', 'NumberTitle', 'off')
    plot(Nx, t, '-o')
    xlabel('$N$')
    ylabel('$t$')
    title('Solution time vs. domain resolution')
    grid minor
end


% Convergence of eigenvalues
for i = 1:length(Nx_vec)
    Nx = Nx_vec(i);
    i_ref_sol = length(Ny_vec)*i;
    for j = Eigenvalues_To_Plot
        count = length(Ny_vec)*(i-1)+1;
        for n = 1:length(Ny_vec)
            domega = Solution(count).Eigenvalues - Solution(i_ref_sol).Eigenvalues(j);
            d = sqrt(real(domega).^2 + imag(domega).^2);
            [dmin, i_dmin] = min(d);
            omega(n,1) = Solution(count).Eigenvalues(i_dmin);
            count = count + 1;
        end
        figure('Name', ['Eigenvalue #' num2str(j) ' convergence, Nx = ' num2str(Nx)], 'NumberTitle', 'off')
        subplot(1,2,1)
        plot(Ny_vec, real(omega'), '-o', 'LineWidth', 0.5)
        xlabel('$N_y$')
        ylabel(['$\mathrm{Re}\{\omega_' num2str(j) '\}$'])
        subplot(1,2,2)
        plot(Ny_vec, imag(omega'), '-o', 'LineWidth', 0.5)
        xlabel('$N_y$')
        ylabel(['$\mathrm{Im}\{\omega_' num2str(j) '\}$'])

        figure('Name', ['Eigenvalue #' num2str(j) ' convergence, Nx = ' num2str(Nx)], 'NumberTitle', 'off')
        title(['$\omega = ' num2str(real(omega(end))) ' + ' num2str(imag(omega(end))) 'i$'])
        plot(Ny_vec, log10(abs(omega-omega(end)))', '-^', 'LineWidth', 0.5)
        xlabel('$N$')
        ylabel(['$\log\left(\varepsilon_{\omega_' num2str(j) '}\right)$'])
    end
end
