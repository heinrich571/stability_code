function errors_plots(Solution, Report, eigenvalues_inds)

if isempty(fieldnames(Report))
    return;
end

Domain = Solution.Domain;
mat_X = Domain.mat_X;
mat_Y = Domain.mat_Y;

omega = Solution.Eigenvalues;
beta  = Solution.Physics.Beta;

default_view = [-45 45];

%% Absolute errors plot - Navier Stokes Check
% Absolute error maps
for i = eigenvalues_inds
    figure('Name', ['LNSE Error Map - Continuity, eigenvalue #' num2str(i)], 'NumberTitle', 'off')
    surf(mat_X, mat_Y, Report.Navier_Stokes_Check.Errors.Continuity.Absolute.Map(:,:,i))
    title(['Continuity absolute error, $\beta = ' num2str(beta) '$, $\omega_' num2str(i) ' = ' num2str(omega(i)) '$'])
    xlabel('$x$')
    ylabel('$y$')
    zlabel('Absolute Error')
    view(default_view)

    figure('Name', ['LNSE Error Map - X Momentum, eigenvalue #' num2str(i)], 'NumberTitle', 'off')
    surf(mat_X, mat_Y, Report.Navier_Stokes_Check.Errors.X_Momentum.Absolute.Map(:,:,i))
    title(['X-Momentum absolute error, $\beta = ' num2str(beta) '$, $\omega_' num2str(i) ' = ' num2str(omega(i)) '$'])
    xlabel('$x$')
    ylabel('$y$')
    view(default_view)

    figure('Name', ['LNSE Error Map - Y Momentum, eigenvalue #' num2str(i)], 'NumberTitle', 'off')
    surf(mat_X, mat_Y, Report.Navier_Stokes_Check.Errors.Y_Momentum.Absolute.Map(:,:,i))
    title(['Y-Momentum absolute error, $\beta = ' num2str(beta) '$, $\omega_' num2str(i) ' = ' num2str(omega(i)) '$'])
    xlabel('$x$')
    ylabel('$y$')
    view(default_view)

    figure('Name', ['LNSE Error Map - Z Momentum, eigenvalue #' num2str(i)], 'NumberTitle', 'off')
    surf(mat_X, mat_Y, Report.Navier_Stokes_Check.Errors.Z_Momentum.Absolute.Map(:,:,i))
    title(['Z-Momentum absolute error, $\beta = ' num2str(beta) '$, $\omega_' num2str(i) ' = ' num2str(omega(i)) '$'])
    xlabel('$x$')
    ylabel('$y$')
    view(default_view)
end

% Maximum absolute error for all eigenvalues - shown for each eigenvalue
figure('Name', 'LNSE Continuity maximum absolute errors', 'NumberTitle', 'off')
s = scatter3(real(omega)', imag(omega)', Report.Navier_Stokes_Check.Errors.Continuity.Absolute.Maximum);
title('Continuity - Maximum Absolute Error')
xlabel('$\mathrm{Re}\{\omega\}$')
ylabel('$\mathrm{Im}\{\omega\}$')
zlabel('Max. Absolute Error')
view(default_view)
s.DataTipTemplate.Interpreter = 'latex';
s.DataTipTemplate.DataTipRows(1).Label = '$\mathrm{Re}\{\omega\}$';
s.DataTipTemplate.DataTipRows(2).Label = '$\mathrm{Im}\{\omega\}$';
s.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Eigenvalue index', 1:length(omega));

figure('Name', 'LNSE X Momentum maximum absolute errors', 'NumberTitle', 'off')
s = scatter3(real(omega)', imag(omega)', Report.Navier_Stokes_Check.Errors.X_Momentum.Absolute.Maximum);
title('X Momentum - Maximum Absolute Error')
xlabel('$\mathrm{Re}\{\omega\}$')
ylabel('$\mathrm{Im}\{\omega\}$')
zlabel('Max. Absolute Error')
view(default_view)
s.DataTipTemplate.Interpreter = 'latex';
s.DataTipTemplate.DataTipRows(1).Label = '$\mathrm{Re}\{\omega\}$';
s.DataTipTemplate.DataTipRows(2).Label = '$\mathrm{Im}\{\omega\}$';
s.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Eigenvalue index', 1:length(omega));

figure('Name', 'LNSE Y Momentum maximum absolute errors', 'NumberTitle', 'off')
s = scatter3(real(omega)', imag(omega)', Report.Navier_Stokes_Check.Errors.Y_Momentum.Absolute.Maximum);
title('Y Momentum - Maximum Absolute Error')
xlabel('$\mathrm{Re}\{\omega\}$')
ylabel('$\mathrm{Im}\{\omega\}$')
zlabel('Max. Absolute Error')
view(default_view)
s.DataTipTemplate.Interpreter = 'latex';
s.DataTipTemplate.DataTipRows(1).Label = '$\mathrm{Re}\{\omega\}$';
s.DataTipTemplate.DataTipRows(2).Label = '$\mathrm{Im}\{\omega\}$';
s.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Eigenvalue index', 1:length(omega));

figure('Name', 'LNSE Z Momentum maximum absolute errors', 'NumberTitle', 'off')
s = scatter3(real(omega)', imag(omega)', Report.Navier_Stokes_Check.Errors.Z_Momentum.Absolute.Maximum);
title('Z Momentum - Maximum Absolute Error')
xlabel('$\mathrm{Re}\{\omega\}$')
ylabel('$\mathrm{Im}\{\omega\}$')
zlabel('Max. Absolute Error')
view(default_view)
s.DataTipTemplate.Interpreter = 'latex';
s.DataTipTemplate.DataTipRows(1).Label = '$\mathrm{Re}\{\omega\}$';
s.DataTipTemplate.DataTipRows(2).Label = '$\mathrm{Im}\{\omega\}$';
s.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Eigenvalue index', 1:length(omega));


% Maximum absolute error for all eigenvalues - shown for each eigenvalue on
% the domain of solution
figure('Name', 'LNSE Continuity maximum absolute errors', 'NumberTitle', 'off')
plot(Domain.mat_X, Domain.mat_Y, 'k', 'LineWidth', 1)
plot(Domain.mat_X', Domain.mat_Y', 'k', 'LineWidth', 1)
s = scatter3(Report.Navier_Stokes_Check.Errors.Continuity.Absolute.Maximum_X, Report.Navier_Stokes_Check.Errors.Continuity.Absolute.Maximum_Y, Report.Navier_Stokes_Check.Errors.Continuity.Absolute.Maximum);
title('Continuity - Maximum Absolute Error')
xlabel('$x$')
ylabel('$y$')
zlabel('Max. Absolute Error')
view(default_view)
s.DataTipTemplate.Interpreter = 'latex';
s.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Eigenvalue index', 1:length(omega));

figure('Name', 'LNSE X Momentum maximum absolute errors', 'NumberTitle', 'off')
plot(Domain.mat_X, Domain.mat_Y, 'k', 'LineWidth', 1)
plot(Domain.mat_X', Domain.mat_Y', 'k', 'LineWidth', 1)
s = scatter3(Report.Navier_Stokes_Check.Errors.Continuity.Absolute.Maximum_X, Report.Navier_Stokes_Check.Errors.Continuity.Absolute.Maximum_Y, Report.Navier_Stokes_Check.Errors.X_Momentum.Absolute.Maximum);
title('X Momentum - Maximum Absolute Error')
xlabel('$x$')
ylabel('$y$')
zlabel('Max. Absolute Error')
view(default_view)
s.DataTipTemplate.Interpreter = 'latex';
s.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Eigenvalue index', 1:length(omega));

figure('Name', 'LNSE Y Momentum maximum absolute errors', 'NumberTitle', 'off')
plot(Domain.mat_X, Domain.mat_Y, 'k', 'LineWidth', 1)
plot(Domain.mat_X', Domain.mat_Y', 'k', 'LineWidth', 1)
s = scatter3(Report.Navier_Stokes_Check.Errors.Continuity.Absolute.Maximum_X, Report.Navier_Stokes_Check.Errors.Continuity.Absolute.Maximum_Y, Report.Navier_Stokes_Check.Errors.Y_Momentum.Absolute.Maximum);
title('Y Momentum - Maximum Absolute Error')
xlabel('$x$')
ylabel('$y$')
zlabel('Max. Absolute Error')
view(default_view)
s.DataTipTemplate.Interpreter = 'latex';
s.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Eigenvalue index', 1:length(omega));

figure('Name', 'LNSE Z Momentum maximum absolute errors', 'NumberTitle', 'off')
plot(Domain.mat_X, Domain.mat_Y, 'k', 'LineWidth', 1)
plot(Domain.mat_X', Domain.mat_Y', 'k', 'LineWidth', 1)
s = scatter3(Report.Navier_Stokes_Check.Errors.Continuity.Absolute.Maximum_X, Report.Navier_Stokes_Check.Errors.Continuity.Absolute.Maximum_Y, Report.Navier_Stokes_Check.Errors.Z_Momentum.Absolute.Maximum);
title('Z Momentum - Maximum Absolute Error')
xlabel('$x$')
ylabel('$y$')
zlabel('Max. Absolute Error')
view(default_view)
s.DataTipTemplate.Interpreter = 'latex';
s.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Eigenvalue index', 1:length(omega));


%% Absolute errors plot - EVP Check
% Absolute error maps
for i = eigenvalues_inds
    figure('Name', ['EVP Error Map - Continuity, eigenvalue #' num2str(i)], 'NumberTitle', 'off')
    surf(mat_X, mat_Y, Report.EVP_Check.Errors.Continuity.Absolute.Map(:,:,i))
    title(['Continuity absolute error, $\beta = ' num2str(beta) '$, $\omega = ' num2str(omega(i)) '$'])
    xlabel('$x$')
    ylabel('$y$')
    zlabel('Absolute Error')
    view(default_view)

    figure('Name', ['EVP Error Map - X Momentum, eigenvalue #' num2str(i)], 'NumberTitle', 'off')
    surf(mat_X, mat_Y, Report.EVP_Check.Errors.X_Momentum.Absolute.Map(:,:,i))
    title(['X-Momentum absolute error, $\beta = ' num2str(beta) '$, $\omega = ' num2str(omega(i)) '$'])
    xlabel('$x$')
    ylabel('$y$')
    view(default_view)

    figure('Name', ['EVP Error Map - Y Momentum, eigenvalue #' num2str(i)], 'NumberTitle', 'off')
    surf(mat_X, mat_Y, Report.EVP_Check.Errors.Y_Momentum.Absolute.Map(:,:,i))
    title(['Y-Momentum absolute error, $\beta = ' num2str(beta) '$, $\omega = ' num2str(omega(i)) '$'])
    xlabel('$x$')
    ylabel('$y$')
    view(default_view)

    figure('Name', ['EVP Error Map - Z Momentum, eigenvalue #' num2str(i)], 'NumberTitle', 'off')
    surf(mat_X, mat_Y, Report.EVP_Check.Errors.Z_Momentum.Absolute.Map(:,:,i))
    title(['Z-Momentum absolute error, $\beta = ' num2str(beta) '$, $\omega = ' num2str(omega(i)) '$'])
    xlabel('$x$')
    ylabel('$y$')
    view(default_view)
end

% Maximum absolute error for all eigenvalues - shown for each eigenvalue
figure('Name', 'EVP Continuity maximum absolute errors', 'NumberTitle', 'off')
s = scatter3(real(omega)', imag(omega)', Report.EVP_Check.Errors.Continuity.Absolute.Maximum);
title('Continuity - Maximum Absolute Error')
xlabel('$\mathrm{Re}\{\omega\}$')
ylabel('$\mathrm{Im}\{\omega\}$')
zlabel('Max. Absolute Error')
view(default_view)
s.DataTipTemplate.Interpreter = 'latex';
s.DataTipTemplate.DataTipRows(1).Label = '$\mathrm{Re}\{\omega\}$';
s.DataTipTemplate.DataTipRows(2).Label = '$\mathrm{Im}\{\omega\}$';
s.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Eigenvalue index', 1:length(omega));

figure('Name', 'EVP X Momentum maximum absolute errors', 'NumberTitle', 'off')
s = scatter3(real(omega)', imag(omega)', Report.EVP_Check.Errors.X_Momentum.Absolute.Maximum);
title('X Momentum - Maximum Absolute Error')
xlabel('$\mathrm{Re}\{\omega\}$')
ylabel('$\mathrm{Im}\{\omega\}$')
zlabel('Max. Absolute Error')
view(default_view)
s.DataTipTemplate.Interpreter = 'latex';
s.DataTipTemplate.DataTipRows(1).Label = '$\mathrm{Re}\{\omega\}$';
s.DataTipTemplate.DataTipRows(2).Label = '$\mathrm{Im}\{\omega\}$';
s.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Eigenvalue index', 1:length(omega));

figure('Name', 'EVP Y Momentum maximum absolute errors', 'NumberTitle', 'off')
s = scatter3(real(omega)', imag(omega)', Report.EVP_Check.Errors.Y_Momentum.Absolute.Maximum);
title('Y Momentum - Maximum Absolute Error')
xlabel('$\mathrm{Re}\{\omega\}$')
ylabel('$\mathrm{Im}\{\omega\}$')
zlabel('Max. Absolute Error')
view(default_view)
s.DataTipTemplate.Interpreter = 'latex';
s.DataTipTemplate.DataTipRows(1).Label = '$\mathrm{Re}\{\omega\}$';
s.DataTipTemplate.DataTipRows(2).Label = '$\mathrm{Im}\{\omega\}$';
s.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Eigenvalue index', 1:length(omega));

figure('Name', 'EVP Z Momentum maximum absolute errors', 'NumberTitle', 'off')
s = scatter3(real(omega)', imag(omega)', Report.EVP_Check.Errors.Z_Momentum.Absolute.Maximum);
title('Z Momentum - Maximum Absolute Error')
xlabel('$\mathrm{Re}\{\omega\}$')
ylabel('$\mathrm{Im}\{\omega\}$')
zlabel('Max. Absolute Error')
view(default_view)
s.DataTipTemplate.Interpreter = 'latex';
s.DataTipTemplate.DataTipRows(1).Label = '$\mathrm{Re}\{\omega\}$';
s.DataTipTemplate.DataTipRows(2).Label = '$\mathrm{Im}\{\omega\}$';
s.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Eigenvalue index', 1:length(omega));


% Maximum absolute error for all eigenvalues - shown for each eigenvalue on
% the domain of solution
figure('Name', 'EVP Continuity maximum absolute errors', 'NumberTitle', 'off')
plot(Domain.mat_X, Domain.mat_Y, 'k', 'LineWidth', 1)
plot(Domain.mat_X', Domain.mat_Y', 'k', 'LineWidth', 1)
s = scatter3(Report.EVP_Check.Errors.Continuity.Absolute.Maximum_X, Report.EVP_Check.Errors.Continuity.Absolute.Maximum_Y, Report.EVP_Check.Errors.Continuity.Absolute.Maximum);
title('Continuity - Maximum Absolute Error')
xlabel('$x$')
ylabel('$y$')
zlabel('Max. Absolute Error')
view(default_view)
s.DataTipTemplate.Interpreter = 'latex';
s.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Eigenvalue index', 1:length(omega));

figure('Name', 'EVP X Momentum maximum absolute errors', 'NumberTitle', 'off')
plot(Domain.mat_X, Domain.mat_Y, 'k', 'LineWidth', 1)
plot(Domain.mat_X', Domain.mat_Y', 'k', 'LineWidth', 1)
s = scatter3(Report.EVP_Check.Errors.Continuity.Absolute.Maximum_X, Report.EVP_Check.Errors.Continuity.Absolute.Maximum_Y, Report.EVP_Check.Errors.X_Momentum.Absolute.Maximum);
title('X Momentum - Maximum Absolute Error')
xlabel('$x$')
ylabel('$y$')
zlabel('Max. Absolute Error')
view(default_view)
s.DataTipTemplate.Interpreter = 'latex';
s.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Eigenvalue index', 1:length(omega));

figure('Name', 'EVP Y Momentum maximum absolute errors', 'NumberTitle', 'off')
plot(Domain.mat_X, Domain.mat_Y, 'k', 'LineWidth', 1)
plot(Domain.mat_X', Domain.mat_Y', 'k', 'LineWidth', 1)
s = scatter3(Report.EVP_Check.Errors.Continuity.Absolute.Maximum_X, Report.EVP_Check.Errors.Continuity.Absolute.Maximum_Y, Report.EVP_Check.Errors.Y_Momentum.Absolute.Maximum);
title('Y Momentum - Maximum Absolute Error')
xlabel('$x$')
ylabel('$y$')
zlabel('Max. Absolute Error')
view(default_view)
s.DataTipTemplate.Interpreter = 'latex';
s.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Eigenvalue index', 1:length(omega));

figure('Name', 'EVP Z Momentum maximum absolute errors', 'NumberTitle', 'off')
plot(Domain.mat_X, Domain.mat_Y, 'k', 'LineWidth', 1)
plot(Domain.mat_X', Domain.mat_Y', 'k', 'LineWidth', 1)
s = scatter3(Report.EVP_Check.Errors.Continuity.Absolute.Maximum_X, Report.EVP_Check.Errors.Continuity.Absolute.Maximum_Y, Report.EVP_Check.Errors.Z_Momentum.Absolute.Maximum);
title('Z Momentum - Maximum Absolute Error')
xlabel('$x$')
ylabel('$y$')
zlabel('Max. Absolute Error')
view(default_view)
s.DataTipTemplate.Interpreter = 'latex';
s.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Eigenvalue index', 1:length(omega));

end
