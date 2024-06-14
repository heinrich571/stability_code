function errors_plots(Solution, eigenvalues_inds)

Domain = Solution.Domain;
mat_X = Domain.mat_X;
mat_Y = Domain.mat_Y;

omega = Solution.Eigenvalues;
beta  = Solution.Physics.Beta;

default_view = [-45 45];

%% Absolute errors plot
% Absolute error maps
for i = eigenvalues_inds
    figure('Name', ['Error Map - Continuity, eigenvalue #' num2str(i)], 'NumberTitle', 'off')
    surf(mat_X, mat_Y, Solution.Report.Errors.Continuity.Absolute.Map(:,:,i))
    title(['Continuity absolute error, $\beta = ' num2str(beta) '$, $\omega = ' num2str(omega(i)) '$'])
    xlabel('$x$')
    ylabel('$y$')
    zlabel('Absolute Error')
    view(default_view)

    figure('Name', ['Error Map - X Momentum, eigenvalue #' num2str(i)], 'NumberTitle', 'off')
    surf(mat_X, mat_Y, Solution.Report.Errors.X_Momentum.Absolute.Map(:,:,i))
    title(['X-Momentum absolute error, $\beta = ' num2str(beta) '$, $\omega = ' num2str(omega(i)) '$'])
    xlabel('$x$')
    ylabel('$y$')
    view(default_view)

    figure('Name', ['Error Map - Y Momentum, eigenvalue #' num2str(i)], 'NumberTitle', 'off')
    surf(mat_X, mat_Y, Solution.Report.Errors.Y_Momentum.Absolute.Map(:,:,i))
    title(['Y-Momentum absolute error, $\beta = ' num2str(beta) '$, $\omega = ' num2str(omega(i)) '$'])
    xlabel('$x$')
    ylabel('$y$')
    view(default_view)

    figure('Name', ['Error Map - Z Momentum, eigenvalue #' num2str(i)], 'NumberTitle', 'off')
    surf(mat_X, mat_Y, Solution.Report.Errors.Z_Momentum.Absolute.Map(:,:,i))
    title(['Z-Momentum absolute error, $\beta = ' num2str(beta) '$, $\omega = ' num2str(omega(i)) '$'])
    xlabel('$x$')
    ylabel('$y$')
    view(default_view)
end

% Maximum absolute error for all eigenvalues - shown for each eigenvalue
figure('Name', 'Continuity maximum absolute errors', 'NumberTitle', 'off')
s = scatter3(real(omega)', imag(omega)', Solution.Report.Errors.Continuity.Absolute.Maximum);
title('Continuity - Maximum Absolute Error')
xlabel('$\mathrm{Re}\{\omega\}$')
ylabel('$\mathrm{Im}\{\omega\}$')
zlabel('Max. Absolute Error')
view(default_view)
s.DataTipTemplate.Interpreter = 'latex';
s.DataTipTemplate.DataTipRows(1).Label = '$\mathrm{Re}\{\omega\}$';
s.DataTipTemplate.DataTipRows(2).Label = '$\mathrm{Im}\{\omega\}$';
s.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Eigenvalue index', 1:length(omega));

figure('Name', 'X Momentum maximum absolute errors', 'NumberTitle', 'off')
s = scatter3(real(omega)', imag(omega)', Solution.Report.Errors.X_Momentum.Absolute.Maximum);
title('X Momentum - Maximum Absolute Error')
xlabel('$\mathrm{Re}\{\omega\}$')
ylabel('$\mathrm{Im}\{\omega\}$')
zlabel('Max. Absolute Error')
view(default_view)
s.DataTipTemplate.Interpreter = 'latex';
s.DataTipTemplate.DataTipRows(1).Label = '$\mathrm{Re}\{\omega\}$';
s.DataTipTemplate.DataTipRows(2).Label = '$\mathrm{Im}\{\omega\}$';
s.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Eigenvalue index', 1:length(omega));

figure('Name', 'Y Momentum maximum absolute errors', 'NumberTitle', 'off')
s = scatter3(real(omega)', imag(omega)', Solution.Report.Errors.Y_Momentum.Absolute.Maximum);
title('Y Momentum - Maximum Absolute Error')
xlabel('$\mathrm{Re}\{\omega\}$')
ylabel('$\mathrm{Im}\{\omega\}$')
zlabel('Max. Absolute Error')
view(default_view)
s.DataTipTemplate.Interpreter = 'latex';
s.DataTipTemplate.DataTipRows(1).Label = '$\mathrm{Re}\{\omega\}$';
s.DataTipTemplate.DataTipRows(2).Label = '$\mathrm{Im}\{\omega\}$';
s.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Eigenvalue index', 1:length(omega));

figure('Name', 'Z Momentum maximum absolute errors', 'NumberTitle', 'off')
s = scatter3(real(omega)', imag(omega)', Solution.Report.Errors.Z_Momentum.Absolute.Maximum);
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
figure('Name', 'Continuity maximum absolute errors', 'NumberTitle', 'off')
plot(Domain.mat_X, Domain.mat_Y, 'k', 'LineWidth', 1)
plot(Domain.mat_X', Domain.mat_Y', 'k', 'LineWidth', 1)
s = scatter3(Solution.Report.Errors.Continuity.Absolute.Maximum_X, Solution.Report.Errors.Continuity.Absolute.Maximum_Y, Solution.Report.Errors.Continuity.Absolute.Maximum);
title('Continuity - Maximum Absolute Error')
xlabel('$x$')
ylabel('$y$')
zlabel('Max. Absolute Error')
view(default_view)
s.DataTipTemplate.Interpreter = 'latex';
s.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Eigenvalue index', 1:length(omega));

figure('Name', 'X Momentum maximum absolute errors', 'NumberTitle', 'off')
plot(Domain.mat_X, Domain.mat_Y, 'k', 'LineWidth', 1)
plot(Domain.mat_X', Domain.mat_Y', 'k', 'LineWidth', 1)
s = scatter3(Solution.Report.Errors.Continuity.Absolute.Maximum_X, Solution.Report.Errors.Continuity.Absolute.Maximum_Y, Solution.Report.Errors.X_Momentum.Absolute.Maximum);
title('X Momentum - Maximum Absolute Error')
xlabel('$x$')
ylabel('$y$')
zlabel('Max. Absolute Error')
view(default_view)
s.DataTipTemplate.Interpreter = 'latex';
s.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Eigenvalue index', 1:length(omega));

figure('Name', 'Y Momentum maximum absolute errors', 'NumberTitle', 'off')
plot(Domain.mat_X, Domain.mat_Y, 'k', 'LineWidth', 1)
plot(Domain.mat_X', Domain.mat_Y', 'k', 'LineWidth', 1)
s = scatter3(Solution.Report.Errors.Continuity.Absolute.Maximum_X, Solution.Report.Errors.Continuity.Absolute.Maximum_Y, Solution.Report.Errors.Y_Momentum.Absolute.Maximum);
title('Y Momentum - Maximum Absolute Error')
xlabel('$x$')
ylabel('$y$')
zlabel('Max. Absolute Error')
view(default_view)
s.DataTipTemplate.Interpreter = 'latex';
s.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Eigenvalue index', 1:length(omega));

figure('Name', 'Z Momentum maximum absolute errors', 'NumberTitle', 'off')
plot(Domain.mat_X, Domain.mat_Y, 'k', 'LineWidth', 1)
plot(Domain.mat_X', Domain.mat_Y', 'k', 'LineWidth', 1)
s = scatter3(Solution.Report.Errors.Continuity.Absolute.Maximum_X, Solution.Report.Errors.Continuity.Absolute.Maximum_Y, Solution.Report.Errors.Z_Momentum.Absolute.Maximum);
title('Z Momentum - Maximum Absolute Error')
xlabel('$x$')
ylabel('$y$')
zlabel('Max. Absolute Error')
view(default_view)
s.DataTipTemplate.Interpreter = 'latex';
s.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Eigenvalue index', 1:length(omega));

end
