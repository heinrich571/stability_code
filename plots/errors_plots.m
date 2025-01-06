function errors_plots(Solution, Report, eigenvalues_inds)

if isempty(fieldnames(Report))
    return;
end

Domain = Solution.Domain;
Nx = length(Domain.vec_X);
Ny = length(Domain.vec_Y);
mat_X = Domain.mat_X;
mat_Y = Domain.mat_Y;

omega = Solution.Eigenvalues;
beta  = Solution.Physics.Beta;

default_view = [-45 45];


%% Matrix A (LHS) equations
A = Report.Navier_Stokes_Check.A;
B = Report.Navier_Stokes_Check.B;
inds.cont = 3*Nx*Ny+1:4*Nx*Ny;
inds.xmom = 0*Nx*Ny+1:1*Nx*Ny;
inds.ymom = 1*Nx*Ny+1:2*Nx*Ny;
inds.zmom = 2*Nx*Ny+1:3*Nx*Ny;
for i = eigenvalues_inds
    q = [Solution.Eigenfunctions.u(:,i) ; Solution.Eigenfunctions.v(:,i) ; Solution.Eigenfunctions.w(:,i) ; Solution.Eigenfunctions.p(:,i)];
    w = Solution.Eigenvalues(i);

    LHS = A*q;
    RHS = w*B*q;
    rLHS = real(LHS);
    iLHS = imag(LHS);
    rRHS = real(RHS);
    iRHS = imag(RHS);
    
    figure('Name', ['continuity, eigenvalue #' num2str(i)], 'NumberTitle', 'off')
    subplot(2,2,1)
    surf(mat_X, mat_Y, reshape(rLHS(inds.cont), Ny, Nx))
    title(['LNSE - Continuity Re(LHS), $\beta = ' num2str(beta) '$, $\omega_' num2str(i) ' = ' num2str(omega(i)) '$'])
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$\mathrm{Re}\{LHS\}$')
    view(default_view)

    subplot(2,2,2)
    surf(mat_X, mat_Y, reshape(rRHS(inds.cont), Ny, Nx))
    title(['LNSE - Continuity Re(RHS), $\beta = ' num2str(beta) '$, $\omega_' num2str(i) ' = ' num2str(omega(i)) '$'])
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$\mathrm{Re}\{RHS\}$')
    view(default_view)

    subplot(2,2,3)
    surf(mat_X, mat_Y, reshape(iLHS(inds.cont), Ny, Nx))
    title(['LNSE - Continuity Im(LHS), $\beta = ' num2str(beta) '$, $\omega_' num2str(i) ' = ' num2str(omega(i)) '$'])
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$\mathrm{Im}\{LHS\}$')
    view(default_view)

    subplot(2,2,4)
    surf(mat_X, mat_Y, reshape(iRHS(inds.cont), Ny, Nx))
    title(['LNSE - Continuity Im(RHS), $\beta = ' num2str(beta) '$, $\omega_' num2str(i) ' = ' num2str(omega(i)) '$'])
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$\mathrm{Im}\{RHS\}$')
    view(default_view)

    figure('Name', ['x-momentum, eigenvalue #' num2str(i)], 'NumberTitle', 'off')
    subplot(2,2,1)
    surf(mat_X, mat_Y, reshape(rLHS(inds.xmom), Ny, Nx))
    title(['LNSE - X-Momentum Re(LHS), $\beta = ' num2str(beta) '$, $\omega_' num2str(i) ' = ' num2str(omega(i)) '$'])
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$\mathrm{Re}\{LHS\}$')
    view(default_view)

    subplot(2,2,2)
    surf(mat_X, mat_Y, reshape(rRHS(inds.xmom), Ny, Nx))
    title(['LNSE - X-Momentum Re(RHS), $\beta = ' num2str(beta) '$, $\omega_' num2str(i) ' = ' num2str(omega(i)) '$'])
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$\mathrm{Re}\{RHS\}$')
    view(default_view)

    subplot(2,2,3)
    surf(mat_X, mat_Y, reshape(iLHS(inds.xmom), Ny, Nx))
    title(['LNSE - X-Momentum Im(LHS), $\beta = ' num2str(beta) '$, $\omega_' num2str(i) ' = ' num2str(omega(i)) '$'])
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$\mathrm{Im}\{LHS\}$')
    view(default_view)

    subplot(2,2,4)
    surf(mat_X, mat_Y, reshape(iRHS(inds.xmom), Ny, Nx))
    title(['LNSE - X-Momentum Im(RHS), $\beta = ' num2str(beta) '$, $\omega_' num2str(i) ' = ' num2str(omega(i)) '$'])
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$\mathrm{Im}\{RHS\}$')
    view(default_view)

    figure('Name', ['y-momentum, eigenvalue #' num2str(i)], 'NumberTitle', 'off')
    subplot(2,2,1)
    surf(mat_X, mat_Y, reshape(rLHS(inds.ymom), Ny, Nx))
    title(['LNSE - Y-Momentum Re(LHS), $\beta = ' num2str(beta) '$, $\omega_' num2str(i) ' = ' num2str(omega(i)) '$'])
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$\mathrm{Re}\{LHS\}$')
    view(default_view)

    subplot(2,2,2)
    surf(mat_X, mat_Y, reshape(rRHS(inds.ymom), Ny, Nx))
    title(['LNSE - Y-Momentum Re(RHS), $\beta = ' num2str(beta) '$, $\omega_' num2str(i) ' = ' num2str(omega(i)) '$'])
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$\mathrm{Re}\{RHS\}$')
    view(default_view)

    subplot(2,2,3)
    surf(mat_X, mat_Y, reshape(iLHS(inds.ymom), Ny, Nx))
    title(['LNSE - Y-Momentum Im(LHS), $\beta = ' num2str(beta) '$, $\omega_' num2str(i) ' = ' num2str(omega(i)) '$'])
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$\mathrm{Im}\{LHS\}$')
    view(default_view)

    subplot(2,2,4)
    surf(mat_X, mat_Y, reshape(iRHS(inds.ymom), Ny, Nx))
    title(['LNSE - Y-Momentum Im(RHS), $\beta = ' num2str(beta) '$, $\omega_' num2str(i) ' = ' num2str(omega(i)) '$'])
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$\mathrm{Im}\{RHS\}$')
    view(default_view)

    figure('Name', ['z-momentum, eigenvalue #' num2str(i)], 'NumberTitle', 'off')
    subplot(2,2,1)
    surf(mat_X, mat_Y, reshape(rLHS(inds.zmom), Ny, Nx))
    title(['LNSE - Z-Momentum Re(LHS), $\beta = ' num2str(beta) '$, $\omega_' num2str(i) ' = ' num2str(omega(i)) '$'])
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$\mathrm{Re}\{LHS\}$')
    view(default_view)

    subplot(2,2,2)
    surf(mat_X, mat_Y, reshape(rRHS(inds.zmom), Ny, Nx))
    title(['LNSE - Z-Momentum Re(RHS), $\beta = ' num2str(beta) '$, $\omega_' num2str(i) ' = ' num2str(omega(i)) '$'])
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$\mathrm{Re}\{RHS\}$')
    view(default_view)

    subplot(2,2,3)
    surf(mat_X, mat_Y, reshape(iLHS(inds.zmom), Ny, Nx))
    title(['LNSE - Z-Momentum Im(LHS), $\beta = ' num2str(beta) '$, $\omega_' num2str(i) ' = ' num2str(omega(i)) '$'])
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$\mathrm{Im}\{LHS\}$')
    view(default_view)

    subplot(2,2,4)
    surf(mat_X, mat_Y, reshape(iRHS(inds.zmom), Ny, Nx))
    title(['LNSE - Z-Momentum Im(RHS), $\beta = ' num2str(beta) '$, $\omega_' num2str(i) ' = ' num2str(omega(i)) '$'])
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$\mathrm{Im}\{RHS\}$')
    view(default_view)
end


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
