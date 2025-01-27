function plotef(Solution, si, evi, Options)

% Solution  standard solution output (array of standardized structures)
% si        solution index from the solution array of structures
% evi       eigenvalue index for which plots will be generated

% Handle special inputs
if nargin < 4
    Options.Solution_Index = length(Solution);
    Options.X_Limit = Solution(Options.Solution_Index).Domain.vec_X(1);
    Options.Y_Limit = 20;
end
x_limit = Options.X_Limit;
y_limit = Options.Y_Limit;

% Get relevant information
Nx = length(Solution(si).Domain.vec_X);
Ny = length(Solution(si).Domain.vec_Y);

mat_X = Solution(si).Domain.mat_X;
mat_Y = Solution(si).Domain.mat_Y;

omega = Solution(si).Eigenvalues(evi);
u = reshape(Solution(si).Eigenfunctions.u(:,evi), Ny, Nx, length(evi));
v = reshape(Solution(si).Eigenfunctions.v(:,evi), Ny, Nx, length(evi));
w = reshape(Solution(si).Eigenfunctions.w(:,evi), Ny, Nx, length(evi));
p = reshape(Solution(si).Eigenfunctions.p(:,evi), Ny, Nx, length(evi));

ur = real(u) ; ui = imag(u) ;
vr = real(v) ; vi = imag(v) ;
wr = real(w) ; wi = imag(w) ;
pr = real(p) ; pi = imag(p) ;

% Plot eigenspectra (all of the computed eigenvalues)
figure('Name', 'Eigenspectra', 'NumberTitle', 'off')
plot(real(Solution(si).Eigenvalues), imag(Solution(si).Eigenvalues), 'o')
xlabel('$\omega_r$')
ylabel('$\omega_i$')

% Plot velocity components (for desired eigenvalues only)
for i = 1:length(evi)
    figure('Name', ['#' num2str(i) ' Velocity components'], 'NumberTitle', 'off')
    subplot(2,3,1)
    surf(mat_X, mat_Y, ur(:,:,i), 'EdgeColor', 'interp', 'FaceColor', 'none')
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$u_r$')
    xlim([-x_limit x_limit])
    ylim([0 y_limit])
    view(-50, 10)
    grid off
    subplot(2,3,2)
    surf(mat_X, mat_Y, vr(:,:,i), 'EdgeColor', 'interp', 'FaceColor', 'none')
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$v_r$')
    xlim([-x_limit x_limit])
    ylim([0 y_limit])
    view(-50, 10)
    grid off
    title(['$\beta = ' num2str(Solution(si).Physics.Beta) '$, $\omega_' num2str(evi(i)) ' = ' num2str(omega(i)) '$'] )
    subplot(2,3,3)
    surf(mat_X, mat_Y, wr(:,:,i), 'EdgeColor', 'interp', 'FaceColor', 'none')
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$w_r$')
    xlim([-x_limit x_limit])
    ylim([0 y_limit])
    view(-50, 10)
    grid off
    subplot(2,3,4)
    surf(mat_X, mat_Y, ui(:,:,i), 'EdgeColor', 'interp', 'FaceColor', 'none')
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$u_i$')
    xlim([-x_limit x_limit])
    ylim([0 y_limit])
    view(-50, 10)
    grid off
    subplot(2,3,5)
    surf(mat_X, mat_Y, vi(:,:,i), 'EdgeColor', 'interp', 'FaceColor', 'none')
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$v_i$')
    xlim([-x_limit x_limit])
    ylim([0 y_limit])
    view(-50, 10)
    grid off
    subplot(2,3,6)
    surf(mat_X, mat_Y, wi(:,:,i), 'EdgeColor', 'interp', 'FaceColor', 'none')
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$w_i$')
    xlim([-x_limit x_limit])
    ylim([0 y_limit])
    view(-50, 10)
    grid off
    
    figure('Name', ['#' num2str(i) ' Pressure'], 'NumberTitle', 'off')
    subplot(1,2,1)
    surf(mat_X, mat_Y, pr(:,:,i), 'EdgeColor', 'interp', 'FaceColor', 'none')
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$p_r$')
    xlim([-x_limit x_limit])
    ylim([0 y_limit])
    view(-50, 10)
    grid off
    title(['$\beta = ' num2str(Solution(si).Physics.Beta) '$, $\omega_' num2str(evi(i)) ' = ' num2str(omega(i)) '$'] )
    subplot(1,2,2)
    surf(mat_X, mat_Y, pi(:,:,i), 'EdgeColor', 'interp', 'FaceColor', 'none')
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$p_i$')
    xlim([-x_limit x_limit])
    ylim([0 y_limit])
    view(-50, 10)
    grid off

    figure('Name', ['#' num2str(i) ' |u|'], 'NumberTitle', 'off')
    plot(mat_X(1,:)', sqrt(ur(mat_Y(:,1) == 0,:,i).^2+ui(mat_Y(:,1) == 0,:,i).^2))
    xlabel('$x$')
    ylabel('$|u|$')
    xlim([-x_limit x_limit])
    title(['$\beta = ' num2str(Solution(si).Physics.Beta) '$, $\omega_' num2str(evi(i)) ' = ' num2str(omega(i)) '$'] )

    figure('Name', ['#' num2str(i) ' |v|'], 'NumberTitle', 'off')
    plot(mat_X(1,:)', sqrt(vr(mat_Y(:,1) == 0,:,i).^2+vi(mat_Y(:,1) == 0,:,i).^2))
    xlabel('$x$')
    ylabel('$|v|$')
    xlim([-x_limit x_limit])
    title(['$\beta = ' num2str(Solution(si).Physics.Beta) '$, $\omega_' num2str(evi(i)) ' = ' num2str(omega(i)) '$'] )

    figure('Name', ['#' num2str(i) ' |w|'], 'NumberTitle', 'off')
    plot(mat_X(1,:)', sqrt(wr(mat_Y(:,1) == 0,:,i).^2+wi(mat_Y(:,1) == 0,:,i).^2))
    xlabel('$x$')
    ylabel('$|w|$')
    xlim([-x_limit x_limit])
    title(['$\beta = ' num2str(Solution(si).Physics.Beta) '$, $\omega_' num2str(evi(i)) ' = ' num2str(omega(i)) '$'] )

    figure('Name', ['#' num2str(i) ' |p|'], 'NumberTitle', 'off')
    plot(mat_X(1,:)', sqrt(pr(mat_Y(:,1) == 0,:,i).^2+pi(mat_Y(:,1) == 0,:,i).^2))
    xlabel('$x$')
    ylabel('$|p|$')
    xlim([-x_limit x_limit])
    title(['$\beta = ' num2str(Solution(si).Physics.Beta) '$, $\omega_' num2str(evi(i)) ' = ' num2str(omega(i)) '$'] )
end


end
