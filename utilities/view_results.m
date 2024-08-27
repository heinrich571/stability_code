function view_results(Results_Filename, Results_Folder, Eigenvalue_Indices, Options)

% Load results from file
if strcmp(Results_Filename(end-3:end), '.mat')
    Results_Path = [Results_Folder '\' Results_Filename];
else
    Results_Path = [Results_Folder '\' Results_Filename '.mat'];
end
load(Results_Path, 'Problem', 'Solution')

% Handle special inputs
if nargin < 4
    Options.Solution_Index = length(Solution)-1;
    Options.X_Limit = Solution(Options.Solution_Index).Domain.vec_X(1);
    Options.Y_Limit = 20;
end
x_limit = Options.X_Limit;
y_limit = Options.Y_Limit;


% Get relevant information
si = Options.Solution_Index;
Nx = length(Solution(si).Domain.vec_X);
Ny = length(Solution(si).Domain.vec_Y);

mat_X = Solution(si).Domain.mat_X;
mat_Y = Solution(si).Domain.mat_Y;

omega = Solution(si).Eigenvalues(Eigenvalue_Indices);
u = reshape(Solution(si).Eigenfunctions.u(:,Eigenvalue_Indices), Ny, Nx, length(Eigenvalue_Indices));
v = reshape(Solution(si).Eigenfunctions.v(:,Eigenvalue_Indices), Ny, Nx, length(Eigenvalue_Indices));
w = reshape(Solution(si).Eigenfunctions.w(:,Eigenvalue_Indices), Ny, Nx, length(Eigenvalue_Indices));
p = reshape(Solution(si).Eigenfunctions.p(:,Eigenvalue_Indices), Ny, Nx, length(Eigenvalue_Indices));

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
for i = 1:length(Eigenvalue_Indices)
    figure('Name', ['#' num2str(i) ' Velocity components'], 'NumberTitle', 'off')
    subplot(2,3,1)
    surf(mat_X, mat_Y, ur(:,:,i))
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$u_r$')
    xlim([-x_limit x_limit])
    ylim([0 y_limit])
    view(-50, 10)
    subplot(2,3,2)
    surf(mat_X, mat_Y, vr(:,:,i))
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$v_r$')
    xlim([-x_limit x_limit])
    ylim([0 y_limit])
    view(-50, 10)
    title(['$\beta = ' num2str(Problem(si).Physics.Beta) '$, $\omega_' num2str(i) ' = ' num2str(omega(i)) '$'] )
    subplot(2,3,3)
    surf(mat_X, mat_Y, wr(:,:,i))
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$w_r$')
    xlim([-x_limit x_limit])
    ylim([0 y_limit])
    view(-50, 10)
    subplot(2,3,4)
    surf(mat_X, mat_Y, ui(:,:,i))
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$u_i$')
    xlim([-x_limit x_limit])
    ylim([0 y_limit])
    view(-50, 10)
    subplot(2,3,5)
    surf(mat_X, mat_Y, vi(:,:,i))
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$v_i$')
    xlim([-x_limit x_limit])
    ylim([0 y_limit])
    view(-50, 10)
    subplot(2,3,6)
    surf(mat_X, mat_Y, wi(:,:,i))
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$w_i$')
    xlim([-x_limit x_limit])
    ylim([0 y_limit])
    view(-50, 10)
    
    figure('Name', ['#' num2str(i) ' Pressure'], 'NumberTitle', 'off')
    subplot(1,2,1)
    surf(mat_X, mat_Y, pr(:,:,i))
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$p_r$')
    xlim([-x_limit x_limit])
    ylim([0 y_limit])
    view(-50, 10)
    title(['$\beta = ' num2str(Problem(si).Physics.Beta) '$, $\omega_' num2str(i) ' = ' num2str(omega(i)) '$'] )
    subplot(1,2,2)
    surf(mat_X, mat_Y, pi(:,:,i))
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$p_i$')
    xlim([-x_limit x_limit])
    ylim([0 y_limit])
    view(-50, 10)
end


end
