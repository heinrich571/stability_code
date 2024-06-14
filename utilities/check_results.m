%% Fresh start

clearvars -except Case_ID Results_Folder


%% Load results

load([Results_Folder Case_ID '.mat']);


%% Plot results

omega_ind = 10;

Nx = length(Solution.Domain.vec_X);
Ny = length(Solution.Domain.vec_Y);

omega = Solution.Eigenvalues(omega_ind);
u = reshape(Solution.Eigenfunctions.u(:,omega_ind), Ny, Nx);
v = reshape(Solution.Eigenfunctions.v(:,omega_ind), Ny, Nx);
w = reshape(Solution.Eigenfunctions.w(:,omega_ind), Ny, Nx);
p = reshape(Solution.Eigenfunctions.p(:,omega_ind), Ny, Nx);

Re_u = real(u) ; Im_u = imag(u) ;
Re_v = real(v) ; Im_v = imag(v) ;
Re_w = real(w) ; Im_w = imag(w) ;
Re_p = real(p) ; Im_p = imag(p) ;


figure('Name', 'Eigenspectra', 'NumberTitle', 'off')
plot(real(Solution.Eigenvalues), imag(Solution.Eigenvalues), 'o')
xlabel('$\mathrm{Re}\{\omega\}$')
ylabel('$\mathrm{Im}\{\omega\}$')

figure('Name', 'Velocity Components', 'NumberTitle', 'off')
subplot(2,3,1)
surf(Solution.Domain.mat_X, Solution.Domain.mat_Y, Re_u)
xlabel('$x$')
ylabel('$y$')
zlabel('$\mathrm{Re}\{u\}$')
view(-50, 10)
subplot(2,3,2)
surf(Solution.Domain.mat_X, Solution.Domain.mat_Y, Re_v)
xlabel('$x$')
ylabel('$y$')
zlabel('$\mathrm{Re}\{v\}$')
view(-50, 10)
title(['$\beta = ' num2str(Problem.Physics.Beta) '$, $\omega = ' num2str(omega) '$'] )
subplot(2,3,3)
surf(Solution.Domain.mat_X, Solution.Domain.mat_Y, Re_w)
xlabel('$x$')
ylabel('$y$')
zlabel('$\mathrm{Re}\{w\}$')
view(-50, 10)
subplot(2,3,4)
surf(Solution.Domain.mat_X, Solution.Domain.mat_Y, Im_u)
xlabel('$x$')
ylabel('$y$')
zlabel('$\mathrm{Im}\{u\}$')
view(-50, 10)
subplot(2,3,5)
surf(Solution.Domain.mat_X, Solution.Domain.mat_Y, Im_v)
xlabel('$x$')
ylabel('$y$')
zlabel('$\mathrm{Im}\{v\}$')
view(-50, 10)
subplot(2,3,6)
surf(Solution.Domain.mat_X, Solution.Domain.mat_Y, Im_w)
xlabel('$x$')
ylabel('$y$')
zlabel('$\mathrm{Im}\{w\}$')
view(-50, 10)

figure('Name', 'Pressure Components', 'NumberTitle' ,'off')
subplot(1,2,1)
surf(Solution.Domain.mat_X, Solution.Domain.mat_Y, Re_p)
xlabel('$x$')
ylabel('$y$')
zlabel('$\mathrm{Re}\{p\}$')
view(-50, 10)
title(['$\beta = ' num2str(Problem.Physics.Beta) '$, $\omega = ' num2str(omega) '$'] )
subplot(1,2,2)
surf(Solution.Domain.mat_X, Solution.Domain.mat_Y, Im_p)
xlabel('$x$')
ylabel('$y$')
zlabel('$\mathrm{Im}\{p\}$')
view(-50, 10)



