%% Fresh start

clearvars -except Case_ID


%% Load results

load(['.\results\' Case_ID '.mat']);


%% Plot results

omega_ind = 1;

% ploteigfun(Domain, Solution_Filtered, 'u', omega_ind)
% ploteigfun(Domain, Solution_Filtered, 'v', omega_ind)
% ploteigfun(Domain, Solution_Filtered, 'w', omega_ind)
% ploteigfun(Domain, Solution_Filtered, 'p', omega_ind)

Nx = length(Domain.vec_X);
Ny = length(Domain.vec_Y);
% var = reshape(Solution.(varname)(:,n_omega), Ny, Nx);

% figure('Name', ['Re{' varname '}'], 'NumberTitle', 'off')
% surf(Domain.mat_X, Domain.mat_Y, real(var))
% title(['Re{' varname '} , \omega = ' num2str(real(Solution.omega(n_omega))) ' + i' num2str(imag(Solution.omega(n_omega)))])
% xlabel('x')
% ylabel('y')
% zlabel(['Re{' varname '}'], 'Interpreter', 'none')
% view(-50, 10)
% light('Position', [1,-1,1])
% light('Position', [0,0,1])
% 
% figure('Name', ['Im{' varname '}'], 'NumberTitle', 'off')
% surf(Domain.mat_X, Domain.mat_Y, imag(var))
% title(['Im{' varname '} , \omega = ' num2str(real(Solution.omega(n_omega))) ' + i' num2str(imag(Solution.omega(n_omega)))])
% xlabel('x')
% ylabel('y')
% zlabel(['Im{' varname '}'], 'Interpreter', 'none')
% view(-50, 10)
% light('Position', [1,-1,1])
% light('Position', [0,0,1])

omega = Solution_Filtered.omega(omega_ind,:);
u = reshape(Solution_Filtered.u(:,omega_ind), Ny, Nx);
v = reshape(Solution_Filtered.v(:,omega_ind), Ny, Nx);
w = reshape(Solution_Filtered.w(:,omega_ind), Ny, Nx);
p = reshape(Solution_Filtered.p(:,omega_ind), Ny, Nx);

figure('Name', 'Velocity Components', 'NumberTitle', 'off')
subplot(2,3,1)
surf(Domain.mat_X, Domain.mat_Y, real(u))
xlabel('x')
ylabel('y')
zlabel('Re{u}', 'Interpreter', 'none')
view(-50, 10)
subplot(2,3,2)
surf(Domain.mat_X, Domain.mat_Y, real(v))
xlabel('x')
ylabel('y')
zlabel('Re{v}', 'Interpreter', 'none')
view(-50, 10)
title(['\omega = ' num2str(omega)])
subplot(2,3,3)
surf(Domain.mat_X, Domain.mat_Y, real(w))
xlabel('x')
ylabel('y')
zlabel('Re{w}', 'Interpreter', 'none')
view(-50, 10)
subplot(2,3,4)
surf(Domain.mat_X, Domain.mat_Y, imag(u))
xlabel('x')
ylabel('y')
zlabel('Im{u}', 'Interpreter', 'none')
view(-50, 10)
subplot(2,3,5)
surf(Domain.mat_X, Domain.mat_Y, imag(v))
xlabel('x')
ylabel('y')
zlabel('Im{v}', 'Interpreter', 'none')
view(-50, 10)
subplot(2,3,6)
surf(Domain.mat_X, Domain.mat_Y, imag(w))
xlabel('x')
ylabel('y')
zlabel('Im{w}', 'Interpreter', 'none')
view(-50, 10)

figure('Name', 'Pressure Components', 'NumberTitle' ,'off')
subplot(1,2,1)
surf(Domain.mat_X, Domain.mat_Y, real(p))
xlabel('x')
ylabel('y')
zlabel('Re{p}', 'Interpreter', 'none')
view(-50, 10)
title(['\omega = ' num2str(omega)])
subplot(1,2,2)
surf(Domain.mat_X, Domain.mat_Y, imag(p))
xlabel('x')
ylabel('y')
zlabel('Im{p}', 'Interpreter', 'none')
view(-50, 10)




