%% Fresh start

close all
clearvars
clc


%% Define scalar parameters

x_limit  = 2;
y_limit = 2;
y_median = 0.5;
d = 0.5;
c = 1-d;


%% Define vector parameters

Nx_vec = 10:10:200;


%% Compute the differentiation matrices

initializer = NaN([1 numel(Nx_vec)]);
max_err_dx = initializer;
max_err_dy = initializer;
for i = 1:numel(Nx_vec)
    Nx = Nx_vec(i);
    Ny = Nx/2;
    Ix = eye(Nx+1);
    Iy = eye(Ny+1);
    
    % Compute differentiation matrices
    [Dx_cheb, xhat] = cheb(Nx);
    [Dy_cheb, yhat] = cheb(Ny);
    
    % Calculate the wall-normal differentiation matrix
    a_y = y_median*y_limit/(y_limit-2*y_median);
    b_y = 1+2*a_y/y_limit;
    y = (a_y*(1+yhat)./(b_y-yhat));
    Dy_physical_domain = diag(a_y*(1+b_y)./((y+a_y).^2))*Dy_cheb;
    Dy = kron(Ix, Dy_physical_domain);
    
    % Calculate the chordwise differentiation matrix
    x = c*xhat.^3 + d*xhat;
    Dx_physical_domain = diag(1./(3*c*xhat.^2+d))*Dx_cheb;
    Dx = kron(Dx_physical_domain, Iy);

    % Calculate the exact derivative of the test function
    [X, Y] = meshgrid(x, y);
%     f        = X.*Y;
%     fx_exact = Y;
%     fy_exact = X;
    f        = sin(X).*cos(Y);
    fx_exact =  cos(X).*cos(Y);
    fy_exact = -sin(X).*sin(Y);
    
    % Calculate the approximated derivative of the test function
    fx_numerical = reshape(Dx*f(:), Ny+1, Nx+1);
    fy_numerical = reshape(Dy*f(:), Ny+1, Nx+1);

    % Calculate the maximum error of the approximation
    max_err_dx(i) = max(abs(fx_numerical - fx_exact), [], 'all');
    max_err_dy(i) = max(abs(fy_numerical - fy_exact), [], 'all');

end


%% Plot results

figure('Name', 'Error of first derivative', 'NumberTitle', 'off')
subplot(2,1,1)
plot(Nx_vec, max_err_dx, '-o')
xlabel('$N_x$', 'Interpreter', 'LaTeX')
ylabel('$\epsilon_{\frac{\partial f}{\partial x}}$', 'Interpreter', 'LaTeX')
subplot(2,1,2)
plot(Nx_vec, max_err_dy, '-s')
xlabel('$N_x$', 'Interpreter', 'LaTeX')
ylabel('$\epsilon_{\frac{\partial f}{\partial y}}$', 'Interpreter', 'LaTeX')