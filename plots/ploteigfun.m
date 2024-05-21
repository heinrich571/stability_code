function ploteigfun(Domain, Solution, varname, n_omega)

Nx = length(Domain.vec_X);
Ny = length(Domain.vec_Y);
var = reshape(Solution.(varname)(:,n_omega), Ny, Nx);

figure('Name', ['Re{' varname '}'], 'NumberTitle', 'off')
surf(Domain.mat_X, Domain.mat_Y, real(var))
title(['Re{' varname '} , \omega = ' num2str(real(Solution.omega(n_omega))) ' + i' num2str(imag(Solution.omega(n_omega)))])
xlabel('x')
ylabel('y')
zlabel(['Re{' varname '}'], 'Interpreter', 'none')
view(-50, 10)
light('Position', [1,-1,1])
light('Position', [0,0,1])

figure('Name', ['Im{' varname '}'], 'NumberTitle', 'off')
surf(Domain.mat_X, Domain.mat_Y, imag(var))
title(['Im{' varname '} , \omega = ' num2str(real(Solution.omega(n_omega))) ' + i' num2str(imag(Solution.omega(n_omega)))])
xlabel('x')
ylabel('y')
zlabel(['Im{' varname '}'], 'Interpreter', 'none')
view(-50, 10)
light('Position', [1,-1,1])
light('Position', [0,0,1])

end
