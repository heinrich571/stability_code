function [Solution_Filtered, Solution_Raw] = biglobalsolver(mat_A, mat_B, Domain)

Nx = length(Domain.vec_X);
Ny = length(Domain.vec_Y);

% [eigenfunctions_matrix, eigenvalues_matrix] = eig(mat_A, mat_B);
[eigenfunctions_matrix, eigenvalues_matrix] = eigs(sparse(mat_A), sparse(mat_B), 20, 'SM');

Solution_Raw.omega = diag(eigenvalues_matrix);

Solution_Raw.u = get_eigenfunction_of(eigenfunctions_matrix, 'u', Nx, Ny);
Solution_Raw.v = get_eigenfunction_of(eigenfunctions_matrix, 'v', Nx, Ny);
Solution_Raw.w = get_eigenfunction_of(eigenfunctions_matrix, 'w', Nx, Ny);
Solution_Raw.p = get_eigenfunction_of(eigenfunctions_matrix, 'p', Nx, Ny);

omega_magnitude_threshold = 100;
inds = find(sqrt(real(Solution_Raw.omega).^2+imag(Solution_Raw.omega).^2) <= omega_magnitude_threshold);
Solution_Filtered.omega = Solution_Raw.omega(inds);
Solution_Filtered.u = Solution_Raw.u(:,inds);
Solution_Filtered.v = Solution_Raw.v(:,inds);
Solution_Filtered.w = Solution_Raw.w(:,inds);
Solution_Filtered.p = Solution_Raw.p(:,inds);

end

% Supporting functions
function var_eigfun = get_eigenfunction_of(eigenfunctions_matrix, varname, Nx, Ny)

varnames = {'u' , 'v' , 'w' , 'p'};
varorder = find(strcmp(varnames, varname), 1, 'first');
var_inds = (1:(Nx*Ny)) + (varorder-1)*Nx*Ny;

var_eigfun = eigenfunctions_matrix(var_inds,:);

end