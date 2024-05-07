function Solution = biglobalsolver(mat_A, mat_B, Domain)

Nx = length(Domain.vec_X);
Ny = length(Domain.vec_Y);

% [eigenfunctions_matrix, eigenvalues_matrix] = eig(mat_A, mat_B);
[eigenfunctions_matrix, eigenvalues_matrix] = eigs(sparse(mat_A), sparse(mat_B));

Solution.omega = diag(eigenvalues_matrix);

Solution.u = get_eigenfunction_of(eigenfunctions_matrix, 'u', Nx, Ny);
Solution.v = get_eigenfunction_of(eigenfunctions_matrix, 'v', Nx, Ny);
Solution.w = get_eigenfunction_of(eigenfunctions_matrix, 'w', Nx, Ny);
Solution.p = get_eigenfunction_of(eigenfunctions_matrix, 'p', Nx, Ny);

end

% Supporting functions
function var_eigfun = get_eigenfunction_of(eigenfunctions_matrix, varname, Nx, Ny)

varnames = {'u' , 'v' , 'w' , 'p'};
varorder = find(strcmp(varnames, varname), 1, 'first');
var_inds = (1:(Nx*Ny)) + (varorder-1)*Nx*Ny;

var_eigfun = eigenfunctions_matrix(var_inds,:);

end