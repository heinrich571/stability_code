function i_eqn_R = get_eqn_right_inds(equation_name, Nx, Ny)

% DESCRIPTION
%   This function returns the indices of the specified variable at the right boundary of the domain in the eigenvalue problem matrices.
% INPUT:
%   equation_name       name of the equation                                                                                        [char]
% OUTPUT:
%   i_eqn_R             row index of specified equation at the right boundary of the domain in the eigenvalue problem matrices      [scalar]

i_eqn_TR = get_eqn_top_right_ind(equation_name, Nx, Ny);
i_eqn_BR = get_eqn_bottom_right_ind(equation_name, Nx, Ny);
i_eqn_R  = i_eqn_TR : i_eqn_BR;

end