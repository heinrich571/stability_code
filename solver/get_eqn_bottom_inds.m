function i_eqn_B = get_eqn_bottom_inds(equation_name, Nx, Ny)

% DESCRIPTION
%   This function returns the indices of the specified equation at the bottom boundary of the domain in the eigenvalue problem matrices.
% INPUT:
%   equation_name       name of the equation                                                                                        [char]
% OUTPUT:
%   i_eqn_B             row index of specified equation at the bottom boundary of the domain in the eigenvalue problem matrices     [scalar]

i_eqn_BR = get_eqn_bottom_right_ind(equation_name, Nx, Ny);
i_eqn_BL = get_eqn_bottom_left_ind(equation_name, Nx, Ny);
i_eqn_B  = i_eqn_BR : Ny : i_eqn_BL;

end