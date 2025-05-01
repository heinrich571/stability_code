function i_eqn_L = get_eqn_left_inds(equation_name, Nx, Ny)

% DESCRIPTION
%   This function returns the indices of the specified variable at the left boundary of the domain in the eigenvalue problem matrices.
% INPUT:
%   equation_name       name of the equation                                                                                    [char]
% OUTPUT:
%   i_eqn_L             row index of specified equation at the left boundary of the domain in the eigenvalue problem matrices   [scalar]

i_eqn_TL = get_eqn_top_left_ind(equation_name, Nx, Ny);
i_eqn_BL = get_eqn_bottom_left_ind(equation_name, Nx, Ny);
i_eqn_L  = i_eqn_TL : i_eqn_BL;

end