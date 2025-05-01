function i_eqn_BL = get_eqn_bottom_left_ind(equation_name, Nx, Ny)

% DESCRIPTION
%   This function returns the index of the specified equation at the bottom-left corner of the domain in the eigenvalue problem matrices.
% INPUT:
%   equation_name       name of the equation                                                                                        [char]
% OUTPUT:
%   i_eqn_BL            row index of specified equation at the bottom-left corner of the domain in the eigenvalue problem matrices  [scalar]

i_eqn_TR = get_eqn_top_right_ind(equation_name, Nx, Ny);
i_eqn_BL = i_eqn_TR + Nx*Ny - 1;

end