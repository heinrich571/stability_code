function j_var_BL = get_var_bottom_left_ind(variable_name, Nx, Ny)

% DESCRIPTION
%   This function returns the index of the specified variable at the bottom-left corner of the domain in the eigenvalue problem matrices.
% INPUT:
%   variable_name       name of the variable                                                                                            [char]
% OUTPUT:
%   i_eqn_BL            column index of specified variable at the bottom-left corner of the domain in the eigenvalue problem matrices   [scalar]

j_var_TR = get_var_top_right_ind(variable_name, Nx, Ny);
j_var_BL = j_var_TR + Nx*Ny - 1;

end