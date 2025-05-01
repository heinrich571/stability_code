function j_var_B = get_var_bottom_inds(variable_name, Nx, Ny)

% DESCRIPTION
%   This function returns the indices of the specified variable at the bottom boundary of the domain in the eigenvalue problem matrices.
% INPUT:
%   variable_name       name of the variable                                                                                        [char]
% OUTPUT:
%   j_var_B             column index of specified variable at the bottom boundary of the domain in the eigenvalue problem matrices  [scalar]

j_var_BR = get_var_bottom_right_ind(variable_name, Nx, Ny);
j_var_BL = get_var_bottom_left_ind(variable_name, Nx, Ny);
j_var_B  = j_var_BR : Ny : j_var_BL;

end