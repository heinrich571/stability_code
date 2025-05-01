function j_var_R = get_var_right_inds(variable_name, Nx, Ny)

% DESCRIPTION
%   This function returns the indices of the specified variable at the right boundary of the domain in the eigenvalue problem matrices.
% INPUT:
%   variable_name       name of the variable                                                                                        [char]
% OUTPUT:
%   j_var_R             column index of specified variable at the right boundary of the domain in the eigenvalue problem matrices     [scalar]

j_var_TR = get_var_top_right_ind(variable_name, Nx, Ny);
j_var_BR = get_var_bottom_right_ind(variable_name, Nx, Ny);
j_var_R = j_var_TR : j_var_BR;

end