function j_var_BR = get_var_bottom_right_ind(variable_name, Nx, Ny)

% DESCRIPTION
%   This function returns the index of the specified variable at the bottom-right corner of the domain in the eigenvalue problem matrices.
% INPUT:
%   variable_name       name of the variable                                                                                            [char]
% OUTPUT:
%   j_var_BR            column index of specified variable at the bottom-right corner of the domain in the eigenvalue problem matrices  [scalar]

j_var_TR = get_var_top_right_ind(variable_name, Nx, Ny);
j_var_BR = j_var_TR + (Ny-1);

end