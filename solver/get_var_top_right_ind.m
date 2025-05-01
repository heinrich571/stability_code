function j_var_TR = get_var_top_right_ind(variable_name, Nx, Ny)

% DESCRIPTION
%   This function returns the index of the specified variable at the top-right corner of the domain in the eigenvalue problem matrices.
% INPUT:
%   variable_name       name of the variable                                                                                        [char]
% OUTPUT:
%   j_var_TR            column index of specified variable at the top-right corner of the domain in the eigenvalue problem matrices [scalar]

variable_number = var2num(variable_name);
j_var_TR = 1 + Nx*Ny*(variable_number - 1);

end