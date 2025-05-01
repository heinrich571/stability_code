function j_var_TL = get_var_top_left_ind(variable_name, Nx, Ny)

% DESCRIPTION
%   This function returns the index of the specified variable at the top-left corner of the domain in the eigenvalue problem matrices.
% INPUT:
%   variable_name       name of the variable                                                                                            [char]
% OUTPUT:
%   j_var_TL            column index of specified variable at the top-left corner of the domain in the eigenvalue problem matrices      [scalar]

j_var_BL = get_var_bottom_left_ind(variable_name, Nx, Ny);
j_var_TL = j_var_BL - (Ny-1);

end