function j_var_T = get_var_top_inds(variable_name, Nx, Ny)

% DESCRIPTION
%   This function returns the indices of the specified variable at the top boundary of the domain in the eigenvalue problem matrices.
% INPUT:
%   variable_name       name of the variable                                                                                        [char]
% OUTPUT:
%   j_var_T             column index of specified variable at the top boundary of the domain in the eigenvalue problem matrices     [scalar]

j_var_TR = get_var_top_right_ind(variable_name, Nx, Ny);
j_var_TL = get_var_top_left_ind(variable_name, Nx, Ny);
j_var_T  = j_var_TR : Ny : j_var_TL;

end