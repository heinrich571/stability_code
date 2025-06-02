function j_var_L = get_var_whole_left_inds(variable_name, Nx, Ny)

% DESCRIPTION
%   This function returns the indices of the specified variable at the whole left side of the domain in the eigenvalue problem matrices.
% INPUT:
%   variable_name       name of the variable                                                                                        [char]
%   Nx                  Number of points along the chordwise direction                                                              [scalar]
%   Ny                  Number of points along the vertical direction                                                               [scalar]
% OUTPUT:
%   j_var_L             column index of specified variable at the whole left side of the domain in the eigenvalue problem matrices  [scalar]

jm = (Nx - 1) / 2 + 1; % Index at the middle of the domain along the x axis

% Check whether "jm" is an integer, throw an error in case it isn't
if jm ~= floor(jm)
    error(['Error in "get_opr_whole_left_inds": jm equals ' num2str(jm) ', but it should be a whole number.'])
end

% Calculate the indices at the whole left side of the domain
j_var_TR = get_var_top_right_ind(variable_name, Nx, Ny);
j_var_L  = (j_var_TR - 1) + ((jm + 1 - 1) * Ny + 1 : Nx * Ny);

end
