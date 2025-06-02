function i_eqn_L = get_eqn_whole_left_inds(equation_name, Nx, Ny)

% DESCRIPTION
%   This function returns the indices of the specified equation at the whole left side of the domain in the eigenvalue problem matrices.
% INPUT:
%   variable_name       name of the equation                                                                                        [char]
%   Nx                  Number of points along the chordwise direction                                                              [scalar]
%   Ny                  Number of points along the vertical direction                                                               [scalar]
% OUTPUT:
%   j_var_L             column index of specified equation at the whole left side of the domain in the eigenvalue problem matrices  [scalar]

jm = (Nx - 1) / 2 + 1; % Index at the middle of the domain along the x axis

% Check whether "jm" is an integer, throw an error in case it isn't
if jm ~= floor(jm)
    error(['Error in "get_opr_whole_left_inds": jm equals ' num2str(jm) ', but it should be a whole number.'])
end

% Calculate the indices at the whole left side of the domain
i_eqn_TR = get_eqn_top_right_ind(equation_name, Nx, Ny);
i_eqn_L  = (i_eqn_TR - 1) + ((jm + 1 - 1) * Ny + 1 : Nx * Ny);

end
