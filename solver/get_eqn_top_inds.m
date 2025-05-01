function i_eqn_T = get_eqn_top_inds(equation_name, Nx, Ny)

% DESCRIPTION
%   This function returns the indices of the specified equation at the top boundary of the domain in the eigenvalue problem matrices.
% INPUT:
%   equation_name       name of the equation                                                                                    [char]
% OUTPUT:
%   i_eqn_T             row index of specified equation at the top boundary of the domain in the eigenvalue problem matrices    [scalar]

i_eqn_TR = get_eqn_top_right_ind(equation_name, Nx, Ny);
i_eqn_TL = get_eqn_top_left_ind(equation_name, Nx, Ny);
i_eqn_T  = i_eqn_TR : Ny : i_eqn_TL;

end