function i_eqn_BR = get_eqn_bottom_right_ind(equation_name, Nx, Ny)

% DESCRIPTION
%   This function returns the index of the specified equation at the bottom-right corner of the domain in the eigenvalue problem matrices.
% INPUT:
%   equation_name       name of the equation                                                                                            [char]
% OUTPUT:
%   i_eqn_BR            row index of specified equation at the bottom-right corner of the domain in the eigenvalue problem matrices     [scalar]

i_eqn_TR = get_eqn_top_right_ind(equation_name, Nx, Ny);
i_eqn_BR = i_eqn_TR + (Ny-1);

end