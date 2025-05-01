function i_eqn_TR = get_eqn_top_right_ind(equation_name, Nx, Ny)

% DESCRIPTION
%   This function returns the index of the specified equation at the top-right corner of the domain in the eigenvalue problem matrices.
% INPUT:
%   equation_name       name of the equation                                                                                        [char]
% OUTPUT:
%   i_eqn_TR            row index of specified equation at the top-right corner of the domain in the eigenvalue problem matrices    [scalar]

equation_number = eqn2num(equation_name);
i_eqn_TR = 1 + Nx*Ny*(equation_number - 1);

end