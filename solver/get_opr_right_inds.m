function i_opr_R = get_opr_right_inds(Nx, Ny)

% DESCRIPTION
%   This function returns the indices of the operator matrix at the right boundary of the domain.
% INPUT:
%   Nx                  Number of points along the chordwise direction                          [scalar]
%   Ny                  Number of points along the vertical direction                           [scalar]
% OUTPUT:
%   i_opr_R             row indices of the operator matrix at the right boundary of the domain  [vector]

i_opr_TR = get_opr_top_right_ind(Nx, Ny);
i_opr_BR = get_opr_bottom_right_ind(Nx, Ny);
i_opr_R  = i_opr_TR : i_opr_BR;

end