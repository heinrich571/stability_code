function i_opr_B = get_opr_bottom_inds(Nx, Ny)

% DESCRIPTION
%   This function returns the indices of the operator matrix at the bottom boundary of the domain.
% INPUT:
%   Nx                  Number of points along the chordwise direction                              [scalar]
%   Ny                  Number of points along the vertical direction                               [scalar]
% OUTPUT:
%   i_opr_B             row indices of the operator matrix at the bottom boundary of the domain     [vector]

i_opr_BR = get_opr_bottom_right_ind(Nx, Ny);
i_opr_BL = get_opr_bottom_left_ind(Nx, Ny);
i_opr_B  = i_opr_BR : Ny : i_opr_BL;

end