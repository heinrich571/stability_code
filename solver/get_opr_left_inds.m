function i_opr_L = get_opr_left_inds(Nx, Ny)

% DESCRIPTION
%   This function returns the indices of the operator matrix at the left boundary of the domain.
% INPUT:
%   Nx                  Number of points along the chordwise direction                          [scalar]
%   Ny                  Number of points along the vertical direction                           [scalar]
% OUTPUT:
%   i_opr_L             row indices of the operator matrix at the left boundary of the domain   [vector]

i_opr_TL = get_opr_top_left_ind(Nx, Ny);
i_opr_BL = get_opr_bottom_left_ind(Nx, Ny);
i_opr_L  = i_opr_TL : i_opr_BL;

end