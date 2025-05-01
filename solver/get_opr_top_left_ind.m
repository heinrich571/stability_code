function i_opr_TL = get_opr_top_left_ind(Nx, Ny)

% DESCRIPTION
%   This function returns the index of the operator matrix at the top-left corner of the domain.
% INPUT:
%   Nx                  Number of points along the chordwise direction                          [scalar]
%   Ny                  Number of points along the vertical direction                           [scalar]
% OUTPUT:
%   i_opr_TL            row index of the operator matrix at the top-left corner of the domain   [vector]

i_opr_BL = get_opr_bottom_left_ind(Nx, Ny);
i_opr_TL = i_opr_BL - (Ny-1);

end