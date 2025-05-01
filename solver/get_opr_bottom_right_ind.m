function i_opr_BR = get_opr_bottom_right_ind(Nx, Ny)

% DESCRIPTION
%   This function returns the index of the operator matrix at the bottom-right corner of the domain.
% INPUT:
%   Nx                  Number of points along the chordwise direction                              [scalar]
%   Ny                  Number of points along the vertical direction                               [scalar]
% OUTPUT:
%   i_opr_BR            row index of the operator matrix at the bottom-right corner of the domain   [vector]

i_opr_TR = get_opr_top_right_ind(Nx, Ny);
i_opr_BR = i_opr_TR + (Ny-1);

end