function i_opr_L = get_opr_whole_left_inds(Nx, Ny)

% DESCRIPTION
%   This function returns the indices of the operator matrix on the whole left side of the domain.
% INPUT:
%   Nx                  Number of points along the chordwise direction                          [scalar]
%   Ny                  Number of points along the vertical direction                           [scalar]
% OUTPUT:
%   i_opr_L             row indices of the operator matrix at the whole left side of the domain [vector]

jm = (Nx - 1) / 2 + 1; % Index at the middle of the domain along the x axis

% Check whether "jm" is an integer, throw an error in case it isn't
if jm ~= floor(jm)
    error(['Error in "get_opr_whole_left_inds": jm equals ' num2str(jm) ', but it should be a whole number.'])
end

% Calculate the indices at the whole left side of the domain
i_opr_L = (jm + 1 - 1) * Ny + 1 : Nx * Ny;

end
