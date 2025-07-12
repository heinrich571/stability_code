function Left_Indices = get_operator_whole_left_indices(Nx, Ny)

% Description:
%   This function returns the row indices of the operator matrix which
%   correspond to the whole left side of the domain.
% 
% Input:
%   Nx              [int]       Number of x points.
%   Ny              [int]       Number of y points.
% 
% Output:
%   Left_Indices    [int]       Row indices of the operator at the whole
%                               left boundary of the domain.

% Input validation
if ~isnumeric(Nx)
    error('get_operator_whole_left_indices:NxNotNumeric', ...
          'Nx must be numeric.')
end
if ~isnumeric(Ny)
    error('get_operator_whole_left_indices:NyNotNumeric', ...
          'Ny must be numeric.')
end
if mod(Nx, 1) ~= 0 % Check wether Nx is an integer
    error('get_operator_whole_left_indices:NxNotInteger', ...
          ['Nx equals ' num2str(Nx) ', but it must be an integer.'])
end
if mod(Ny, 1) ~= 0 % Check wether Ny is an integer
    error('get_operator_whole_left_indices:NyNotInteger', ...
          ['Ny equals ' num2str(Ny) ', but it must be an integer.'])
end

jm = (Nx - 1) / 2 + 1; % Index at the middle of the domain along the x axis

% Check whether "jm" is an integer, throw an error in case it isn't
if jm ~= floor(jm)
    error(['Error in "get_opr_whole_left_inds": jm equals ' num2str(jm) ', but it should be a whole number.'])
end

% Calculate the indices at the whole left side of the domain
Left_Indices = (jm + 1 - 1) * Ny + 1 : Nx * Ny;

end
