function Bottom_Left_Index = get_operator_bottom_left_index(Nx, Ny)

% Description:
%   This function returns the row index of the operator matrix that
%   corresponds to the bottom-left corner of the domain.
% 
% Input:
%   Nx                  [int]       Number of x points.
%   Ny                  [int]       Number of y points.
% 
% Output:
%   Bottom_Left_Index   [int]       Row index of the operator at the
%                                   bottom-left corner of the domain.

% Input validation
if ~isnumeric(Nx)
    error('get_operator_bottom_left_index:NxNotNumeric', ...
          'Nx must be numeric.')
end
if ~isnumeric(Ny)
    error('get_operator_bottom_left_index:NyNotNumeric', ...
          'Ny must be numeric.')
end
if mod(Nx, 1) ~= 0 % Check wether Nx is an integer
    error('get_operator_bottom_left_index:NxNotInteger', ...
          ['Nx equals ' num2str(Nx) ', but it must be an integer.'])
end
if mod(Ny, 1) ~= 0 % Check wether Ny is an integer
    error('get_operator_bottom_left_index:NyNotInteger', ...
          ['Ny equals ' num2str(Ny) ', but it must be an integer.'])
end

% Calculate the bottom-left index of an operator
Bottom_Left_Index = Nx * Ny;

end
