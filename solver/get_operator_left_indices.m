function Left_Indices = get_operator_left_indices(Nx, Ny)

% Description:
%   This function returns the row indices of the operator matrix which
%   correspond to the left boundary of the domain.
% 
% Input:
%   Nx              [int]       Number of x points.
%   Ny              [int]       Number of y points.
% 
% Output:
%   Left_Indices    [int]       Row indices of the operator at the left
%                               boundary of the domain.

% Input validation
if ~isnumeric(Nx)
    error('get_operator_left_indices:NxNotNumeric', ...
          'Nx must be numeric.')
end
if ~isnumeric(Ny)
    error('get_operator_left_indices:NyNotNumeric', ...
          'Ny must be numeric.')
end
if mod(Nx, 1) ~= 0 % Check wether Nx is an integer
    error('get_operator_left_indices:NxNotInteger', ...
          ['Nx equals ' num2str(Nx) ', but it must be an integer.'])
end
if mod(Ny, 1) ~= 0 % Check wether Ny is an integer
    error('get_operator_left_indices:NyNotInteger', ...
          ['Ny equals ' num2str(Ny) ', but it must be an integer.'])
end

% Calculate the left indices of an operator
Top_Left_Index    = get_operator_top_left_index   (Nx, Ny);
Bottom_Left_Index = get_operator_bottom_left_index(Nx, Ny);

Left_Indices  = Top_Left_Index : Bottom_Left_Index;

end