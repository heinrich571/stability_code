function Right_Indices = get_operator_right_indices(Nx, Ny)

% Description:
%   This function returns the row indices of the operator matrix which
%   correspond to the right boundary of the domain.
% 
% Input:
%   Nx               [int]      Number of x points.
%   Ny               [int]      Number of y points.
% 
% Output:
%   Right_Indices   [int]      Row indices of the operator at the right
%                              boundary of the domain.

% Input validation
if ~isnumeric(Nx)
    error('get_operator_right_indices:NxNotNumeric', ...
          'Nx must be numeric.')
end
if ~isnumeric(Ny)
    error('get_operator_right_indices:NyNotNumeric', ...
          'Ny must be numeric.')
end
if mod(Nx, 1) ~= 0 % Check wether Nx is an integer
    error('get_operator_right_indices:NxNotInteger', ...
          ['Nx equals ' num2str(Nx) ', but it must be an integer.'])
end
if mod(Ny, 1) ~= 0 % Check wether Ny is an integer
    error('get_operator_right_indices:NyNotInteger', ...
          ['Ny equals ' num2str(Ny) ', but it must be an integer.'])
end

% Calculate the right indices of an operator
Top_Right_Index    = get_operator_top_right_index   (Nx, Ny);
Bottom_Right_Index = get_operator_bottom_right_index(Nx, Ny);

Right_Indices  = Top_Right_Index : Bottom_Right_Index;

end