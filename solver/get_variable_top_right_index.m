function Top_Right_Index = get_variable_top_right_index(Variable_Name, Nx, Ny)

% Description:
%   This function returns the index of the specified variable at the
%   top-right corner of the domain in the eigenvalue problem matrices.
% 
% Input:
%   Variable_Name       [char]      Name of the variable.
%   Nx                  [int]       Number of x points.
%   Ny                  [int]       Number of y points.
% 
% Output:
%   Top_Right_Index     [int]       Column index of specified variable at
%                                   the top-right corner of the domain in 
%                                   the eigenvalue problem matrices.

% Input validation
if ~isnumeric(Nx)
    error('get_variable_top_right_index:NxNotNumeric', ...
          'Nx must be numeric.')
end
if ~isnumeric(Ny)
    error('get_variable_top_right_index:NyNotNumeric', ...
          'Ny must be numeric.')
end
if mod(Nx, 1) ~= 0 % Check wether Nx is an integer
    error('get_variable_top_right_index:NxNotInteger', ...
          ['Nx equals ' num2str(Nx) ', but it must be an integer.'])
end
if mod(Ny, 1) ~= 0 % Check wether Ny is an integer
    error('get_variable_top_right_index:NyNotInteger', ...
          ['Ny equals ' num2str(Ny) ', but it must be an integer.'])
end

% Calculate the top-right index of the specified variable
Variable_Number = get_variable_block_number(Variable_Name);

Top_Right_Index = 1 + (Nx * Ny) * (Variable_Number - 1);

end