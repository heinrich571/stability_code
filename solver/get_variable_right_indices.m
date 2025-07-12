function Right_Indices = get_variable_right_indices(Variable_Name, Nx, Ny)

% Description:
%   This function returns the indices of the specidied variable at the
%   right boundary of the domain in the eigenvalue problem matrices.
% 
% Input:
%   Variable_Name   [char]      Name of the variable.
%   Nx              [int]       Number of x points.
%   Ny              [int]       Number of y points.
% 
% Output:
%   Right_Indices   [scalar]    Row indices of specified variable at the
%                               right boundary of the domain in the
%                               eigenvalue problem matrices.

% Input validation
if ~isnumeric(Nx)
    error('get_variable_right_indices:NxNotNumeric', ...
          'Nx must be numeric.')
end
if ~isnumeric(Ny)
    error('get_variable_right_indices:NyNotNumeric', ...
          'Ny must be numeric.')
end
if mod(Nx, 1) ~= 0 % Check wether Nx is an integer
    error('get_variable_right_indices:NxNotInteger', ...
          ['Nx equals ' num2str(Nx) ', but it must be an integer.'])
end
if mod(Ny, 1) ~= 0 % Check wether Ny is an integer
    error('get_variable_right_indices:NyNotInteger', ...
          ['Ny equals ' num2str(Ny) ', but it must be an integer.'])
end

% Calculate the right indices of the specified variable
Top_Right_Index    = get_variable_top_right_index   (Variable_Name, Nx, Ny);
Bottom_Right_Index = get_variable_bottom_right_index(Variable_Name, Nx, Ny);

Right_Indices = Top_Right_Index : Bottom_Right_Index;

end