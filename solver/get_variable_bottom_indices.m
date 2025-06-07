function Bottom_Indices = get_variable_bottom_indices(Variable_Name, Nx, Ny)

% Description:
%   This function returns the indices of the specified variable at the
%   bottom boundary of the domain in the eigenvalue problem matrices.
% 
% Input:
%   Variable_Name   [char]      Name of the variable.
%   Nx              [int]       Number of x points.
%   Ny              [int]       Number of y points.
% 
% Output:
%   Bottom_Indices  [scalar]    Column index of specified variable at the
%                               bottom boundary of the domain in the
%                               eigenvalue problem matrices.

% Input validation
if ~isnumeric(Nx)
    error('get_variable_bottom_indices:NxNotNumeric', ...
          'Nx must be numeric.')
end
if ~isnumeric(Ny)
    error('get_variable_bottom_indices:NyNotNumeric', ...
          'Ny must be numeric.')
end
if mod(Nx, 1) ~= 0 % Check wether Nx is an integer
    error('get_variable_bottom_indices:NxNotInteger', ...
          ['Nx equals ' num2str(Nx) ', but it must be an integer.'])
end
if mod(Ny, 1) ~= 0 % Check wether Ny is an integer
    error('get_variable_bottom_indices:NyNotInteger', ...
          ['Ny equals ' num2str(Ny) ', but it must be an integer.'])
end

% Calculate the bottom indices of the specified variable
Bottom_Right_Index = get_variable_bottom_right_index(Variable_Name, Nx, Ny);
Bottom_Left_Index  = get_variable_bottom_left_index (Variable_Name, Nx, Ny);

Bottom_Indices = Bottom_Right_Index : Ny : Bottom_Left_Index;

end
