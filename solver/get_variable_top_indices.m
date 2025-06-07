function Top_Indices = get_variable_top_indices(Variable_Name, Nx, Ny)

% Description:
%   This function returns the indices of the specified variable at the
%   top boundary of the domain in the eigenvalue problem matrices.
% 
% Input:
%   Variable_Name   [char]      Name of the variable.
%   Nx              [int]       Number of x points.
%   Ny              [int]       Number of y points.
% 
% Output:
%   Top_Indices     [scalar]    Column index of specified variable at the
%                               top boundary of the domain in the
%                               eigenvalue problem matrices.

% Input validation
if ~isnumeric(Nx)
    error('get_variable_top_indices:NxNotNumeric', ...
          'Nx must be numeric.')
end
if ~isnumeric(Ny)
    error('get_variable_top_indices:NyNotNumeric', ...
          'Ny must be numeric.')
end
if mod(Nx, 1) ~= 0 % Check wether Nx is an integer
    error('get_variable_top_indices:NxNotInteger', ...
          ['Nx equals ' num2str(Nx) ', but it must be an integer.'])
end
if mod(Ny, 1) ~= 0 % Check wether Ny is an integer
    error('get_variable_top_indices:NyNotInteger', ...
          ['Ny equals ' num2str(Ny) ', but it must be an integer.'])
end

% Calculate the top indices of the specified variable
Top_Right_Index = get_variable_top_right_index(Variable_Name, Nx, Ny);
Top_Left_Index  = get_variable_top_left_index (Variable_Name, Nx, Ny);

Top_Indices  = Top_Right_Index : Ny : Top_Left_Index;

end
