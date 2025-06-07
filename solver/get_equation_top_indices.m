function Top_Indices = get_equation_top_indices(Equation_Name, Nx, Ny)

% Description:
%   This function returns the indices of the specidied equation at the
%   top boundary of the domain in the eigenvalue problem matrices.
% 
% Input:
%   Equation_Name   [char]      Name of the equation.
%   Nx              [int]       Number of x points.
%   Ny              [int]       Number of y points.
% 
% Output:
%   Top_Indices     [scalar]    Row indices of specified equation at the
%                               top boundary of the domain in the
%                               eigenvalue problem matrices.

% Input validation
if ~isnumeric(Nx)
    error('get_equation_top_indices:NxNotNumeric', ...
          'Nx must be numeric.')
end
if ~isnumeric(Ny)
    error('get_equation_top_indices:NyNotNumeric', ...
          'Ny must be numeric.')
end
if mod(Nx, 1) ~= 0 % Check wether Nx is an integer
    error('get_equation_top_indices:NxNotInteger', ...
          ['Nx equals ' num2str(Nx) ', but it must be an integer.'])
end
if mod(Ny, 1) ~= 0 % Check wether Ny is an integer
    error('get_equation_top_indices:NyNotInteger', ...
          ['Ny equals ' num2str(Ny) ', but it must be an integer.'])
end

% Calculate the top indices of the specified equaion
Top_Right_Index = get_equation_top_right_index(Equation_Name, Nx, Ny);
Top_Left_Index  = get_equation_top_left_index (Equation_Name, Nx, Ny);

Top_Indices = Top_Right_Index : Ny : Top_Left_Index;

end
