function Top_Left_Index = get_equation_top_left_index(Equation_Name, Nx, Ny)

% Description:
%   This function returns the index of the specidied equation at the
%   top-left boundary of the domain in the eigenvalue problem matrices.
% 
% Input:
%   Equation_Name   [char]      Name of the equation.
%   Nx              [int]       Number of x points.
%   Ny              [int]       Number of y points.
% 
% Output:
%   Top_Left_Index  [scalar]    Row indices of specified equation at the
%                               top boundary of the domain in the
%                               eigenvalue problem matrices.

% Input validation
if ~isnumeric(Nx)
    error('get_equation_top_left_index:NxNotNumeric', ...
          'Nx must be numeric.')
end
if ~isnumeric(Ny)
    error('get_equation_top_left_index:NyNotNumeric', ...
          'Ny must be numeric.')
end
if mod(Nx, 1) ~= 0 % Check wether Nx is an integer
    error('get_equation_top_left_index:NxNotInteger', ...
          ['Nx equals ' num2str(Nx) ', but it must be an integer.'])
end
if mod(Ny, 1) ~= 0 % Check wether Ny is an integer
    error('get_equation_top_left_index:NyNotInteger', ...
          ['Ny equals ' num2str(Ny) ', but it must be an integer.'])
end

% Calculate the top left index of the specified equaion
Bottom_Left_Index = get_equation_bottom_left_index(Equation_Name, Nx, Ny);

Top_Left_Index = Bottom_Left_Index - (Ny - 1);

end
