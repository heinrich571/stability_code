function Top_Right_Index = get_equation_top_right_index(Equation_Name, Nx, Ny)

% Description
%   This function returns the index of the specified equation at the
%   top-right corner of the domain in the eigenvalue problem matrices.
% 
% Input:
%   Equation_Name       [char]      Name of the equation.
%   Nx                  [int]       Number of x points.
%   Ny                  [int]       Number of y points.
% 
% Output:
%   Top_Right_Index     [scalar]    Row index of specified equation at the
%                                   top-right boundary of the domain in the
%                                   eigenvalue problem matrices.

% Input validation
if ~isnumeric(Nx)
    error('get_equation_top_right_index:NxNotNumeric', ...
          'Nx must be numeric.')
end
if ~isnumeric(Ny)
    error('get_equation_top_right_index:NyNotNumeric', ...
          'Ny must be numeric.')
end
if mod(Nx, 1) ~= 0 % Check wether Nx is an integer
    error('get_equation_top_right_index:NxNotInteger', ...
          ['Nx equals ' num2str(Nx) ', but it must be an integer.'])
end
if mod(Ny, 1) ~= 0 % Check wether Ny is an integer
    error('get_equation_top_right_index:NyNotInteger', ...
          ['Ny equals ' num2str(Ny) ', but it must be an integer.'])
end

% Calculate the top right index of the specified equaion
Equation_Block_Number = get_equation_block_number(Equation_Name);

Top_Right_Index = 1 + (Nx * Ny) * (Equation_Block_Number - 1);

end
