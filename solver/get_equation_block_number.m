function Equation_Block_Number = get_equation_block_number(Equation_Name)

% Description:
%   This function returns the equation number according to their order in
%   the eigenvalue matrix formulation.
% 
% Input:
%   Equation_Name           [char]      Name of the equation.
% 
% Output:
%   Equation_Block_Number   [scalar]    Number of the equation block
%                                       according to its position in the
%                                       eigenvalue problem formulation.

EQUATION_NAMES  = {'x momentum'     % --> 1
                   'y momentum'     % --> 2
                   'z momentum'     % --> 3
                   'continuity'};   % --> 4

% Input validation
if ~ischar(Equation_Name) && ~isstring(Equation_Name)
    error('get_equation_block_number:EquationNameInvalidType', ...
          '"Equation_Name" input is of invalid type.')
end

% Convert equation name to a number
Equation_Block_Number = find(strcmp(EQUATION_NAMES, Equation_Name));

% Make sure a match has been found. If not, then throw an error.
if isempty(Equation_Block_Number)
    error('get_equation_block_number:EquationNameInvalid', ...
          'Equation name is invalid.')
end

end
