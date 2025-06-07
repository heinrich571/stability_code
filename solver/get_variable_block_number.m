function Variable_Block_Number = get_variable_block_number(Variable_Name)

% Description:
%   This function returns the variable number according to their order in
%   the eigenvalue matrix formulation.
% 
% Input:
%   Variable_Name           [char]      Name of the variable.
% 
% Output:
%   Variable_Block_Number   [scalar]    Number of the variable block
%                                       according to its position in the
%                                       eigenvalue problem formulation.

VARIABLE_NAMES  = {'u'              % --> 1
                   'v'              % --> 2
                   'w'              % --> 3
                   'p'};            % --> 4

% Input validation
if ~ischar(Variable_Name) && ~isstring(Variable_Name)
    error('get_variable_block_number:VariableNameInvalidType', ...
          '"Variable_Name" input is of invalid type.')
end

% Convert variable name to a number
Variable_Block_Number = find(strcmp(VARIABLE_NAMES, Variable_Name));

% Make sure a match has been foun. If not, then throw an error.
if isempty(Variable_Block_Number)
    error('get_variable_block_number:VariableNameInvalid', ...
          'Variable name is invalid')
end

end
