function variable_number = var2num(variable_name)

% DESCRIPTION
%   This function returns the variable number according to their order in the eigenvalue vector formulation
% INPUT
%   variable_name       name of the variable                                                                    [char]
% OUTPUT
%   variable_number     number of the variable according to its position in the eigenvalue problem formulation  [scalar]

variable_names  = {'u' , 'v' , 'w' , 'p'};
variable_number = find(strcmp(variable_names, variable_name));

if isempty(variable_number)
    error('Variable name is invalid')
end

end
