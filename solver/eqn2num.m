function equation_number = eqn2num(equation_name)

% DESCRIPTION
%   This function returns the equation number according to their order in the eigenvalue matrix formulation
% INPUT
%   equation_name       name of the equation                                                                    [char]
% OUTPUT
%   equation_number     number of the equation according to its position in the eigenvalue problem formulation  [scalar]

equation_names  = {'x momentum' , 'y momentum' , 'z momentum' , 'continuity'};
equation_number = find(strcmp(equation_names, equation_name));

if isempty(equation_number)
    error('Equation name is invalid')
end

end
