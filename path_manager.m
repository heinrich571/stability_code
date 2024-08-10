function path_manager(command, additional_paths)

if nargin < 2, additional_paths = {}; end
path_list = {'./baseflow_definition'
             './grid_generation'
             './solver'
             './BatchX'
             './utilities'
             './plots'};

for i = 1:length(path_list)
    switch command
        case {'a','-a','-add','add','Add','ADD','ad','AD','Ad'}
            addpath(genpath(path_list{i}))
        case {'r','-r','-rm','remove','rm','RM','Remove','REMOVE'}
            rmpath(genpath(path_list{i}))
        otherwise
            error('Unsupported command')
    end
end

if ~isempty(additional_paths)
    for i = 1:length(additional_paths)
        switch command
            case {'a','-a','-add','add','Add','ADD','ad','AD','Ad'}
                addpath(genpath(additional_paths{i}))
            case {'r','-r','-rm','remove','rm','RM','Remove','REMOVE'}
                rmpath(genpath(additional_paths{i}))
            otherwise
                error('Unsupported command')
        end
    end
end

end
