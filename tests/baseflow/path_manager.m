function path_manager(command)

folder_list = {'../../baseflow_definition/'
               '../../grid_generation/'};                

switch command
    case {'add' 'a' 'Add' 'ADD'}
        for i = 1 : length(folder_list)
            addpath(folder_list{i})
        end
    case {'rm' 'r' 'remove' 'Remove' 'REMOVE'}
        for i = 1 : length(folder_list)
            rmpath(folder_list{i})
        end
    otherwise
        error(['path_manager - command' command ' is not supported'])
end

end
