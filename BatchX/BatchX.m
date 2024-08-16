function BatchX(Folders, Filenames, Output_Folder, N_Workers)

% Input checks
if length(Folders) ~= length(Filenames)
    if length(Folders) == 1 % Allow for compact input (in case all the files are situated in the same folder)
        Folders = repmat(Folders, [length(Filenames) 1]);
    else
        error('Number of paths is not equal to the number of file names')
    end
end

% Generate paths to files
Paths = cell([length(Filenames) 1]);
for i = 1:length(Folders)
    Paths{i} = [Folders{i} Filenames{i}];
end

% Create an output folder (in case it does not exist already)
if ~isfolder(Output_Folder)
    mkdir(Output_Folder)
end

% % Create a "parpool"
% parpool(N_Workers)

% Batch execution of runs
for i = 1:length(Paths)
    % Load case
    load(Paths{i}, 'Case_Name', 'Problem')

    dispstatus(['CASE: ' Case_Name], 0)

    % Perform 1st run for initialization
    [Solution, Report] = BiGlobalTemporalSolver(Problem(1));

    Solution = repmat(Solution, [length(Problem) 1]);
    Report   = repmat(Report,   [length(Problem) 1]);

    % Parallel execution
    parfor j = 2:length(Problem)
        [Solution(j), Report(j)] = BiGlobalTemporalSolver(Problem(j));
    end

    % Save results
    output_folder = [Output_Folder '/' Case_Name];
    if ~isfolder(output_folder)
        mkdir(output_folder)
    end
    output_path = [output_folder '/' Case_Name '.mat'];
    save(output_path, 'Problem', 'Solution', 'Report', '-v7.3')
    
    dispstatus(['CASE: ' Case_Name], 1)
    dispstatus()
    
    % Clear memory
    clearvars Problem Solution Report
end

dispstatus('Batch execution finished!', -2)

end
