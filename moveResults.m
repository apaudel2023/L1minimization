function destination = moveResults(folderName)
    fprintf('Moving results to a separate directory ...\n')
    % Get the current directory
    currentDir = pwd;

    % Get the parent directory of the current directory
    [parentDir, ~, ~] = fileparts(currentDir);

    % Create the destination path using the provided folder name
    destination = fullfile(parentDir, folderName);

    % Check if the destination directory exists, if not, create it
    if ~exist(destination, 'dir')
        mkdir(destination);
    end

    % Define file extensions to be copied
    fileExtensions = {'*.mat', '*.fig', '*.png', '*.jpeg', '*.jpg'};

    files = [];

    % Get a list of files for each extension and concatenate them
    for i = 1:length(fileExtensions)
        files = [files; dir(fullfile(currentDir, fileExtensions{i}))];
    end

    % Copy each file to the destination directory
    for i = 1:length(files)
        srcFile = fullfile(files(i).folder, files(i).name);
        destFile = fullfile(destination, files(i).name);
        movefile(srcFile, destFile);
    end

    fprintf('... Done !!!\n')
end
