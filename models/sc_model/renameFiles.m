% Set the directory path
dir_path = 'noiseModelOutput/';

% Get a list of all files in the directory
file_list = dir(fullfile(dir_path, '*.mat'));

% Loop through each file in the directory
for i = 1:length(file_list)
    % Get the file name
    file_name = file_list(i).name;
    
    % Check if the file name matches the pattern
    if ~isempty(regexp(file_name, '^\w+_\d{2}_\w+', 'once'))
        % Extract the two-digit number from the file name
        old_number = regexp(file_name, '(\d{2})', 'match', 'once');
        
        % Convert the two-digit number to a three-digit number
        new_number = sprintf('%03d', str2double(old_number));
        
        % Replace the two-digit number with the three-digit number in the file name
        new_file_name = regexprep(file_name, '(\d{2})', new_number, 'once');
        
        % Rename the file
        movefile(fullfile(dir_path, file_name), fullfile(dir_path, new_file_name));
    end
end
