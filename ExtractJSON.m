% ExtractJSON extracts BrainSense Streaming data from a JSON file
% Author: Albert Guillemette, 02.06.2025

% Select a JSON file to extract
[filename, data_pathname] = uigetfile('*.json', 'Select a .json file');

if isequal(filename, 0)
    error('No file selected. Extraction aborted.');
end

% Load JSON data
fullpath = fullfile(data_pathname, filename);
data = jsondecode(fileread(fullpath));

% Save the output to a .mat file
[~, name, ~] = fileparts(filename);  % Get the file name without extension
save(fullfile(data_pathname, [name '_output.mat']), 'data');