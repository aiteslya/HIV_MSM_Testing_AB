function Num_parts_distr=prepareCasualPropensity(dataset_name, source_folder, file_name)
% this function asks user to provide source of the data (ACS, EMIS, 
% The Network Study),source folder for data file, name of the file
% If no input is provided default data set for EMIS-2017 will be used 
% Number of partners becomes non-dimensionalized to be converted to the
% relative likelihood to acquire partners
% the output has granularity of 52 different entries and its sum is
% normalized to 1

% Set default values if input arguments are not provided
if nargin < 1 || isempty(dataset_name)
    dataset_name = 'EMIS2017';
end
if nargin < 2 || isempty(source_folder)
    source_folder = 'Data';
end
if nargin < 3 || isempty(file_name)
    file_name = 'EMIS_2017.csv';
end

% Create the full path of the file
file_path = fullfile(source_folder, file_name);

% Read the CSV file into a table
try
    data = readtable(file_path);
    fprintf('Dataset "%s" loaded successfully from the file "%s" in the folder "%s".\n', dataset_name, file_name, source_folder);
catch ME
    warning('Unable to read the file "%s" in the folder "%s".\nError: %s', file_name, source_folder, ME.message);
    data = [];
end

myname=mfilename;

switch dataset_name
    case 'EMIS2017'
        Num_parts_distr=EMIS_2017_non_steady_acq_Rate(data);
    case 'ACS'
        warning([myname,': ACS data set is not coded yet.']);
        Num_parts_distr=[];
    case 'NetworkStudy'
        warning([myname,': the Network Study data set is not coded yet.']);
        Num_parts_distr=[];
    otherwise
        warning([myname,': the name of the schema of the data set provided is invalid.']);
        Num_parts_distr=[];
end

end