classdef TestingLog < LargeData
    properties
        EventCounter % Keeps track of the number of logged events
    end

    methods
        function obj = TestingLog(max_events)
            obj.Data = nan(max_events, 8); % Preallocate numeric matrix
            obj.EventCounter = 0; % Initialize event counter
        end

        function logEvent(obj, id, index, current_time, ...
                          HIV_status, age, propensity, age_inf, ...
                          casual_partners)
            % Increment the event counter before logging
            obj.EventCounter = obj.EventCounter + 1;

            obj.Data(obj.EventCounter, :) = [id, index, ...
                                          current_time, HIV_status, ...
                                          age, propensity, age_inf, ...
                                          casual_partners];
        end
    
        function trimData(obj)
            % Remove all rows where the first column (InfectorID) is still NaN
            obj.Data(any(isnan(obj.Data), 2), :) = [];
        end

        function lastNEvents = getLastNEvents(obj, N)
            % Retrieve the last N infection events efficiently
            if obj.EventCounter == 0
                warning('No events logged yet.');
                lastNEvents = [];
                return;
            end

            % Ensure N does not exceed the number of logged events
            N = min(N, obj.EventCounter);

            % Extract the last N rows
            lastNEvents = obj.Data(obj.EventCounter - N + 1 : obj.EventCounter, :);
        end

        function saveToFile(obj, filename)
            % Save the object to a .mat file
            tempStruct.Data = obj.Data;
            tempStruct.EventCounter = obj.EventCounter;
            save(filename, '-struct', 'tempStruct'); % Save as struct
        end

    end

    methods (Static)
        function obj = loadFromFile(filename)
            % Load the object from a .mat file and reconstruct it
            if exist(filename, 'file') ~= 2
                error('File not found: %s', filename);
            end
            
            tempStruct = load(filename); % Load as struct
            
            % Ensure required fields exist
            if ~isfield(tempStruct, 'Data') || ~isfield(tempStruct, 'EventCounter')
                error('Invalid file format: missing required fields.');
            end
            
            % Create a new instance and restore properties
            obj = TestingLog(size(tempStruct.Data, 1));
            obj.Data = tempStruct.Data;
            obj.EventCounter = tempStruct.EventCounter;
        end
    end
end
