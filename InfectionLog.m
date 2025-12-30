classdef InfectionLog < LargeData
    properties
        EventCounter % Keeps track of the number of logged events
    end

    methods
        function obj = InfectionLog(max_events)
            obj.Data = nan(max_events, 17); % Preallocate numeric matrix
            obj.EventCounter = 0; % Initialize event counter
        end

        function logEvent(obj, infector_id, infectee_id, ...
                          infector_index, infectee_index, current_time, ...
                          infector_HIV_status, infectee_HIV_status, ...
                          infector_age, infectee_age, ...
                          infector_propensity, infectee_propensity, ...
                          infector_infection_age, ...
                          infector_casual_partners, infectee_casual_partners, partnershipType, infector_testing_rate, infectee_testing_rate)
            % Increment the event counter before logging
            obj.EventCounter = obj.EventCounter + 1;

            obj.Data(obj.EventCounter, :) = [infector_id, infectee_id, ...
                                          infector_index, infectee_index, ...
                                          current_time, ...
                                          infector_HIV_status, infectee_HIV_status, ...
                                          infector_age, infectee_age, ...
                                          infector_propensity, infectee_propensity, ...
                                          infector_infection_age, ...
                                          infector_casual_partners, infectee_casual_partners, partnershipType, infector_testing_rate, infectee_testing_rate];
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
    %     function obj = loadFromFile(filename)
    %         % Load the object from a .mat file and reconstruct it
    %         if exist(filename, 'file') ~= 2
    %             error('File not found: %s', filename);
    %         end
    % 
    %         tempStruct = load(filename); % Load as struct
    % 
    %         % Ensure required fields exist
    %         if ~isfield(tempStruct, 'Data') || ~isfield(tempStruct, 'EventCounter')
    %             error('Invalid file format: missing required fields.');
    %         end
    % 
    %         % Create a new instance and restore properties
    %         obj = InfectionLog(size(tempStruct.Data, 1));
    %         obj.Data = tempStruct.Data;
    %         obj.EventCounter = tempStruct.EventCounter;
    %     end
    % end
    
        function obj = loadFromFile(filename)
            if exist(filename, 'file') ~= 2
                error('File not found: %s', filename);
            end
        
            tempStruct = load(filename);
        
            % If it was saved as an object
            if isfield(tempStruct, 'infectionLog') && isa(tempStruct.infectionLog, 'InfectionLog')
                obj = tempStruct.infectionLog;
                return;
            end
        
            % If it was saved as a struct (legacy or alternate style)
            if isfield(tempStruct, 'Data') && isfield(tempStruct, 'EventCounter')
                obj = InfectionLog(size(tempStruct.Data, 1));
                obj.Data = tempStruct.Data;
                obj.EventCounter = tempStruct.EventCounter;
            else
                error('Invalid file format: missing required fields.');
            end
        end
    end

end
