function [casual_dur, all_casual_parts]=break_up_casual(casual_dur, all_casual_parts, tcounter) 
% this function breaks up a casual pair that is listed in Rels_casual table
% appearing in the list pairs_break, and updates Rels_casual and Pop_casual and
% casual_dur

    Rels_breakup_temp = [];

    % remove partnerships from the map
    keysToRemove = [];
    
    % Extract keys into a variable for clarity
    keysList = keys(all_casual_parts);
    
    for k = 1:length(keysList)
        currentKey = keysList{k};
        updated_partnerships = [];
        
        for p = all_casual_parts(currentKey)
            if p.end_date > tcounter
                updated_partnerships = [updated_partnerships, p];
            else
                % If partnership is removed, calculate its duration and add to the durations array
                partnership_duration = p.end_date - p.start_date;
                Rels_breakup_temp = [Rels_breakup_temp; sort([currentKey p.id]) p.start_date p.end_date];
            end
        end
        
        if isempty(updated_partnerships) % no partnerships left, remove the entry for this person from the registry
            keysToRemove(end+1) = currentKey;
        else
            all_casual_parts(currentKey) = updated_partnerships;
        end
    end
    
   % Remove keys with no partnerships
    for i = 1:length(keysToRemove)
        if isKey(all_casual_parts, keysToRemove(i))
            remove(all_casual_parts, keysToRemove(i));
        end
    end
    
    % break up the pairs
    if ~isempty(Rels_breakup_temp)
        % Find unique partnerships that broke up
        [uniqueRels, ~] = unique(Rels_breakup_temp, 'rows', 'stable');
       
        durations = uniqueRels(:, 4) - uniqueRels(:, 3); % end_date - start_date
        
        % Store the durations in casual_dur2.Data
        casual_dur.Data = [casual_dur.Data ; durations];
    end
end