function [all_partnerships, rel_dur] = handleDeath(all_partnerships, deceased_id, tcounter)
% this function handles death of individuals in the context of cleaning up
% partnerships including the diseased person
    rel_dur = []; % Initialize the durations array

    % Check if deceased person had partnerships
    if isKey(all_partnerships, deceased_id)
        deceased_partnerships = all_partnerships(deceased_id); % Extract partnerships of the deceased person
        
        % Remove the deceased person from each of their partner's lists and calculate durations
        for p = deceased_partnerships
            all_partnerships(p.id) = removeDeceasedPartnership(all_partnerships(p.id), deceased_id);
            
            % Calculate the duration for the deceased person's partnership with this partner
            partnership_duration = min(tcounter - p.start_date, p.end_date - p.start_date);
            rel_dur(end+1) = partnership_duration;

            % Check and remove empty partnerships of the partner
            if isempty(all_partnerships(p.id))
                remove(all_partnerships, p.id);
            end
        end
        
        % Finally, remove the deceased person's entry from all_partnerships
        remove(all_partnerships, deceased_id);
    end
    
    % Helper function to remove deceased from a partnership list
    function updated_list = removeDeceasedPartnership(partnership_list, deceased_id)
        updated_list = partnership_list([partnership_list.id] ~= deceased_id);
    end
end

