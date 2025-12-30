function exists = partnershipExistsPast(all_partnerships, i, j,start_date)

   if ~isKey(all_partnerships, i)
        exists = false;
        return;
    end

    % Loop through the array of structures to check for id and start_date
    partnership_array = all_partnerships(i);
    
    for idx = 1:length(partnership_array)
        if partnership_array(idx).id == j && partnership_array(idx).start_date == start_date
            exists = true;
            return;
        end
    end
    
    exists = false; % If loop completes without finding a match
end