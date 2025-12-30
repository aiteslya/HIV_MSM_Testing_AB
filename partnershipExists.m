function exists = partnershipExists(all_partnerships, i, j)

   if ~isKey(all_partnerships, i)
        exists = false;
        return;
    end

    % Loop through the array of structures to extract the ids
    partnership_array = all_partnerships(i);
    exists = false; % Default value

    for idx = 1:length(partnership_array)
        if partnership_array(idx).id == j
            exists = true;
            break;
        end
    end

end