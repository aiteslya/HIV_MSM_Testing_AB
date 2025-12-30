function [all_partnerships, s_fl] = addPartnership(all_partnerships, i, j, start_date, end_date)
    

    % Check if partnership already exists
    if partnershipExists(all_partnerships, i, j)
        s_fl = 0;
        return; % Exit the function if partnership already exists
    end

    % If the partnership does not exist, then proceed to add it
    s_fl = 1;
    if ~isKey(all_partnerships, i)
        all_partnerships(i) = [];
    end
    if ~isKey(all_partnerships, j)
        all_partnerships(j) = [];
    end
    
    partnership_i = createPartnership(j, start_date, end_date);
    partnership_j = createPartnership(i, start_date, end_date);

    current_i_partnerships = all_partnerships(i);
    all_partnerships(i) = [current_i_partnerships, partnership_i];

    current_j_partnerships = all_partnerships(j);
    all_partnerships(j) = [current_j_partnerships, partnership_j];
end
