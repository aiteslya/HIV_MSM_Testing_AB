function partners = retrievePartners(all_partnerships, i)
    if isKey(all_partnerships, i)
        partners = all_partnerships(i);
    else
        partners = [];
    end
end