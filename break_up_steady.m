function [Population, Rels_steady,rel_dur,ind_new_bracket]=break_up_steady(Population,pop_index,Rels_steady,pairs_break,rel_dur) 
% this function breaks up a steady pair that is listed in Rels_steady table
% under the number_pair counter, and updates Rels_steady and Population and
% rel_dur
    Rels_breakup=Rels_steady.Data(pairs_break,:);
    Rels_steady.Data(pairs_break,:)=[];
    
    % Concatenate IDs from column 1 and 2 of Rels_breakup
    ids = unique(reshape(Rels_breakup(:,1:2)', [], 1));
    
    % Map IDs to indices
    [~, ind_new_bracket] = ismember(ids, Population.Data(pop_index.id,:));
    
    % Remove invalid indices (where ismember could not find a match)
    ind_new_bracket(ind_new_bracket == 0) = [];
    
    % Decrease nsteady for each individual
    Population.Data(pop_index.nsteady, ind_new_bracket) = Population.Data(pop_index.nsteady, ind_new_bracket) - 1;
    
    % Record the duration of each partnership
    rel_dur.Data = [rel_dur.Data; Rels_breakup(:,4) - Rels_breakup(:,3)];
    
    % Remove relationships from the records
    Population.Data(pop_index.sp1_id, ind_new_bracket) = -1;
    Population.Data(pop_index.sp1_stdate, ind_new_bracket) = -40000;
    
end