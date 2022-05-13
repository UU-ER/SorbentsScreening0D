function [T_list_new] = delete_temperatures_from_list (T_list_old)
    %========================================================
    % Takes a temperature list. If the list is too long, it deletes some
    % entries to make it a maximum length of 6.
    %-------------------------------------------------------
    %Input:     - T_list_old The old temperature list.
    %Output:    - T_list_new The new temperature list: no longer than 6
    %               temperatures.
    %========================================================
    if length(T_list_old) > 6 %if it is still larger than 6
        difference = length(T_list_old) - 6;
        T_list_new = sort(T_list_old);
        index_to_delete = 2;
        while difference > 0
            T_list_new (index_to_delete) =[];
            index_to_delete = index_to_delete + 2;
            if index_to_delete > 6
                index_to_delete = 3;
            end
            difference = difference - 1;
        end
    else
        T_list_new = T_list_old;
    end
end