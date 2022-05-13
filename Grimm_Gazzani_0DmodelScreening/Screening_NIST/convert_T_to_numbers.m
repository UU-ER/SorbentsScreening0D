function [T_list_new,T_list_num_new] = convert_T_to_numbers_2(T_list_old)
    %========================================================
    % Takes a temperature list or array and converts this to a list of
    % numbers .
    % (E.g. {TaAa273num356} --> T_list_new = [273], T_list_num_new = 356)
    %-------------------------------------------------------
    %Input:     - T_list_old The old temperature list or cell array, 
    %             consisting of strings.
    %           - list_or_array String: if "list", T_list_old is a list.
    %Output:    - T_list_new The new temperature list, consisting of
    %             integers.
    %           - T_list_num_new List of the numbers of the temperatures
    %========================================================
    
    list_or_array = 'list'; %% I included this file, since it was not found (18.08.21)

    clearvars T_list_num_new specific_num
    if list_or_array == "list"
        for every_temp =[1:length(T_list_old)]
                every_temp_name = strip(string(T_list_old(every_temp)), "T");
                every_temp_name = strip(every_temp_name, "a");
                every_temp_name = strip(every_temp_name, "A");
                every_temp_name = strip(every_temp_name, "a");
                index_num = strfind(every_temp_name,"num");
                T_list_num_new(every_temp) = 0;
                if isempty(index_num) == 0
                    if index_num>0
                        every_temp_name = replaceBetween(every_temp_name,index_num, length(char(every_temp_name)),'');
                        if index_num > 1
                            old_name_T = T_list_old(every_temp);
                            specific_num =  strsplit(old_name_T,'num');
                            T_list_num_new(every_temp) = specific_num (2);
                        end
                    end
                end
                T_list_new (every_temp) = str2num(every_temp_name);
        end
    else
       T_list_old = cellfun(@(x) strip(string(x), "T") ,T_list_old,'UniformOutput',false);
       T_list_old = cellfun(@(x) strip(string(x), "a") ,T_list_old,'UniformOutput',false);
       T_list_old = cellfun(@(x) strip(string(x), "A") ,T_list_old,'UniformOutput',false);
       T_list_old = cellfun(@(x) strip(string(x), "a") ,T_list_old,'UniformOutput',false);
       T_list_old = vertcat(T_list_old{:});
       for every_temp =[1:length(T_list_old)]
                every_temp_name = string(T_list_old(every_temp));
                index_num = strfind(every_temp_name,"num");
                T_list_num_new(every_temp) = 0;
                if isempty(index_num) == 0
                    if index_num>0
                        every_temp_name = replaceBetween(every_temp_name,index_num, length(char(every_temp_name)),'');
                     if index_num > 1
                            old_name_T = T_list_old(every_temp);
                            specific_num =  strsplit(old_name_T,'num');
                            T_list_num_new(every_temp) = specific_num (2);
                        end
                    end
                end
                T_list_new (every_temp) = str2num(every_temp_name);
       end
    end
end