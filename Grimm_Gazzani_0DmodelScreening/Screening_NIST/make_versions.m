%=====================================
% If flag_with_versions == 1, makes multiple versions:
% E.g. T273, T298, T298num1 --> v1 T273 T298 and v2 T273 T298.
% 
% NOTE: Advised to only run when looking at one specific material;
%   because many versions can be made and running can take very long.
%--------------------------------------


% --- make versions ---
 clearvars struct_T_version_specific  T_list_counts  T_list_num
 T_for_all_versions = [];   T_version_specific = [];
 [T_list_num, T_num]= convert_T_to_numbers_2(T_list,"list");
 T_list_counts = histc(T_list_num(:), unique(T_list_num));
 T_unique = unique(T_list_num);

 if sum(T_list_counts) ~= length(T_list_counts) % if one of the counts is >1 
     count_iterations = 0;
     for every_count_in_T_list = [1:length(T_list_counts)]
            % don't know what is meant here -> I changed it
%             temp = edges(every_count_in_T_list);
         temp = T_unique(every_count_in_T_list);
         index_T = find(T_list_num ==temp);

         % get all temperatures that should be in all versions
         if T_list_counts(every_count_in_T_list) == 1 % only once in list
             if T_num(index_T) ~= 0 % there was a num originally
                add_temp = strcat("TaAa", string(temp), "num", string(T_num(index_T)));
             else
                add_temp = strcat("TaAa",string(temp));
             end
             T_for_all_versions =[T_for_all_versions add_temp] 
             for every_index_T = [1:length(index_T)]
                struct_T_version_specific.(strcat("TaAa", string(temp))) = add_temp
             end
         else % if a T is more than once in the list, versions need to be made
             if length (index_T) >1
                 for every_index_T = [1:length(index_T)]
                     count_iterations = count_iterations +1;
                     if T_num(index_T(every_index_T)) ~= 0 % there was a num originally
                        add_temp = strcat("TaAa", string(temp), "num", string(T_num(index_T(every_index_T))));
                     else
                        add_temp = strcat("TaAa", string(temp));
                     end
                     if count_iterations > 1 
                         if temp ~= old_temp
                            T_version_specific = [];
                         end
                     end
                     T_version_specific = [T_version_specific add_temp] ;
                     struct_T_version_specific.(strcat("TaAa", string(temp))) = T_version_specific;
                     old_temp = temp;
                 end
             end
         end
     end
     
     %now make the versions:
     list_fieldnames = fieldnames(struct_T_version_specific);
     final_variable_count = 1;
     field_count = 0;
     
     for i = [1:length(T_list_counts)]
        final_variable_count = final_variable_count*T_list_counts(i);
     end
         count_old = 1;
         for fieldname_i = [1:length(list_fieldnames)]
             field_count = field_count + 1;
             count_new = T_list_counts(fieldname_i) + count_old-1;
             T_fieldname_in_struct_version = string(list_fieldnames(fieldname_i));
             list_temperatures = struct_T_version_specific.(T_fieldname_in_struct_version) 
             eval(['list_temperatures_' num2str(fieldname_i) '= list_temperatures']);
             if final_variable_count <= final_variable_count %by changing final_variable_count in e.g. 3, the maximum amount of combinations per DOI can be determined
                list_temperatures_for_comb = [count_old : count_new]
                eval(['list_temperatures_for_comb_' num2str(fieldname_i) '= list_temperatures_for_comb']);
             else %dpeending on value final_variable_count. If chosen to change to 3, if more than 3 combinations, take 1 combination only
                list_temperatures_for_comb = [count_old]
                eval(['list_temperatures_for_comb_' num2str(fieldname_i) '= list_temperatures_for_comb']);
             end
             count_old = count_new+1;
         end

         if field_count == 1
            combos = combvec(list_temperatures_for_comb_1);
         elseif field_count == 2
            combos = combvec(list_temperatures_for_comb_1,list_temperatures_for_comb_2);
         elseif field_count == 3
             combos = combvec(list_temperatures_for_comb_1,list_temperatures_for_comb_2,list_temperatures_for_comb_3);
         elseif field_count == 4
             combos = combvec(list_temperatures_for_comb_1,list_temperatures_for_comb_2,list_temperatures_for_comb_3,list_temperatures_for_comb_4) ;        
         elseif field_count == 5
             combos = combvec(list_temperatures_for_comb_1,list_temperatures_for_comb_2,list_temperatures_for_comb_3,list_temperatures_for_comb_4,list_temperatures_for_comb_5);
         elseif field_count == 6
             combos = combvec(list_temperatures_for_comb_1,list_temperatures_for_comb_2,list_temperatures_for_comb_3,list_temperatures_for_comb_4,list_temperatures_for_comb_5,list_temperatures_for_comb_6);
         elseif field_count == 7
             combos = combvec(list_temperatures_for_comb_1,list_temperatures_for_comb_2,list_temperatures_for_comb_3,list_temperatures_for_comb_4,list_temperatures_for_comb_5,list_temperatures_for_comb_6,list_temperatures_for_comb_7);
         elseif field_count == 8
             combos = combvec(list_temperatures_for_comb_1,list_temperatures_for_comb_2,list_temperatures_for_comb_3,list_temperatures_for_comb_4,list_temperatures_for_comb_5,list_temperatures_for_comb_6,list_temperatures_for_comb_7,list_temperatures_for_comb_8);
         elseif field_count == 9
             combos = combvec(list_temperatures_for_comb_1,list_temperatures_for_comb_2,list_temperatures_for_comb_3,list_temperatures_for_comb_4,list_temperatures_for_comb_5,list_temperatures_for_comb_6,list_temperatures_for_comb_7,list_temperatures_for_comb_8,list_temperatures_for_comb_9);
         elseif field_count == 10
             combos = combvec(list_temperatures_for_comb_1,list_temperatures_for_comb_2,list_temperatures_for_comb_3,list_temperatures_for_comb_4,list_temperatures_for_comb_5,list_temperatures_for_comb_6,list_temperatures_for_comb_7,list_temperatures_for_comb_8,list_temperatures_for_comb_9,list_temperatures_for_comb_10);
         elseif field_count == 11
             combos = combvec(list_temperatures_for_comb_1,list_temperatures_for_comb_2,list_temperatures_for_comb_3,list_temperatures_for_comb_4,list_temperatures_for_comb_5,list_temperatures_for_comb_6,list_temperatures_for_comb_7,list_temperatures_for_comb_8,list_temperatures_for_comb_9,list_temperatures_for_comb_10,list_temperatures_for_comb_11);
         elseif field_count == 12
             combos = combvec(list_temperatures_for_comb_1,list_temperatures_for_comb_2,list_temperatures_for_comb_3,list_temperatures_for_comb_4,list_temperatures_for_comb_5,list_temperatures_for_comb_6,list_temperatures_for_comb_7,list_temperatures_for_comb_8,list_temperatures_for_comb_9,list_temperatures_for_comb_10,list_temperatures_for_comb_11,list_temperatures_for_comb_12);
         elseif field_count == 13
             combos = combvec(list_temperatures_for_comb_1,list_temperatures_for_comb_2,list_temperatures_for_comb_3,list_temperatures_for_comb_4,list_temperatures_for_comb_5,list_temperatures_for_comb_6,list_temperatures_for_comb_7,list_temperatures_for_comb_8,list_temperatures_for_comb_9,list_temperatures_for_comb_10,list_temperatures_for_comb_11,list_temperatures_for_comb_12,list_temperatures_for_comb_13);
         elseif field_count == 14
             combos = combvec(list_temperatures_for_comb_1,list_temperatures_for_comb_2,list_temperatures_for_comb_3,list_temperatures_for_comb_4,list_temperatures_for_comb_5,list_temperatures_for_comb_6,list_temperatures_for_comb_7,list_temperatures_for_comb_8,list_temperatures_for_comb_9,list_temperatures_for_comb_10,list_temperatures_for_comb_11,list_temperatures_for_comb_12,list_temperatures_for_comb_13,list_temperatures_for_comb_14);
         elseif field_count == 15
             combos = combvec(list_temperatures_for_comb_1,list_temperatures_for_comb_2,list_temperatures_for_comb_3,list_temperatures_for_comb_4,list_temperatures_for_comb_5,list_temperatures_for_comb_6,list_temperatures_for_comb_7,list_temperatures_for_comb_8,list_temperatures_for_comb_9,list_temperatures_for_comb_10,list_temperatures_for_comb_11,list_temperatures_for_comb_12,list_temperatures_for_comb_13,list_temperatures_for_comb_14,list_temperatures_for_comb_15);
         elseif field_count == 16
             combos = combvec(list_temperatures_for_comb_1,list_temperatures_for_comb_2,list_temperatures_for_comb_3,list_temperatures_for_comb_4,list_temperatures_for_comb_5,list_temperatures_for_comb_6,list_temperatures_for_comb_7,list_temperatures_for_comb_8,list_temperatures_for_comb_9,list_temperatures_for_comb_10,list_temperatures_for_comb_11,list_temperatures_for_comb_12,list_temperatures_for_comb_13,list_temperatures_for_comb_14,list_temperatures_for_comb_15,list_temperatures_for_comb_16);
         elseif field_count == 17
             combos = combvec(list_temperatures_for_comb_1,list_temperatures_for_comb_2,list_temperatures_for_comb_3,list_temperatures_for_comb_4,list_temperatures_for_comb_5,list_temperatures_for_comb_6,list_temperatures_for_comb_7,list_temperatures_for_comb_8,list_temperatures_for_comb_9,list_temperatures_for_comb_10,list_temperatures_for_comb_11,list_temperatures_for_comb_12,list_temperatures_for_comb_13,list_temperatures_for_comb_14,list_temperatures_for_comb_15,list_temperatures_for_comb_16,list_temperatures_for_comb_17);
          elseif field_count == 18
             combos = combvec(list_temperatures_for_comb_1,list_temperatures_for_comb_2,list_temperatures_for_comb_3,list_temperatures_for_comb_4,list_temperatures_for_comb_5,list_temperatures_for_comb_6,list_temperatures_for_comb_7,list_temperatures_for_comb_8,list_temperatures_for_comb_9,list_temperatures_for_comb_10,list_temperatures_for_comb_11,list_temperatures_for_comb_12,list_temperatures_for_comb_13,list_temperatures_for_comb_14,list_temperatures_for_comb_15,list_temperatures_for_comb_16,list_temperatures_for_comb_17,list_temperatures_for_comb_18);
         elseif field_count == 19
             combos = combvec(list_temperatures_for_comb_1,list_temperatures_for_comb_2,list_temperatures_for_comb_3,list_temperatures_for_comb_4,list_temperatures_for_comb_5,list_temperatures_for_comb_6,list_temperatures_for_comb_7,list_temperatures_for_comb_8,list_temperatures_for_comb_9,list_temperatures_for_comb_10,list_temperatures_for_comb_11,list_temperatures_for_comb_12,list_temperatures_for_comb_13,list_temperatures_for_comb_14,list_temperatures_for_comb_15,list_temperatures_for_comb_16,list_temperatures_for_comb_17,list_temperatures_for_comb_18,list_temperatures_for_comb_19);
         elseif field_count == 20
             combos = combvec(list_temperatures_for_comb_1,list_temperatures_for_comb_2,list_temperatures_for_comb_3,list_temperatures_for_comb_4,list_temperatures_for_comb_5,list_temperatures_for_comb_6,list_temperatures_for_comb_7,list_temperatures_for_comb_8,list_temperatures_for_comb_9,list_temperatures_for_comb_10,list_temperatures_for_comb_11,list_temperatures_for_comb_12,list_temperatures_for_comb_13,list_temperatures_for_comb_14,list_temperatures_for_comb_15,list_temperatures_for_comb_16,list_temperatures_for_comb_17,list_temperatures_for_comb_18,list_temperatures_for_comb_19,list_temperatures_for_comb_20);
         elseif field_count == 21
             combos = combvec(list_temperatures_for_comb_1,list_temperatures_for_comb_2,list_temperatures_for_comb_3,list_temperatures_for_comb_4,list_temperatures_for_comb_5,list_temperatures_for_comb_6,list_temperatures_for_comb_7,list_temperatures_for_comb_8,list_temperatures_for_comb_9,list_temperatures_for_comb_10,list_temperatures_for_comb_11,list_temperatures_for_comb_12,list_temperatures_for_comb_13,list_temperatures_for_comb_14,list_temperatures_for_comb_15,list_temperatures_for_comb_16,list_temperatures_for_comb_17,list_temperatures_for_comb_18,list_temperatures_for_comb_19,list_temperatures_for_comb_20,list_temperatures_for_comb_21);
         elseif field_count == 22
             combos = combvec(list_temperatures_for_comb_1,list_temperatures_for_comb_2,list_temperatures_for_comb_3,list_temperatures_for_comb_4,list_temperatures_for_comb_5,list_temperatures_for_comb_6,list_temperatures_for_comb_7,list_temperatures_for_comb_8,list_temperatures_for_comb_9,list_temperatures_for_comb_10,list_temperatures_for_comb_11,list_temperatures_for_comb_12,list_temperatures_for_comb_13,list_temperatures_for_comb_14,list_temperatures_for_comb_15,list_temperatures_for_comb_16,list_temperatures_for_comb_17,list_temperatures_for_comb_18,list_temperatures_for_comb_19,list_temperatures_for_comb_20,list_temperatures_for_comb_21,list_temperatures_for_comb_22);
         elseif field_count == 23
             combos = combvec(list_temperatures_for_comb_1,list_temperatures_for_comb_2,list_temperatures_for_comb_3,list_temperatures_for_comb_4,list_temperatures_for_comb_5,list_temperatures_for_comb_6,list_temperatures_for_comb_7,list_temperatures_for_comb_8,list_temperatures_for_comb_9,list_temperatures_for_comb_10,list_temperatures_for_comb_11,list_temperatures_for_comb_12,list_temperatures_for_comb_13,list_temperatures_for_comb_14,list_temperatures_for_comb_15,list_temperatures_for_comb_16,list_temperatures_for_comb_17,list_temperatures_for_comb_18,list_temperatures_for_comb_19,list_temperatures_for_comb_20,list_temperatures_for_comb_21,list_temperatures_for_comb_22,list_temperatures_for_comb_23);    
         end
         clearvars combos_strings
         sum_T_list_counts = cumsum( T_list_counts); %Example: Temperature 1 appears 6 times, etc: 6 5 5 5 --> 6 11 17 21
         for j = [1:size(combos,2)] %for every specific combination
             specific_variant = combos(:,j); % e.g. 1 7 12 17
             clearvars specific_variant_strings
             for fieldname_i = [1:length(specific_variant)] %for every temperature set in the combination (e.g. 7)
                 clearvars list_temperatures
                 T_fieldname_in_struct_version = string(list_fieldnames(fieldname_i));
                 list_temperatures = struct_T_version_specific.(T_fieldname_in_struct_version) ;
                 difference = sum_T_list_counts(fieldname_i) - specific_variant(fieldname_i); % 11-7 =4 
                 specific_variant_strings(fieldname_i) = list_temperatures(end-difference); %take the right temperature form the list        
             end
             combos_strings(:,j) = specific_variant_strings;
         end
        clearvars combos
        combos_strings ;
        flag_multiple_versions = 1;
 end

 %--- make a loop over every version ---
 if sum(T_list_counts) ~= length(T_list_counts) && flag_multiple_versions == 1 
     if size(combos_strings,2) >= size(combos_strings,2) 
         end_j = size(combos_strings,2) 
     else
         end_j = size(combos_strings,2);
     end
     for j = [1:end_j]
          T_list = combos_strings(:,j);
          material_version = strcat(material_version_original, ", variant: ", string(j));
          fitting
     end
     flag_fitted_because_of_variants = 1;
 end