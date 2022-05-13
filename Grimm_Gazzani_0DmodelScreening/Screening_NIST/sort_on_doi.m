%=====================================
% Makes a temperature list T_list with the right temperatures. Calls
% make_versions if flag_with_versions == 1.
%--------------------------------------

%--- initialize flags ---
flag_not_fitted = 0; 
flag_error_fitting_input =0;
flag_doi_unique =0;
flag_multiple_versions =0;
flag_version_has_only_one_T = 0;
flag_doi_not_unique = 0;


%--- Only fit the temperatures belonging to the specific doi ---
for every_temp = [1:length(T_list)]
    temperature_name = T_list(every_temp);
    if isempty(data_NIST.(material_struc_name).(every_gas_name).(temperature_name))
     
    else
    
        if ismember(data_NIST.(material_struc_name).(every_gas_name).(temperature_name).Doi_NIST, material_version) == 1  
            flag_doi_not_unique = 1;                  % doi is not unique --> proceed to fitting
            T_list_2 = [T_list_2 T_list(every_temp)]; % list of Ts with doi not unique
                % check how many datapoints for this specific temperature
                number_list = [number_list length(data_NIST.(material_struc_name).(every_gas_name).(temperature_name).pressure_NIST)];
        end
    end
end

%--- if doi that is not unique new T_list ---
if flag_doi_not_unique == 1
    T_list = T_list_2;
end

%---if all temperatures in T_list are the same, do not fit! ---
T_list_to_check = convert_T_to_numbers(T_list);
unique_list_to_check = unique(T_list_to_check);
if length(unique_list_to_check) == 1
    flag_version_has_only_one_T = 1;
    flag_not_fitted =1;
end

material_version_original = material_version;

%--- check if the temperatures are in Celsius, if so, convert to Kelvin ---
if flag_version_has_only_one_T == 0 % only if more than one T!
    flag_celsius = 0;
    T_list_to_check = convert_T_to_numbers(T_list);
    average_T_list = mean(T_list_to_check);

    % if average of all temperatures is below 70, it is assumed that the 
    % unit is Celsius. The next lines convert them into Kelvin
   if average_T_list < 70 
         flag_celsius = 1;
         T_list_in_Celsius = T_list;
         for every_temp =[1:length(T_list_to_check)]
             old_temp = T_list_to_check(every_temp);
             new_temp = T_list_to_check(every_temp) + 273 ;% convert to Kelvin 
             temp_name = strrep(T_list(every_temp),string(old_temp),string(new_temp));
             T_list_3 = [ T_list_3, temp_name];
         end
         T_list = T_list_3;
   end
    
    % only carry out when screening of databse!!! (takes very long) !!!!!!
        % next time adapt the make_versions file. At the moment there are
        % sometimes more than 50000 possible combinations, because there
        % are many different temperatures -> reduce them in an intelligent
        % way!!!
    flag_with_versions = 0;
   
   if flag_with_versions == 1
        make_versions % only if flag_with_versions == 1 
   end
   
% if there is more data for the same temperature, chose only one
% make index: ones = take temp., zero = don't use
% if length(T_list_to_check) > length(unique_list_to_check)
%     
%     u = unique(T_list_to_check(~isnan(T_list_to_check)));
%     n = histc(T_list_to_check,u);
%     d = u(n > 1);
%     for i = 1:length(d)
%         out(i).T = d(i);
%         out(i).row = find(ismember(T_list_to_check,d(i)));
%         out(i).noEx = number_list(out(i).row);
%         [M,I] = max(number_list(out(i).row));
%         out(i).choose = out(i).row(I);
%         
%         for l=1:length(out(i).row)
%             if l==out(i).choose
%             else
%                 temperature_name = T_list(out(i).row(l)); 
%                 data_NIST.(material_struc_name).(every_gas_name).(temperature_name) = [];
%             end
%         end
%     end
% end
   
end