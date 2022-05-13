function [gas_display, row_length_table_NIST_1, table_number_of_T, material_name, table_mat_temps] = analyze_data (data_NIST)
    %=====================================
    % Analyzes data structure, gives information on what data is in there.
    %--------------------------------------
    %Input:  - data_NIST A data structure with the structure
    %          materials-gasses - temperatures (in which DOI is also stored)
    %Output: - gas_display (list) list of all gasses in the database
    %        - row_length_table_NIST_1 Total number of different materials
    %        - table_number_of_T Table with number of unique temperatures 
    %           for every gas - material combination 
    %           (e.g. material 1, gas 1, 2 (e.g. T273, T293))
    %        - material_name List of all names of the materials
    %        - table_mat_temps Table with number of non-unique temperatures
    %           (e.g. material 1, gas 1, 3 (e.g. T273, T293, T293#2) )
    %--------------------------------------

fprintf("≡ analyze_data ≡\n")

names_materials = fieldnames(data_NIST);
table_data_NIST_1 = names_materials;

size_table_NIST_1 = size(table_data_NIST_1);
row_length_table_NIST_1 = size_table_NIST_1(1);

count_num_gasses_for_mat_1 = 1;
flag_extra_gas = 0; 
count_material_CO2_combinations = 0;

%--- make a list of the gasses ---
for material_number = [1:row_length_table_NIST_1] %for every material (every row in gas)(26)
    material_name(material_number) = string(table_data_NIST_1(material_number,1));
    place_gas_in_structure = data_NIST.(material_name(material_number));
    names_gasses_1_material  = fieldnames(place_gas_in_structure);
    if material_number == 1
        delete_full_name_entry = names_gasses_1_material(2:end,:);
        gas_display = string(delete_full_name_entry);
    end
    number_of_gasses_material(material_number) = length(names_gasses_1_material);
   
   for gas_number = [2:number_of_gasses_material(material_number)] %for every column in gas (7)
        current_gas = string(names_gasses_1_material(gas_number));

        if material_number == 1
             count_num_gasses_for_mat_1 = count_num_gasses_for_mat_1 + 1 ;
             gas(material_number,gas_number) = current_gas;
        elseif ismember(current_gas,gas)==1; % say if mat 2 has O2 on place 4, but mat 1 has O2 on place 3, synchronize these places
            if flag_extra_gas == 0
                synchronize_j = find(gas(1,:) == current_gas);
                 gas(material_number,synchronize_j) = current_gas;
            else
                for k = [1:material_number]
                    if isempty(find(gas(k,:) == current_gas)) == 0
                       synchronize_j = find(gas(k,:) == current_gas);
                       gas(material_number,synchronize_j) = current_gas;
                    end
                end
            end
        else
            count_num_gasses_for_mat_1 = count_num_gasses_for_mat_1 + 1;
            gas(material_number,count_num_gasses_for_mat_1) = current_gas;
            gas_display(end+1) = [string(current_gas)];
            flag_extra_gas = 1;
        end 
        
   end
    
end

%find which gas/material combinations actually have data
size_gas = size(gas);
number_of_columns = size_gas(2);

doi_list = [];
doi_list_CO2 = [];

number_of_rows = size_gas(1);
gas = gas(:,2:number_of_columns); %omit first column (that specified Full name)
number_of_columns = number_of_columns - 1;
for material_number = [1:number_of_columns] 
    for gas_number = [1:number_of_rows]
        TFmissing = ismissing(gas(gas_number,material_number)); %(1 if missing, 0 if not missing)
        if TFmissing ==1
            TFgas(gas_number,material_number) = 0;
        else
            TFgas(gas_number,material_number) =1;
        end
    end
end
debug_catch_of_try = 0;
%%
for every_gas = [1:length(gas_display)] %%%%%%%%%%%% I changed this since I only wanted the CO2
    for material_number = [1:row_length_table_NIST_1] %for every material (every row in gas)(26)
        start_j = 1;
        try
            place_T_in_structure = data_NIST.(material_name(material_number)).(gas_display(every_gas));
            if gas_display(every_gas) == "Carbon_Dioxide"
                count_material_CO2_combinations = count_material_CO2_combinations+1;
            end

            unique_temperatures_1_gas = fieldnames(place_T_in_structure); % get the temperatures
            
            %make doi list of all possible dois
            doi_list_one_mat = [];
            doi_list_CO2_one_mat = [];
            flag_CO2 = 0;
            for temperature = [1: length(unique_temperatures_1_gas)]
                temperature_name = string(unique_temperatures_1_gas(temperature));
                doi = data_NIST.(material_name(material_number)).(gas_display(every_gas)).(temperature_name).Doi_NIST; % get the doi
                doi_list_one_mat = [doi_list_one_mat doi]; % save all the dois in one vector
                if gas_display(every_gas) == "Carbon_Dioxide"
                    flag_CO2 = 1;
                    doi_list_CO2_one_mat = [doi_list_CO2_one_mat doi]; % save all the dois in one vector special for CO2
                end
            end
            unique_doi_list_one_mat =  unique(doi_list_one_mat); % because doi belonging to several materials (e.g. 3) is counted multiple times (e.g. 3), if all dois are the same, then just one is saved here
            if flag_CO2 ==1
                unique_doi_list_CO2_one_mat = unique(doi_list_CO2_one_mat);
            end
            
            doi_list = [doi_list, unique_doi_list_one_mat];
            if flag_CO2 ==1
                doi_list_CO2 = [doi_list_CO2, unique_doi_list_CO2_one_mat];
            end

            %materials for which there is 1 temperature will be deleted.
            %Also if there are multiple verions of this material: e.g.
            %T220#1 T220#2 T220#3. 
            
            temperatures_1_gas_original = unique_temperatures_1_gas;
            
            unique_temperatures_1_gas = convert_T_to_numbers_2(unique_temperatures_1_gas, 'array');
            unique_temperatures_1_gas = unique(unique_temperatures_1_gas);
            unique_temperatures_1_gas = num2cell(unique_temperatures_1_gas);
            
            table_number_of_T(material_number,every_gas) = length(unique_temperatures_1_gas); %number of unique temperatures: e.g. T220, 225, 330
            number_of_non_unique_isotherms(material_number,every_gas) = length(temperatures_1_gas_original); %e.g. T220#1 T220#2 T220#3 T225, T330#1, T330#2 T400

            for j = [1:length(temperatures_1_gas_original)]
                table_mat_temps (every_gas ,material_number, j) = temperatures_1_gas_original(j); % GAS, MATERIAL, T
            end
        catch
            debug_catch_of_try = debug_catch_of_try + 1;
            number_of_non_unique_isotherms(material_number,every_gas) = 0;
        end

    end

end


%printing analysis: total number of isotherms includes T220#1 T220#2 T220#3
%T225, T330#1, T330#2 T400 hence number_of_non_unique_isotherms
number_of_isotherms_per_gas = sum(number_of_non_unique_isotherms,1); 
number_of_iotherms_per_material = sum(number_of_non_unique_isotherms,2);
total_isotherms_comp_gas = 0;


for every_gas = [1:length(gas_display)] %--> again onlye CO2
    total_isotherms_comp_gas = total_isotherms_comp_gas + number_of_isotherms_per_gas(every_gas);
    if gas_display(every_gas) == "Carbon_Dioxide"
        isotherms_CO2 = number_of_isotherms_per_gas(every_gas);
    end
end

doi_list_unique = unique(doi_list) ; %number of unique dois
doi_list_CO2_unique = unique(doi_list_CO2); % number of unique dois for CO2 

fprintf("Total number of isotherms  %d  \n", total_isotherms_comp_gas)
fprintf("Total number of different materials  %d  \n", row_length_table_NIST_1 )
fprintf("Total number of different DOIs  %d  \n", length(doi_list_unique) )
fprintf("Total number of DOIs (counted multiple times if DOI for multiple materials) %d  \n", length(doi_list) )
fprintf("Total number of isotherms for CO2  %d  \n", isotherms_CO2)
fprintf("Total number of materials for CO2  in new structure %d  \n", count_material_CO2_combinations)
fprintf("Total number of different DOIs for CO2  %d  \n", length(doi_list_CO2_unique) )
fprintf("Total number of DOIs for CO2 (counted multiple times if DOI for multiple materials) %d  \n", length(doi_list_CO2) )
end