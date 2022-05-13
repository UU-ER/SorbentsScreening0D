function fit_nist_3 (name_data_nist_document,structure_name_characteristics,structure_name_TP,first_sorting,second_sorting,pCO2_1)
%========================================================
% Process NIST & CCDC - with right structure implementation
% Gets data out of the NIST database, saves data from relevant isotherms
% (CO2) in a file (useful_data_nist.xml). CCDC data can be manually added.
%-------------------------------------------------------
%Input:     - name_data_nist_document, document containing all data
%           - structure_name_characteristics Document in which
%             characteristics will be saved
%           - structure_name_TP Document in which the isotrherms (only
%              Temperature, Pressure and amount adsorbed) will be saved
%           - first_sorting If "yes" carries out first sorting round
%           - second_sorting If "yes" carries out second sorting round
%           - pCO2_1 Partial pressure of adsorption
% Program:  Opens the structure saved in name_data_nist_document
%           Sorts and re-organizes it.  The data in the structures is 
%           fitted, characteristics are computed and plotted.
% Output:  - File named structure_name_characteristics in which a structure
%            with all chracteristics are stored.
%          - File named structure_name_TP in which a structure
%            in which all temperature/pressure/amount adsorbed data is stored.
%          - Multiple plots (see plots function)
%========================================================
fprintf("≡ fit_nist_3 ≡\n")

%---SETTINGS---
plot_isotherms = "no";  %plotting all isotherms
flag_with_versions = 0; % no versions
list_convertible_units = ["mmol/g","mmol/kg","mol/g","kg/mol", "cm3(STP)/g"]; %when changed, change check_units.m accordingly

%--initialize counter ---
Counter_m_step_1_2 = 0;

%--- read the structure from NIST ---
data_nist_open = fopen(name_data_nist_document, 'r');
data_NIST = readstruct(name_data_nist_document);
fclose (data_nist_open);

%--- analyze original NIST data ---
[gas_display, row_length_table_NIST_1, table_number_of_T, material_name,table_mat_temps] = analyze_data (data_NIST);

%% --- data for which there is only one temperature available are deleted ---
if first_sorting =="yes"
    %make structure sorted_data_NIST and delete all entries that are not useful (with only 1 Temperature)
    sorted_data_NIST = data_NIST; 
    for every_gas = [1:length(gas_display)]
        for material_number = [1:row_length_table_NIST_1]
            if table_number_of_T(material_number,every_gas ) <= 1 %if there is only 1 unique temperature available
                if isfield(sorted_data_NIST, (material_name(material_number))) == 1 %if the field exists
                    if isfield(sorted_data_NIST.(material_name(material_number)),(gas_display(every_gas))) ==1
                        temp_to_be_deleted = string(fieldnames(sorted_data_NIST.(material_name(material_number)).(gas_display(every_gas)))); 
                        sorted_data_NIST.(material_name(material_number)).(gas_display(every_gas)) = rmfield(sorted_data_NIST.(material_name(material_number)).(gas_display(every_gas)), (temp_to_be_deleted)); 
                       if numel(fieldnames(sorted_data_NIST.(material_name(material_number)).(gas_display(every_gas)))) == 0 %if structure under gas in empty, delete the gas
                            sorted_data_NIST.(material_name(material_number)) = rmfield(sorted_data_NIST.(material_name(material_number)) , (gas_display(every_gas)) );  
                       end
                       if numel(fieldnames(sorted_data_NIST.(material_name(material_number)))) == 1 %if struct under materials only consists of "full name" and has no gas, delete it
                            sorted_data_NIST = rmfield(sorted_data_NIST , (material_name(material_number)) );
                       end
                    end
                end
            end
        end 
    end
    %--- analyze newly sorted data ---
    [gas_display_sorted, row_length_table_NIST_1_sorted, table_number_of_T_sorted, material_name_sorted,table_mat_temps_sorted] = analyze_data (sorted_data_NIST);
end

%% structuring based on gases, sort on doi, fitting
if second_sorting == "yes"
    %--- set counters ---
    %counters for units
    count_cm3stpg = 0;          count_wt =0;        count_mlg =0;
    count_mmolg =0;             count_molmol =0;    count_mmolkg =0;
    count_molecules_unitcell =0;count_mgg =0;       count_gg =0;
    count_cm3STPcm3 =0;         count_mmolm2 =0;    count_molg =0;
    count_molm3 =0;             count_gcm3 =0;      count_gL =0;
    count_moleculespore =0;     count_mmolcm3 =0;   count_kgmol= 0;
    count_mgm2 = 0;             count_moleculeseightunitcells =0;
    count_bar = 0;              count_Pa =0;        count_pascal =0;

    %counters sorting
    Counter_m_step_2_3 = 0;         Counter_d_step_2_3 =0;  Counter_d_step_3_4 = 0;
    Counter_d_step_4_5 = 0;         Counter_d_step_5_6 = 0; Counter_d_step_6_7 = 0;
    Counter_d_step_8_9_13_14 = 0;	count_dois_toth =  0;   count_dois_s_shape = 0;
    count_dois_langmfr = 0;
    count_plots = 0;                count_isotherms_plotted = 0;
    count_size_combos_strings_larger = 0; Counter_d_step_8_13 =0;

    %--- set empty lists ---
    dH_list = [];   wc_list =[];        gamma_list=[]; 
    m0_list= [];    q0_list=[];         ns_equ_list = []; 
    KH_list=[];     list_of_R2=[];      material_names =[];
    list_all_units_amount_adsorbed = [];list_of_R2 = [];

    clearvars structure_from_NIST
    
    %--- initial data-related values are assigned ---
    if first_sorting =="yes"
        size_table_mat_temps_sorted = size(table_mat_temps_sorted);
        max_Ts = size_table_mat_temps_sorted (3);
        max_material_number = row_length_table_NIST_1_sorted;
    else
        size_table_mat_temps_sorted = size(table_mat_temps);
        max_Ts = size_table_mat_temps_sorted (3);
        max_material_number = row_length_table_NIST_1;
        table_mat_temps_sorted = table_mat_temps;
        sorted_data_NIST = data_NIST;
    end
    
    %--- make a new structure, based on the gasses ---
    dataset.fitted_data = struct();
    for every_gas = 1 %Only CO2
        if first_sorting =="yes"
            every_gas_name =  gas_display_sorted(every_gas); 
        else
            every_gas_name =  gas_display(every_gas); 
        end    
        
        %% new save bad fitting
        no_R2_low = 0;

       for material_number = [1:max_material_number]
           clearvars temperature_name_for_input fitting_input doi_of_T non_unique_dois
           
           % --- assign material field name for structure (struc_name) and full name ---
           if first_sorting =="yes"
               material_struc_name = material_name_sorted (material_number);
           else
               material_list = fieldnames(data_NIST);
               material_struc_name = string(material_list(material_number));
           end
           material_full_name = data_NIST.(material_struc_name).Full_name;
           
            %--- initialize flags, counters and lists ---
            flag_T_present = 0;
            flag_doi_not_unique = 0; 
            count_temp_doi = 0; 
            T_list = []; 

            %--- determine number of temperatures for a material and gas --- 
            for Ts = [1: max_Ts]  
                % table contains if a certain Temperatures is available for a certain material & gas
                place_in_table = table_mat_temps_sorted(every_gas,material_number,Ts ); 
                if isempty(place_in_table{1,1}) == 0 %if this place in the table is not empty (0)
                    if startsWith(string(place_in_table),"TaAa") ==1
                         T = string(table_mat_temps_sorted(every_gas,material_number,Ts ));

                         % save DOI to check later
                         count_temp_doi = count_temp_doi + 1;
                         doi_of_T(count_temp_doi) = sorted_data_NIST.(material_struc_name).(every_gas_name).(T).Doi_NIST;

                         %this gas-material combination exists, temperature is added to list
                         T_list = [T_list T];
                         flag_T_present = 1; 
                    end
                end
            end

            T_list_original = T_list; 

            %-- checking if doi is not unique (at least 2 datasets (temperatures) per DOI) ---
            non_unique_dois =[];
            if flag_T_present == 1  %if gas-material combination exists
                flag_non_unique_doi = 0;
                Counter_m_step_1_2 = Counter_m_step_1_2 + 1;
                for every_doi_of_T = [1:length(doi_of_T)]
                    count_dois_int = sum(ismember (doi_of_T, doi_of_T(every_doi_of_T)),'all');
                    if count_dois_int > 1
                        %only fitting the materials for which there is a non-unique doi
                        flag_non_unique_doi =1; % at least one not-unique
                        non_unique_dois = [non_unique_dois doi_of_T(every_doi_of_T)];
                    end
                end
            end

            %--- filter the data and fit the data ---
            if  flag_non_unique_doi == 1 && flag_T_present == 1 %if gas-material combination exists and at least one DOI has two datasets available
                non_unique_dois = unique(non_unique_dois);
                Counter_m_step_2_3 = Counter_m_step_2_3 +1;
                
                for every_non_unique_doi = [1:length(non_unique_dois)]
                    material_version = string(non_unique_dois(every_non_unique_doi));
                    Counter_d_step_2_3 = Counter_d_step_2_3 +1;

                    %initialize for fitting and filtering temperature lists and doi lists
                    clearvars fitting_input
                    T_list_2 =[];       T_list_3 = [];      doi_list =[];   number_list = [];
                    if first_sorting =="yes"
                        material_version_name = doi_convert_to_name(material_version);
                    else
                        material_version_name = material_version;    
                    end
                    T_list = T_list_original ;
                    list_unit_amount_adsorbed = []; 
                    flag_fitted_because_of_variants = 0; % becomes 1 if versions are made (see make_versions)

                    %--- get the right temperature list per DOI ---
                    sort_on_doi

                    %--- fit the data ---
                    % if no versions are made and if there are two temperatures to fit
                    if flag_fitted_because_of_variants ==0 && flag_version_has_only_one_T == 0
                        fitting
                    end
                end
            end
        end
    end
end
%%
    list_all_units_amount_adsorbed = unique(list_all_units_amount_adsorbed);

%--- save data in files ---
try
    Nist_CHAR_Data_open = fopen(structure_name_characteristics, 'wt');
    dataset_name = "fitted_data";
    writestruct(dataset.(dataset_name),structure_name_characteristics);                    
    fclose(Nist_CHAR_Data_open);  

    Nist_TP_Data_open = fopen(structure_name_TP, 'wt');
    writestruct(structure_from_NIST,structure_name_TP);                    
    fclose(Nist_TP_Data_open);

    fprintf("Data succesfully saved. \n")
catch
    fprintf("Saving data in structures was not succesful. \n")
end

%--- print counters ---
fprintf("=========")
fprintf("Counter_m_step_1_2: %d \n", Counter_m_step_1_2)  
fprintf("Counter_m_step_2_3: %d \n", Counter_m_step_2_3)
fprintf("Counter_d_step_2_3: %d \n", Counter_d_step_2_3)
fprintf("Counter_d_step_3_4: %d \n", Counter_d_step_3_4)
fprintf("Counter_d_step_4_5: %d \n", Counter_d_step_4_5)
fprintf("Counter_d_step_5_6: %d \n", Counter_d_step_5_6)
fprintf("Counter_d_step_6_7: %d \n", Counter_d_step_6_7)
fprintf("Counter_d_step_8_13: %d \n", Counter_d_step_8_13)
fprintf("Counter_d_step_8_9_13_14: %d \n", Counter_d_step_8_9_13_14)
fprintf("Fitted with Toth: %d \n", count_dois_toth)
fprintf("Fitted with Langmuir-Freundlich: %d \n", count_dois_langmfr)
fprintf("Fitted with S-shape : %d \n", count_dois_s_shape)
fprintf("Added Toth+s_shape+langm-fr: %d \n", (count_dois_toth+count_dois_s_shape+count_dois_langmfr))
fprintf("#of plots: %d \n", count_plots)
fprintf("mean R2: %f \n", mean(list_of_R2))
fprintf("Isothems plotted (only non-zero if plot_isotherms == 'yes'): %d \n", count_isotherms_plotted)
%--- make plots ---
plots('R2_yes', 'characteristics_yes',dH_list, wc_list ,gamma_list, m0_list, q0_list, ns_equ_list, KH_list, list_of_R2,pCO2_1,[],[])

save('save_R2_low','save_R2_low')
end
