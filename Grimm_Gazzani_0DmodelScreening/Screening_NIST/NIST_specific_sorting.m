function NIST_specific_sorting(name_data_nist_document, name_char_nist_document,NIST_TP_SORTED)
%========================================================
% Here, materials can be deleted by hand, if needed. Furhermore, 
% the structure is changed.
%-------------------------------------------------------
%Input:   - name_data_nist_document The name of the document in whcih the
%           isotherm data is stored.
%         - name_char_nist_document The name of the document in which the
%           characteristics data is stored
%         - NIST_TP_SORTED Name of the new document, with another structure
%           and where some isotherms are deleted.
%Output:  - File with the name NIST_TP_SORTED
%========================================================
fprintf("≡ NIST specific sorting  ≡\n")

%--- read the structure from NIST ---
data_nist_open  = fopen(name_data_nist_document, 'r');
structure_from_NIST = readstruct(name_data_nist_document);
fclose (data_nist_open);

%--- Sort NIST specifically ---
%  %NUM 8
%  structure_from_NIST.Carbon_Dioxide.Carbon.ozdvzow___b__k__c__s_dwzofdhvdwdhss  = rmfield(structure_from_NIST.Carbon_Dioxide.Carbon.ozdvzow___b__k__c__s_dwzofdhvdwdhss  , "TaAa196" ) ;
%  structure_from_NIST.Carbon_Dioxide.Carbon.ozdvzow___b__k__c__s_dwzofdhvdwdhss   = rmfield(structure_from_NIST.Carbon_Dioxide.Carbon.ozdvzow___b__k__c__s_dwzofdhvdwdhss  , "TaAa196num6457" );

 
%--- change order in structure ---
% structure_from_NIST has a structure Gas-material-doi-temperature
% but in order to eventually put it back in the fit_nist_3 program, it
% needs: material-gas-temperature


%first, reduce Gas-material-doi-temperature to Gas-material-temperature
%where the DOI is saved in every structure with the name Doi_NIST
count_dois = 0;
material_list = fieldnames(structure_from_NIST.Carbon_Dioxide);

char_nist_open  = fopen(name_char_nist_document, 'r');
char_struct_from_NIST = readstruct(name_char_nist_document);
fclose (char_nist_open);

newtest_count_debug = 0;
newtest_count_debug_2 = 0;
count_structures_inconsistent= 0;
count_doi_deleted_char_NaN = 0;
for every_material =  [1:length(material_list)]
    material_name = string(material_list(every_material));
    structure_from_NIST_new.(material_name).Full_name = material_name; %because original data structure from NiST also contains full name, and leaving the full name out gives problems using analyze_data  
    doi_list = fieldnames(structure_from_NIST.Carbon_Dioxide.(material_name));
    for every_doi = [1:length(doi_list)]
        newtest_count_debug = newtest_count_debug + 1;
        doi_name = string(doi_list(every_doi));
        try
            char_list = fieldnames(char_struct_from_NIST.Carbon_Dioxide.(material_name).(doi_name));
            flag_char_exist = 1;
            count_strings = 0;
            for every_char = [1:length(char_list)]
                characteristic = string(char_list(every_char));
                value_characteristic = char_struct_from_NIST.Carbon_Dioxide.(material_name).(doi_name).(characteristic);
                if isstring(value_characteristic) == 1 %class(value_characteristic) == string %""
                    count_strings = count_strings + 1; %only the model characteristic is a string. If more characteristics are strimngs, it means that some values are NaN: disgard doi
                    if count_strings >= 3 % the model and full name are always strings. if any other char is also a string, something went wrong
                        flag_char_exist = 0;
                    end
                end
            end
            newtest_count_debug_2= newtest_count_debug_2 +1;
            if flag_char_exist == 1 % if all characteristics are numbers
                count_dois = count_dois + 1;
                temperature_list = fieldnames(structure_from_NIST.Carbon_Dioxide.(material_name).(doi_name));
                for every_temp = [1:length(temperature_list)]
                    temperature_name = string(temperature_list(every_temp));
                    structure_from_NIST_new.(material_name).Carbon_Dioxide.(temperature_name) = structure_from_NIST.Carbon_Dioxide.(material_name).(doi_name).(temperature_name);
                    structure_from_NIST_new.(material_name).Carbon_Dioxide.(temperature_name).Doi_NIST = doi_name;
                end
            else
                 count_doi_deleted_char_NaN = count_doi_deleted_char_NaN + 1;
            end
        catch
            fprintf("Material name in list of TP but not in list char")
            count_structures_inconsistent = count_structures_inconsistent + 1 ; % needed because of the wc < 10e14 line
        end
    end
end

fprintf("DOIs counted: %d \n", count_dois)
fprintf("DOIs deleted because a character is NaN: %d \n", count_doi_deleted_char_NaN)
fprintf("DOIs deleted because working capacity > 10^14: %d \n", count_structures_inconsistent)

Nist_TP_Data_sorted_open = fopen(NIST_TP_SORTED, 'wt');
writestruct(structure_from_NIST_new,'NIST_TP_DATA_SORTED.xml');                    
fclose(Nist_TP_Data_sorted_open);