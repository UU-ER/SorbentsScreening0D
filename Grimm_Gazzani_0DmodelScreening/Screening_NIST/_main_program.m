%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% Utrecht University, Netherlands 
% Copernicus Institute of Sustainable Development, Section of Energy & Resources 
%
% Year: 	2022
% MATLAB: 	R2021a
% Authors:	Alexa Grimm (AG) (group of Dr. Matteo Gazzani) 
% 
%
% Purpose: 
% With this code the isotherm parameters of the materials provided in the 
% online NIST/APRA-E database can be read and saved. 
% In a next step, these materials are filtered, sorted, fitted and specific 
% isotherm metrics are calculated. In addition, a 0D model can be run for the 
% resulting materials to calculate the process performance parameters, i.e. 
% productivity, specific energy, purity and recovery. More details can be 
% found in <Title of paper, doi>.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MAIN PROGRAMM
% The main program has the following functions:
%        1) Read out the NIST/APRA-E database
%        2) Filtering the data from the NIST/APRA-E database and the data 
%           by Khurana and Farooq(2016); computing and
%           displaying their characteristics.
%        3) Combining the databases, filtering on positive working capacity
%        4) Computing the energy consumption for specific materials
%        5) Plotting isotherms for specific materials in one graph

%% 1) Reading the NIST/APRA-E database
% Data from the NIST/APRA_E database is read out and saved in a file name
% specified by name_to_save
% NOTE: Running this sub-section can take long! (+- 6 hours)

scan_NIST = false;

if scan_NIST
    initial_data = "database_100_DOIs.xml"; % choose a database name: e.g. database_100_DOIs.xml
    i_max = 1;                            % if i_max = 1, save all DOIs (+-6 hrs!)
    nist_ccdc_4 (i_max,initial_data );  
end

%% 2) Filtering the data from the NIST/APRA-E database and the data by Khurana and Farooq (2016)
% The NIST/APRA-E data is read out from the file name_to_save. The
% characteristis are computed for the pressures specified by
% adsorption_pressures.
clear all
initial_data  = "NIST_database.xml"; 
adsorption_pressures = 40; 

for every_pressure = adsorption_pressures
    % The characteristics of the data by Khurana and Farooq are computed and
    % saved in a filename specified by string_name_KhuFar
    string_name_KhuFar = strcat("KHUFAR_CHAR_", string(every_pressure),".xml")
    Khu_Far (every_pressure,string_name_KhuFar) %compute and plot characteristics
       
    % The database NIST/APRA-E is filtered in fit_nist_3. Some data for 
    % which the error is large is adapted in NIST_specific_sorting.
    % Then, the data is fitted again. The characteristics are stored in the
    % file string_name_NIST_CHAR_2.
    string_name_NIST_CHAR_1  = strcat("NIST_CHAR_DATA1_",string(every_pressure),".xml");
    string_name_NIST_TP_1   = strcat("NIST_TP_DATA1_",string(every_pressure),".xml");
    fit_nist_3 (initial_data,string_name_NIST_CHAR_1,string_name_NIST_TP_1,"yes","yes",every_pressure)
    delete NIST_TP_DATA_SORTED.xml
    NIST_specific_sorting(string_name_NIST_TP_1, string_name_NIST_CHAR_1,'NIST_TP_DATA_SORTED.xml')
    
    %% FITTING of all the sorbents ----------------------------------------
    string_name_NIST_CHAR_2 = strcat("NIST_CHAR_DATA2_", string(every_pressure),".xml");
    string_name_NIST_TP_2 = strcat("NIST_TP_DATA2_", string(every_pressure),".xml");
    fit_nist_3 ('NIST_TP_DATA_SORTED.xml',string_name_NIST_CHAR_2,string_name_NIST_TP_2,"no","yes",every_pressure)
end

%% 3) Combining the databases, filtering on positive working capacity
% The two databases obtained above are combined and filtered on a 
% positive working capacity.

dataset = [];
every_pressure = 40;
string_name_KhuFar = strcat("KHUFAR_CHAR_", string(every_pressure),".xml");
string_name_NIST_CHAR_2 = strcat("NIST_CHAR_DATA2_", string(every_pressure),".xml");
string_combined_CHAR = strcat("data_", string(every_pressure), ".xml");
string_sorted_file = strcat("data_wc_sorted_", string(every_pressure), ".xml");
combine_databases(string_name_KhuFar, string_name_NIST_CHAR_2,string_combined_CHAR); % name of dataset = string_combined_CHAR
for every_temperature = [100, 150] % celsius
    for every_pressure = [40, 100, 1000, 10000] % Pa
        dataset = Analyze_all_sorbents(string_combined_CHAR,every_pressure,every_temperature,dataset);
    end
end
% save dataset
save('data_var_yCO2_Tdes','dataset')

