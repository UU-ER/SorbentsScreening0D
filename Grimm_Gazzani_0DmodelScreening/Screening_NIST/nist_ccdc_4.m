function nist_ccdc_4 (i_max, name_to_save)
%========================================================
% Gets data out of the NIST/APRA-E database,saves data from relevant isotherms
% (CO2) in a file (name_to_save). CCDC data can be manually added.
%-------------------------------------------------------
%Input:     - i_max number of data_documents scanned. If i_max =1, all
%             documents that can be founcd are scanned.
%           - name_to_save name of the .xml document in which data is
%             saved.
%Program:   Opens all data in biography that inlcudes daya on CO2 adsorbents,
%           opens the ones that are adsorbents for CO2. 
%            The data is being saved in a structure called materials.
%           In case include_CCDC == yes: the T, density and formula of the
%           CCDC database can be manually stored as well.
% Output:  - The file bibliography_nist.csv which contains a bibliography of 
%            all academic works that measured isotherm of at least one CO2 adsorbent
%          - The file url_isotherms.txt which contains urls of all data
%            from the bibliography
%          - The file name_to_save which contains all the data as a
%            structure
% NOTES: - Change deleting files preference: otherwise bin on computer can 
%          be overflown by documents as the program makes a document for 
%          each isotherm(which is about 16200) and deletes them all again.
%          How to change deleting files preferences: Go to Home tab,
%          click preferences, select MATLAB >general. Choose 'Delete
%          permanently'
%========================================================

fprintf("-------nist_ccdc_4------- \n");

%---SETTINGS PART 1/2 ---
process_nist_data_from_website = "yes"; %can be deleted
take_already_existing_bib = "no"; % if file  bibliography_nist.csv exists, change to "yes"
take_already_existing_url_isotherms = "no"; % if file  url_isotherms.txt exists, change to "yes" 
include_CCDC = "no" ; % if yes: allows for manually inserting CCDC data. If no: only takes NIST data

bibliography = 'bibliography_nist.csv';
url_nist = 'https://adsorbents.nist.gov/isodb/';
              options = weboptions('RequestMethod','get','ArrayFormat', ...
              'csv','ContentType','text','Timeout',120); 
url_ccdc ='https://www.ccdc.cam.ac.uk/structures/Search?Doi=';
                  
%--- save bibliography of NIST ----
if take_already_existing_bib == "no"
    gas = 'CO2';
    api = sprintf('api/biblio-gas/%s.json',gas);
    url_data_nist = strcat(url_nist,api); 
    outfilename = websave(bibliography,url_data_nist,options);
end

if process_nist_data_from_website == "yes"
    %--- open NIST bibliography ---
    if take_already_existing_url_isotherms == "no"
        bibliography_open = fopen(bibliography,'rt');
            if bibliography_open < 0
                error('error opening file %s\n\n',bibliography)
            end
        %--- get and save urls NIST ---
        txt = readtable(bibliography,'HeaderLines',0,'ReadVariableNames',false);
        txt_height = height(txt);
        TF_filename = startsWith(txt.Var1,"filename");

        count_filenames = 0;
        isotherm_open = fopen('url_isotherms.txt', 'wt');
        for i = [1: txt_height]
          if TF_filename(i) == 1
              count_filenames = count_filenames+1;
              string_txt = string(txt.Var1(i));
              cell_erased_filename = erase( string_txt,"filename: ");
              cell_erased_quotations = strip( cell_erased_filename,"""");
              filenames_txt{count_filenames} = cell_erased_quotations;
              api = sprintf('api/isotherm/%s.json', filenames_txt{count_filenames});

              url{count_filenames,1} = strcat(url_nist,api);
              fprintf(isotherm_open, '%s ;', url{count_filenames});
          end
        end

        fclose(isotherm_open);
        fclose(bibliography_open);
    end

    
%--- set up counters ---
count_number_of_materials_w_same_name = 0;
count_deleted_isotherms = 0;
count_useful_data_NIST = 0;
count_useful_data_NIST_and_CCDC = 0;
count_number_of_same_isotherms = 0;

%--- open file with isotherm urls and read ---

url_isotherms_open= fopen('url_isotherms.txt', 'r')    ;
string_isotherm = fgetl(url_isotherms_open);
fclose(url_isotherms_open);

isotherms_urls = split(string_isotherm, ';');

%--- SETTINGS PART 2/2 ---
total_number_of_data = length(isotherms_urls);
fprintf("%d isotherms were taken from the NIST Databank. \n",total_number_of_data-1 ); %-1 because last data point is ''
if i_max == 1;
  i_max =  length(isotherms_urls) ; % about 16200 isotherms 
end


%--- start of reading data of every isotherm url -------
pressure_name_nist = "pressure_NIST";
amount_adsorbed_name_nist = "amount_adsorbed_NIST";
url_name_nist = "URL_NIST";
doi_name_nist = "Doi_NIST";
isotherm_doc_name_nist = "Isotherm_doc";
full_name = "Full_name";
unit_pressure_name = "Unit_pressure";
unit_adsorption_name = "Unit_adsorption";

count_data_skipped = 0;
data_is_useful = 1;
        for i = [1:i_max] 
           if ismember(i, [0:100:16000]) ==1
              fprintf("\n +iteration: %d \n",i)
           end
           data_is_useful = 1;
           clearvars txt_isotherm
           while data_is_useful == 1  %this loop is broken if data is not CO2 specific/not present in CCDC data
              %for every isotherm url, a document is made and read
              isotherm_name_csv = sprintf('isotherm_%d.csv', i);
              url = char(isotherms_urls(i));
              try % to avoid errors
                  outfilename = websave(isotherm_name_csv,url,options);
              catch
                  try
                    outfilename = websave(isotherm_name_csv,url,options);
                  catch 
                     fprintf("Isotherm %d. skipped because protocol specified in URL, '', is not supported. \n", i )
                    count_data_skipped = count_data_skipped+1;
                    data_is_useful = 0;
                  end
              end
              open_specific_isotherm = fopen(isotherm_name_csv,'r');
              if open_specific_isotherm < 0
                 fprintf("Isotherm %d. Isotherm skipped because .csv file could not be openend. \n", i )
                 count_data_skipped = count_data_skipped+1;
                 data_is_useful = 0;
                 break
              end
            
              txt_isotherm = readtable(isotherm_name_csv,'HeaderLines',0,'ReadVariableNames',false);

               [size_row, size_column] = size(txt_isotherm);
               if size_row < 12 % I changed this, before it was 8 but I got an error with a material containing only 11 rows
                   fprintf("Isotherm %d. Isotherm skipped because sizerow (%d) is smaller than 8. \n", i, size_row )
                   count_data_skipped = count_data_skipped+1;
                   data_is_useful = 0;   
                   break
               end
              %TF arrays are made in order to see where the data is stored in the files.
              TF_pressure = startsWith(txt_isotherm.Var1,"pressure");
              TF_temperature = startsWith(txt_isotherm.Var1,"temperature:");
              TF_name = startsWith(txt_isotherm.Var1,"name:");
              TF_doi = startsWith(txt_isotherm.Var1,"DOI");
              TF_unit_adsorption = startsWith(txt_isotherm.Var1,"adsorptionUnits");
              
              %initialize counters
              time_names = 0; % there are two data entries with name, first one is name of material, second one name of the gas adsorbed
              count_useful_data_NIST = count_useful_data_NIST + 1;
              deleted_file = 0; 
              count_number_of_pressures = 0;
              
              for j = [1:height(txt_isotherm)] %for all lines in isotherm file
                 % check the temperature
                 if TF_temperature (j) ==1
                     temperature_nist = str2num(strip(string (txt_isotherm.Var2(j)), ',') );
                 end  
                 
                 % check the name (one for gas, one of material)
                 if TF_name(j) == 1
                   time_names = time_names +1;  
                   if time_names == 2 %check second name
                      gas = txt_isotherm.Var2(j);
                   elseif time_names ==1 
                      name_material_final = txt_isotherm.Var2(j);
                      string_name_material = convert_string_name ( name_material_final, 'name' );
                   end
                 end

                 % check the doi
                 if TF_doi (j) ==1
                     doi = txt_isotherm.Var2(j);
                     doi = strip( doi, ',');
                     doi_store{i} = doi;
                 end
                 
                 %check unit of adsorption
                 if TF_unit_adsorption (j) ==1
                     unit_adsorption = string(txt_isotherm.Var2(j));
                     unit_adsorption = strip( unit_adsorption, ',');
                 end

                 %check the pressure
                 if TF_pressure(j) == 1  % if line starts with pressure
                      count_number_of_pressures = count_number_of_pressures +1;
                      if count_number_of_pressures ==1
                          unit_pressure = strip(string (txt_isotherm.Var2(j)), ',');
                      else
                          pressure_to_compare(count_number_of_pressures - 1) = str2num(strip(string (txt_isotherm.Var2(j)), ',') );      
                          %adsorbed amount is usually on the next line: if
                          %this gives a problem then take adsorbed_amount =0.00
                          try
                            adsorbed_amount(count_number_of_pressures - 1) = str2num(strip(string (txt_isotherm.Var2(j+1)), ',') );
                          catch
                            adsorbed_amount(count_number_of_pressures - 1) = 0.00;
                          end
                      end

                 end
              end % end for all lines
              
              % check if pressure vector is there
              if exist('pressure_to_compare') == 0
                  data_is_useful = 0;
                  fprintf("The pressure vector is empty.")
                  count_data_skipped = count_data_skipped+1;
              end
              % store data in materials structure
                  if include_CCDC == "yes" 
                      %get url to CCDC
                      url_doi = string(strcat(doi, '&DatabaseToSearch=Published'));
                      url_data_ccdc = strcat(url_ccdc,url_doi); 
                      sprintf("Go to: %s",url_data_ccdc )
                      %manually give CCDC input
                      is_your_input_ok = 0; %while input not ok (0) repeat
                      while  is_your_input_ok == 0
                          what_to_do = input("y =continue, n = not found. ", 's');
                          what_to_do = string(what_to_do);
                          if what_to_do == 'n'
                              materials.(string_name_material);
                              materials = rmfield(materials, string_name_material);
                              data_is_useful = 0;
                              is_your_input_ok = 1;
                          elseif what_to_do == 'y'
                              count_useful_data_NIST_and_CCDC  = count_useful_data_NIST_and_CCDC+1;
                              formula_ccdc = input("What is the formula? Formula: ", 's');
                              temperature_ccdc = input("What is the temperature? Temperature (K): ");
                              density_ccdc = input("What is the density? Density: ") ;
                              your_input_ok = strcat("Check your input, type 1 if it is right, 0 if wrong. " );%, your_input_total, ". ");
                              is_your_input_ok = input(your_input_ok);  % in case of mistake in entering data
                          end
                      end
                  end
                 if data_is_useful ==1 
                     %check if the exact same isotherm already exists (that is same material, same gas, same temperature)
                     temperature_name_nist = convert_string_name (temperature_nist, 'temperature');
                     temperature_name_nist = strcat ("T",temperature_name_nist);
                     gas_name_nist = string(convert_string_name (gas, 'gas'));
                     if count_useful_data_NIST >1
                         if isfield(materials,(string_name_material)) ==1 && isfield(materials.(string_name_material),(gas_name_nist))== 1 && isfield(materials.(string_name_material).(gas_name_nist),(temperature_name_nist))== 1;
                             count_number_of_same_isotherms = count_number_of_same_isotherms +1;
                             number_behind_temperature = sprintf ( 'num%d', count_number_of_same_isotherms);
                             temperature_name_nist = strcat (temperature_name_nist, number_behind_temperature);
                         end
                     end

                     % debug: check whether number of actual data points matces
                     % number of data points that this program found.
                     number_of_pressures_acc_to_TF =sum(TF_pressure);
                     if number_of_pressures_acc_to_TF ~= count_number_of_pressures
                         fprint ("Error: number of data points given in NIST data does not match number of datapoints that this program sees")
                         break
                     end
                      materials.(string_name_material).(full_name) = string_name_material; %'Pressure NIST', pressure_to_compare);
                      materials.(string_name_material).(gas_name_nist).(temperature_name_nist).(pressure_name_nist) = pressure_to_compare; %'Pressure NIST', pressure_to_compare);
                      materials.(string_name_material).(gas_name_nist).(temperature_name_nist).(amount_adsorbed_name_nist) = adsorbed_amount;
                      materials.(string_name_material).(gas_name_nist).(temperature_name_nist).(url_name_nist) = url;
                      materials.(string_name_material).(gas_name_nist).(temperature_name_nist).(doi_name_nist) = doi;
                      materials.(string_name_material).(gas_name_nist).(temperature_name_nist).(isotherm_doc_name_nist) = isotherm_name_csv;
                      materials.(string_name_material).(gas_name_nist).(temperature_name_nist).(unit_pressure_name) = unit_pressure;
                      materials.(string_name_material).(gas_name_nist).(temperature_name_nist).(unit_adsorption_name) = unit_adsorption;
                      if include_CCDC == "yes" 
                          formula_name_ccdc = 'Formula_CCDC';
                          temperature_name_ccdc = 'Temperature_CCDC';
                          denisty_name_ccdc = 'Density_CCDC';
                          url_name_ccdc = 'URL_CCDC';
                          materials.(string_name_material).(formula_name_ccdc) =  formula_ccdc ;
                          materials.(string_name_material).(temperature_name_ccdc) =  temperature_ccdc ;
                          materials.(string_name_material).(denisty_name_ccdc) =  density_ccdc ;  
                          materials.(string_name_material).(url_name_ccdc) =  url_data_ccdc ;
                      end
                 end
              clearvars  TF_pressure  TF_doi  TF_name  TF_temperature TF_unit_adsorption count_number_of_pressures pressure_to_compare adsorbed_amount
              data_is_useful = 0;
           end  
           if deleted_file == 0;  % delete the file if it has not been deleted yet
               fclose('all');
               delete(isotherm_name_csv);     
           end
        end

%---save structure in .xml file ---   
    try
        useful_data_nist_open = fopen(name_to_save, 'wt'); %'useful_data_nist_2.xml'
        writestruct(materials,name_to_save);                    
        fclose(useful_data_nist_open);
        fprintf("Data succesfully saved in %s. \n", name_to_save)
    catch
        fprintf("Saving data in %s was not succesful. \n", name_to_save)
    end
 end
close all
end