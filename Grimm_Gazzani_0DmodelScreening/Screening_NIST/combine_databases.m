function combine_databases(database_1, database_2,new_database)
%========================================================
% Combines database_1 and database_2 and saves the result in a file called
% new_database.
%-------------------------------------------------------
%Input:   - database_1 String name of file in which database 1 is saved
%         - database_2 String name of file in which database 2 is saved
%         - new_database String name of file in which resulting combined
%           database is saved
%Output:  - File with the name new_database that contains combined database
%========================================================

fprintf("≡ combine_databases  ≡\n")

d1_open  = fopen(database_1, 'r');
data_1 = readstruct(database_1);
fclose (d1_open);

d2_open  = fopen(database_2, 'r');
data_2 = readstruct(database_2);
fclose (d2_open);


%--- combining database process ---
material_list_1 = fieldnames(data_1.Carbon_Dioxide);
material_list_2 = fieldnames(data_2.Carbon_Dioxide);
material_list = unique([material_list_1; material_list_2]);
doi_list_for_count = [];

count_dois_1 = 0;
for every_mat_1 = [1:length(material_list_1)]
    mat = string(material_list_1(every_mat_1));
    count_dois_1 = count_dois_1+length(fieldnames(data_1.Carbon_Dioxide.(mat)));
end
fprintf("dois 1: %d",count_dois_1)
count_dois_2 = 0;
for every_mat_2 = [1:length(material_list_2)]
    mat = string(material_list_2(every_mat_2));
    count_dois_2 = count_dois_2+length(fieldnames(data_2.Carbon_Dioxide.(mat)));
end
fprintf("dois 2: %d",count_dois_2)

for every_mat = [1:length(material_list)] %for every material
    mat = string(material_list(every_mat));
    doi_list_for_count_one_mat  = [];
    if isfield(data_1.Carbon_Dioxide, (mat)) == 1 && isfield(data_2.Carbon_Dioxide, (mat))== 0 %only in data1
        new_structure.Carbon_Dioxide.(mat) = data_1.Carbon_Dioxide.(mat);
        doi_list_1 = fieldnames(data_1.Carbon_Dioxide.(mat));
        for every_doi = [1 : length(doi_list_1)]
            doi = string(doi_list_1(every_doi));
            doi_list_for_count_one_mat = [doi_list_for_count_one_mat, doi];
        end
    elseif isfield(data_1.Carbon_Dioxide, (mat)) == 0 && isfield(data_2.Carbon_Dioxide, (mat)) == 1 %only in data2
        new_structure.Carbon_Dioxide.(mat) = data_2.Carbon_Dioxide.(mat);
        doi_list_2 = fieldnames(data_2.Carbon_Dioxide.(mat));
        for every_doi = [1 : length(doi_list_2)]
            doi = string(doi_list_2(every_doi));
            doi_list_for_count_one_mat = [doi_list_for_count_one_mat, doi];
        end
    else % both in data_1 & data_2!
        doi_list_1 = fieldnames(data_1.Carbon_Dioxide.(mat));
        doi_list_2 = fieldnames(data_2.Carbon_Dioxide.(mat));
        doi_list_combined = [doi_list_1; doi_list_2];
        for every_doi = [1 : length(doi_list_combined)]
            doi = string(doi_list_combined(every_doi));
            doi_list_for_count_one_mat = [doi_list_for_count_one_mat, doi];
        end
        doi_list = unique(doi_list_combined);
        for every_doi = [1 : length(doi_list)]
            doi = string(doi_list(every_doi));
           
            if isfield(data_1.Carbon_Dioxide.(mat),(doi)) == 1 && isfield(data_2.Carbon_Dioxide.(mat),(doi)) == 0
                new_structure.Carbon_Dioxide.(mat).(doi) = data_1.Carbon_Dioxide.(mat).(doi);
            elseif isfield(data_1.Carbon_Dioxide.(mat),(doi)) == 0 && isfield(data_2.Carbon_Dioxide.(mat),(doi)) == 1
                new_structure.Carbon_Dioxide.(mat).(doi) = data_2.Carbon_Dioxide.(mat).(doi);
            else
                T_list_1 = fieldnames(data_1.Carbon_Dioxide.(mat).(doi));
                T_list_2 = fieldnames(data_2.Carbon_Dioxide.(mat).(doi));
                T_list   = unique([T_list_1; T_list_2]);
                for every_T = [1:length(T_list)]
                    temp = string(T_list(every_T));
                    if isfield(data_1.Carbon_Dioxide.(mat).(doi),(temp)) == 1 && isfield(data_2.Carbon_Dioxide.(mat).(doi),(temp)) == 0
                        new_structure.Carbon_Dioxide.(mat).(doi).(temp) = data_1.Carbon_Dioxide.(mat).(doi).(temp);
                    elseif isfield(data_1.Carbon_Dioxide.(mat).(doi),(temp))== 0 && isfield(data_2.Carbon_Dioxide.(mat).(doi),(temp)) == 1
                         new_structure.Carbon_Dioxide.(mat).(doi).(temp) = data_2.Carbon_Dioxide.(mat).(doi).(temp);
                    else
                        error Error merging databases: some materials have the same temperatures!
                    end
                end
            end
        end
    end
    doi_list_for_count = [doi_list_for_count doi_list_for_count_one_mat];
end
 
fprintf("\n %s has %d materials \n", database_1,length(fieldnames(data_1.Carbon_Dioxide)));
fprintf("%s has %d materials \n", database_2,length(fieldnames(data_2.Carbon_Dioxide)));
fprintf("Combined, there are %d DOIs \n", length(doi_list_for_count) );

%--- save new database ---
new_combined_database = fopen(new_database, 'wt');
writestruct(new_structure,new_database);                    
fclose(new_combined_database);

end