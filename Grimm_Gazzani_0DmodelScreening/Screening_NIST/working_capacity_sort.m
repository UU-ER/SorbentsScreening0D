function dataset = working_capacity_sort (name_datafile,pCO2_1,sorted_wc_structure, specific_material_list)
%========================================================
% Opens name_datafile, sorts on positive working capacity, saves structure
% in specific_material_list and calls plots.m
%-------------------------------------------------------
%Input:   - name_datafile String name of file in which database to be sorted is saved
%         - pCO2_1 Partial pressure of CO2,input needed for
%           plots.m
%         - sorted_wc_structure  String name of file in which resulting combined
%           database is saved
%         - specific_material_list Specific materials, input needed for
%           plots.m
%Output:  - Multiple plots (by calling plots.m)
%         - File named sorted_wc_structure in which the sorted structure is
%           stored.
%========================================================
fprintf("≡ working_capacity_sort  ≡\n")

data_open = fopen(name_datafile, 'r');
data = readstruct(name_datafile); 
fclose (data_open);

count_dois_fitted = 0;
dH_list = []; wc_list =[]; gamma_list=[]; m0_list= []; q0_list=[]; ns_equ_list = []; KH_list=[]; list_of_R2=[]; material_list_final =[]; p =[];
KH_save = 0; ns_equ_save = 0; q0_save = 0; m0_save = 0; gamma_save = 0; wc_save = 0;dH_save = 0;dH_list = 0; R2 = 0; model = "";
model_list =[];
total_doi_list = [];
Pr_list = []; Q_list = []; purity_list = []; recovery_list = [];

count_dois_fitted_WC_neg = 0;
dH_list_WC_neg = []; wc_list_WC_neg =[]; gamma_list_WC_neg=[]; m0_list_WC_neg= []; q0_list_WC_neg=[]; ns_equ_list_WC_neg = []; KH_list_WC_neg=[]; list_of_R2_WC_neg=[]; material_list_final_WC_neg =[]; p_WC_neg =[];
KH_save_WC_neg = 0; ns_equ_save_WC_neg = 0; q0_save_WC_neg = 0; m0_save_WC_neg = 0; gamma_save_WC_neg = 0; wc_save_WC_neg = 0;dH_save_WC_neg = 0;dH_list_WC_neg = 0; R2_WC_neg = 0; model_WC_neg = "";
model_list_WC_neg =[];

 p_list = sprintf("concentration_%.0f",pCO2_1);

material_list = fieldnames(data.Carbon_Dioxide);
    for every_material =  [1:length(material_list)]
        material_name = string(material_list(every_material));
        doi_list = fieldnames(data.Carbon_Dioxide.(material_name));
        for every_doi = [1:length(doi_list)]
            doi_name = string(doi_list(every_doi));
            total_doi_list =[total_doi_list doi_name];
            characteristics_list = fieldnames(data.Carbon_Dioxide.(material_name).(doi_name));
            for every_char = [1:length(characteristics_list)]
                characteristic_name = string(characteristics_list(every_char));
                if characteristic_name == "wc_save"
                    specific_working_capacity = data.Carbon_Dioxide.(material_name).(doi_name).wc_save;
                    if specific_working_capacity >= 0.05
                        material_name
                        doi_name
                        count_dois_fitted = count_dois_fitted +1;
                        dataset.new_structure.Carbon_Dioxide.(material_name).(doi_name) = data.Carbon_Dioxide.(material_name).(doi_name);
                        material_full_name = data.Carbon_Dioxide.(material_name).(doi_name).Full_name;
                        [dH_list, wc_list, gamma_list, m0_list, q0_list, ns_equ_list, ...
                            KH_list, dataset, list_of_R2, material_list_final, p, model, full_name ,T0,rho,rhoMat,cp] ...
                            = assign_characteristics ("new_structure", dataset, "Carbon_Dioxide", material_name, material_full_name, ...
                            doi_name, count_dois_fitted, KH_save, ns_equ_save , q0_save, m0_save, ...
                            gamma_save, wc_save,dH_save,dH_list, wc_list,gamma_list, ...
                            m0_list, q0_list, ns_equ_list,KH_list,R2,list_of_R2,1,...
                            material_list_final,p,model,[],[],[],[]);
                        model_list = [model_list, model];
                        
                        % run 0D model ------------------------------------
                        addpath('./0D_model/');
                        addpath('C:\Program Files\NAG\MB25\mbw6i25ddl');
                        % flag, E_spec, Q_spec, Pr, purity, recovery        
                        out = run_0D_model_screening2(material_name,doi_name,data,pCO2_1);
                        % save calculated process performance:
                            dataset.new_structure.Carbon_Dioxide.(material_name).(doi_name).(p_list).Pr = out(4);
                            dataset.new_structure.Carbon_Dioxide.(material_name).(doi_name).(p_list).Qth = out(3);
                            dataset.new_structure.Carbon_Dioxide.(material_name).(doi_name).(p_list).purity = out(5);
                            dataset.new_structure.Carbon_Dioxide.(material_name).(doi_name).(p_list).recovery = out(6);

                        % performance positive WC
                        Pr_list = [Pr_list out(4)]; 
                        Q_list = [Q_list out(3)]; 
                        purity_list = [purity_list out(5)]; 
                        recovery_list = [recovery_list out(6)];
                        
                    else
                        % new: (I also want the characteristics for the other
                        % materials with WC<0
                        count_dois_fitted_WC_neg = count_dois_fitted_WC_neg +1;
                        dataset.new_structure.Carbon_Dioxide.(material_name).(doi_name) = data.Carbon_Dioxide.(material_name).(doi_name);
                        material_full_name = data.Carbon_Dioxide.(material_name).(doi_name).Full_name;
                        [dH_list_WC_neg, wc_list_WC_neg,gamma_list_WC_neg, m0_list_WC_neg, q0_list_WC_neg, ns_equ_list_WC_neg, ...
                            KH_list_WC_neg, dataset, list_of_R2_WC_neg,material_list_final_WC_neg,p_WC_neg,model_WC_neg,full_name_WC_neg,T0,rho,rhoMat,cp]...
                            = assign_characteristics ("new_structure", dataset, "Carbon_Dioxide", material_name, material_full_name, ...
                            doi_name, count_dois_fitted_WC_neg, KH_save_WC_neg, ns_equ_save_WC_neg , q0_save_WC_neg, m0_save_WC_neg, ...
                            gamma_save_WC_neg, wc_save_WC_neg, dH_save_WC_neg, dH_list_WC_neg, wc_list_WC_neg, gamma_list_WC_neg, ...
                            m0_list_WC_neg, q0_list_WC_neg, ns_equ_list_WC_neg, KH_list_WC_neg, R2_WC_neg, list_of_R2_WC_neg, 1, ...
                            material_list_final_WC_neg, p_WC_neg, model_WC_neg,[],[],[],[]);
                        model_list_WC_neg = [model_list_WC_neg, model_WC_neg];
                        
                        % empty process performance:
                        dataset.new_structure.Carbon_Dioxide.(material_name).(doi_name).(p_list).Pr = NaN;
                        dataset.new_structure.Carbon_Dioxide.(material_name).(doi_name).(p_list).Qth = NaN;
                        dataset.new_structure.Carbon_Dioxide.(material_name).(doi_name).(p_list).purity = NaN;
                        dataset.new_structure.Carbon_Dioxide.(material_name).(doi_name).(p_list).recovery = NaN;
                    end
                end
            end
        end
    end
    
fprintf("Total number of DOIs: %d \n", length(total_doi_list))    
fprintf("Number of DOIs with a positive working capacity: %d \n", count_dois_fitted)
fprintf("Number of DOIs with a negative working capacity: %d \n", count_dois_fitted_WC_neg)

charac_WC_pos = [dH_list; wc_list; gamma_list; m0_list; q0_list; ns_equ_list; ...
    KH_list; list_of_R2; material_list_final]';
charac_WC_neg = [dH_list_WC_neg; wc_list_WC_neg; gamma_list_WC_neg; m0_list_WC_neg; q0_list_WC_neg; ns_equ_list_WC_neg;...
    KH_list_WC_neg; list_of_R2_WC_neg; material_list_final_WC_neg]';
performance_WC_pos = [Pr_list;Q_list;purity_list;recovery_list;material_list_final];

if count_dois_fitted ~= 0
%     plots("R2_yes", "characteristics_yes", dH_list, wc_list ,gamma_list, m0_list, q0_list, ns_equ_list, KH_list, list_of_R2,pCO2_1,material_list_final,specific_material_list )
    plotType = [true, true, true, true]; % distribution R2, characteristics all, characteristics positive, performance 
    plots2(plotType, charac_WC_pos, charac_WC_neg, pCO2_1, specific_material_list, performance_WC_pos)
                        
    amount_of_dois = length(wc_list) ;
    max_length = length(wc_list);
    datalist = sortrows([wc_list; [1:max_length]]',1, 'descend')';
    model_list = model_list(datalist(2,:));
    for i = [1:amount_of_dois]
        %fprintf (" %s & %d & %f & %s\n", material_list_final(datalist(2,i)) , datalist(1,i), list_of_R2(datalist(2,i)), model_list(i) )
        fprintf (" %s & %.3f & %.3f & %s\n", material_list_final(datalist(2,i)) , round(datalist(1,i),3), round(list_of_R2(datalist(2,i)),3), model_list(i) );
    end
    fprintf("Total number of points: %d \n",length(wc_list))
end

existing_models = categories( categorical(model_list));
model_number = countcats( categorical (model_list));
for i = [1:length(model_number)]
    fprintf(" %s  %d times \n", existing_models{i}, model_number(i) )
end

gamma_list = gamma_list(datalist(2,:));
gamma_list_toth =[];gamma_list_s =[];gamma_list_DSL = [];
for i =[1:length(model_list)]
    if model_list(i) =="toth_cp"
        gamma_list_toth = [gamma_list_toth gamma_list(i)];
    elseif model_list(i) =="DSL"
        gamma_list_DSL = [gamma_list_DSL gamma_list(i)];
    else
        gamma_list_s = [gamma_list_s gamma_list(i)];
    end
end

fprintf("Mean value non-linearity for Toth-cp: %.3f \n", round(mean(gamma_list_toth),3));
fprintf("The maximum value of the non-linearity for Totch-cp:%.3f \n ", round(max(gamma_list_toth),3) );
fprintf("Mean value non-linearity for S-shape: %.3f \n",round(mean(gamma_list_s),3) );
fprintf("The maximum value of the non-linearity for S-shape:%.3f \n ", round(max(gamma_list_s),3) );
fprintf("Mean value non-linearity for DSL: %.3f \n", round(mean(gamma_list_DSL),3));
fprintf("The maximum value of the non-linearity for DSL:%.3f \n ", round(max(gamma_list_DSL),3) );
fprintf("Mean value non-linearity for all models: %.3f \n", round(mean(gamma_list),3));

new_sorted_structure_open = fopen(sorted_wc_structure, 'wt');
writestruct(dataset.new_structure,sorted_wc_structure);                    
fclose(new_sorted_structure_open);
end