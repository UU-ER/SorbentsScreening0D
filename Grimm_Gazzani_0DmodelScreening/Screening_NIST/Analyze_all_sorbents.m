function dataset = Analyze_all_sorbents (name_datafile,pCO2_1,Tdes,dataset)
%=====================================
% Compute characteristics for all sorbents
% they already have characteristic data, but here these parameters are
% calculated again for different input parameters
fprintf("≡ Analyzing all sorbents  ≡\n")

p_list = sprintf("concentration_%.0f",pCO2_1);
T_list = sprintf("desorption_temperature_%.0f",Tdes);

name_analysis = strcat("data_", string(pCO2_1), "_Pa_",string(Tdes),"_C");

% calculate characteristics
data_open = fopen(name_datafile, 'r');
dataset.(name_analysis) = readstruct(name_datafile); 
fclose (data_open);

% calculated characteristics
characteristics_list_short = {'KH_save','ns_equ_save','q0_save','m0_save','gamma_save','wc_save','dH_save'};

count_dois_fitted = 0;
dH_list = []; wc_list =[]; gamma_list=[]; m0_list= []; q0_list=[]; ns_equ_list = []; KH_list=[]; list_of_R2=[]; material_list_final =[]; p =[];
KH_save = 0; ns_equ_save = 0; q0_save = 0; m0_save = 0; gamma_save = 0; wc_save = 0;dH_save = 0;dH_list = 0; R2 = 0; model = "";
model_list =[];
total_doi_list = [];
Pr_list = []; Q_list = []; purity_list = []; recovery_list = [];
material_names =[];

count_dois_fitted_WC_neg = 0;
dH_list_WC_neg = []; wc_list_WC_neg =[]; gamma_list_WC_neg=[]; m0_list_WC_neg= []; q0_list_WC_neg=[]; ns_equ_list_WC_neg = []; KH_list_WC_neg=[]; list_of_R2_WC_neg=[]; material_list_final_WC_neg =[]; p_WC_neg =[];
KH_save_WC_neg = 0; ns_equ_save_WC_neg = 0; q0_save_WC_neg = 0; m0_save_WC_neg = 0; gamma_save_WC_neg = 0; wc_save_WC_neg = 0;dH_save_WC_neg = 0;dH_list_WC_neg = 0; R2_WC_neg = 0; model_WC_neg = "";
model_list_WC_neg =[];
number_of_materials_OS = 0;
Counter_d_step = 0;

every_gas_name = "Carbon_Dioxide";
%%
material_list = fieldnames(dataset.(name_analysis).Carbon_Dioxide);
    for every_material =  [1:length(material_list)]
        material_name = string(material_list(every_material));
        material_number_name_normal = convert_string_name_to_normal ( material_name, "name");
        doi_list = fieldnames(dataset.(name_analysis).Carbon_Dioxide.(material_name));
        for every_doi = [1:length(doi_list)]
            doi_name = string(doi_list(every_doi));

            % obtain the parameters, parameter 1 and 4 need to be divided by the
            % materials' density
            rhoMat = dataset.(name_analysis).(every_gas_name).(material_name).(doi_name).rhoMat;
            rho = dataset.(name_analysis).(every_gas_name).(material_name).(doi_name).rho;
            cp = dataset.(name_analysis).(every_gas_name).(material_name).(doi_name).cp;

            p = dataset.(name_analysis).(every_gas_name).(material_name).(doi_name).p;
            model = string(dataset.(name_analysis).(every_gas_name).(material_name).(doi_name).model);
            T0 = dataset.(name_analysis).(every_gas_name).(material_name).(doi_name).T0;

            % compute the characteristics for concentration and temperature
           [KH_save, ns_equ_save, q0_save, m0_save, gamma_save, wc_save, dH_save] = characteristics_4(p,model,pCO2_1,rhoMat,T0,Tdes) ;
           
            Counter_d_step = Counter_d_step + 1; 
            structure_exists = 0;
            %--- clear existing characteristics ---
            for every_char = [1:length(characteristics_list_short)]
                characteristic_name = string(characteristics_list_short(every_char));
                dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).(characteristic_name) = [];
            end
            Counter_d_step = 1;
            %--- save the characteristics ---     
            [dH_list, wc_list,gamma_list, m0_list, q0_list, ns_equ_list,KH_list,dataset,list_of_R2,material_names,p,model,material_number_name_normal,T0 ,rho,rhoMat,cp] = assign_characteristics (name_analysis, dataset, every_gas_name, material_name, material_number_name_normal ,doi_name, Counter_d_step, KH_save, ns_equ_save , q0_save, m0_save, gamma_save, wc_save,dH_save,dH_list, wc_list,gamma_list, m0_list, q0_list, ns_equ_list,KH_list,R2,list_of_R2,structure_exists,material_names,p,model,T0,rho,rhoMat,cp);

            for every_char = [1:length(characteristics_list_short)]
                characteristic_name = string(characteristics_list_short(every_char));
                if characteristic_name == "wc_save"
                    specific_working_capacity = dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).wc_save;
                    if specific_working_capacity >= 0.05
                        count_dois_fitted = count_dois_fitted +1;
                        
                        % run 0D model ------------------------------------
                        % here the 0D model needs to be added
                        addpath('./0D_model/');
                        out = run_0D_model_screening2(material_name,doi_name,dataset,pCO2_1,Tdes,name_analysis);
                        % save calculated process performance:
                            dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).(p_list).Pr = out(4);
                            dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).(p_list).Qth = out(3);
                            dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).(p_list).purity = out(5);
                            dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).(p_list).recovery = out(6);
                            
                            % TEST
                            if pCO2_1>50 && out(6)<0.90
                                % empty process performance:
                                dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).(p_list).Pr = NaN;
                                dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).(p_list).Qth = NaN;
                                dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).(p_list).purity = NaN;
                                dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).(p_list).recovery = NaN;
                            elseif out(6)>1.01
                                % empty process performance:
                                dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).(p_list).Pr = NaN;
                                dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).(p_list).Qth = NaN;
                                dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).(p_list).purity = NaN;
                                dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).(p_list).recovery = NaN;
                            else
                                % performance positive WC
                                Pr_list = [Pr_list out(4)]; 
                                Q_list = [Q_list out(3)]; 
                                purity_list = [purity_list out(5)]; 
                                recovery_list = [recovery_list out(6)];
                            end
                    else
                        % new: (I also want the characteristics for the other
                        % materials with WC<0
                        count_dois_fitted_WC_neg = count_dois_fitted_WC_neg +1;
                        
                        % empty process performance:
                        dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).(p_list).Pr = NaN;
                        dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).(p_list).Qth = NaN;
                        dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).(p_list).purity = NaN;
                        dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).(p_list).recovery = NaN;
                    end
                end
            end
        end
    end
end