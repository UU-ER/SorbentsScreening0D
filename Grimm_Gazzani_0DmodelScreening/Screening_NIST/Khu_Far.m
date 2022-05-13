function Khu_Far (pCO2_1,structure_name_characteristics)
%=====================================
% Compute characteristics for the 74 materials analysed by Khurana, M., 
% & Farooq, S. (2016) https://doi.org/10.1021/acs.iecr.5b04531
%--------------------------------------
%Input:  - pCO2_1 Adsorption pressure
%        - structure_name_characteristics Name of file in which data will
%          be saved
%        - The structure 'materialProperties.mat', which contains the
%           parameters for the materials
%Output: - File named structure_name_characteristics in which structure 
%          with characteristics of the materials is saved. 
%        - Multiple plots (see plots function)
%--------------------------------------
fprintf("≡ Khurana & Farooq  ≡\n")
load('materialProperties.mat'); % qs1 (mol/m3), b0 1/(mol/m3), dU1 (kJ/mol), qs2 (mol/m3), d0 1/(mol/m3), dU2 (kJ/mol),

%initialize variables
clearvars    dH_list wc_list  gamma_list m0_list q0_list ns_equ_list KH_list material_names 
dH_list = []; wc_list =[]; gamma_list=[]; m0_list= []; q0_list=[]; ns_equ_list = []; KH_list=[]; list_of_R2=[]; material_list = [];
dataset_name = "fitted_data_Khu_Far";
material_version_name = doi_convert_to_name('10.1021/acs.iecr.5b04531');
model = 'DSL';
number_of_materials_KF = 0;
number_of_materials_KF_wc_check = 0;
% qs1 (mol/m3 -> mol/kg), b0 (m3/mol), dH1 (J/mol), qs2 (mol/m3 -> mol/kg), d0 (m3/mol), dH2 (J/mol)
p_CO2 = materialProperties.isothermParametersCO2;
material_number_name_list = materialProperties.materialList;
material_density_list = materialProperties.materialDensity;
every_gas_name = "Carbon_Dioxide";

for every_material = [1:length(p_CO2)]
    %obtain the parameters, parameter 1 and 4 need to be divided by the
    %materials' density
    number_of_materials_KF = number_of_materials_KF + 1;
    material_number_name_normal = string(material_number_name_list(every_material));
    density_material = material_density_list(every_material);
    for every_parameter = [1: length(p_CO2(every_material,:))] % save the fitted parameters and convert the units
        p(every_parameter) = p_CO2(every_material,every_parameter);
        if every_parameter ==1 
            p(every_parameter) = p(every_parameter) / density_material; % mol/m3 -> mol/kg
        elseif every_parameter == 4
            p(every_parameter) = p(every_parameter) / density_material; % mol/m3 -> mol/kg
        end
    end
    
    % add material properties
    rhoMat = density_material;
    void_fraction = 0.35;
    rho = density_material*(1-void_fraction);
    cp = 1070; 
    T0 = [];
    
    % compute the characteristics
   [KH_save, ns_equ_save, q0_save, m0_save, gamma_save, wc_save, dH_save, T0] = characteristics_4(p,model,pCO2_1,rhoMat,T0,120) ;
   %check if all characteristics have a plausible value (no NaNs)
   if class(wc_save) ~= string
       number_of_materials_KF_wc_check = number_of_materials_KF_wc_check+1;
       material_number_name = convert_string_name (string(material_number_name_list(every_material)),"material");
       fitted_data_Khu_Far.(every_gas_name).(material_number_name).(material_version_name).p = p;
       if number_of_materials_KF == 1
            dataset.fitted_data_Khu_Far = fitted_data_Khu_Far;
       end
       %save characteristics in data list
       [dH_list, wc_list, gamma_list, m0_list, q0_list, ns_equ_list, KH_list, dataset,list_of_R2,material_list,p,model, full_name ,T0,rho,rhoMat,cp] = assign_characteristics ("fitted_data_Khu_Far", dataset, every_gas_name, material_number_name, material_number_name_normal, material_version_name, every_material,  KH_save, ns_equ_save , q0_save, m0_save, gamma_save, wc_save,dH_save,dH_list, wc_list,gamma_list, m0_list, q0_list, ns_equ_list,KH_list, 1, list_of_R2,0,material_list,p,model,T0,rho,rhoMat,cp);
   end
end

%save structure
Nist_CHAR_Data_open = fopen(structure_name_characteristics, 'wt');
writestruct(dataset.(dataset_name),structure_name_characteristics);                    
fclose(Nist_CHAR_Data_open);  

% make plots
plots('R2_no', 'characteristics_yes',dH_list, wc_list ,gamma_list, m0_list, q0_list, ns_equ_list, KH_list, [],pCO2_1, [],[] )
end