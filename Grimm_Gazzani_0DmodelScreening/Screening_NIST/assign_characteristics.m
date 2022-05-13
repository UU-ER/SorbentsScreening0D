 
function [dH_list, wc_list,gamma_list, m0_list, q0_list, ns_equ_list,KH_list, dataset, list_of_R2,material_list,p, model,full_name,T0,rho,rhoMat,cp] = assign_characteristics (dataset_name, dataset, gas, material, full_name,version, count_dois, KH_save, ns_equ_save , q0_save, m0_save, gamma_save, wc_save,dH_save,dH_list, wc_list,gamma_list, m0_list, q0_list, ns_equ_list,KH_list,R2,list_of_R2, structure_exists, material_list,p,model,T0,rho,rhoMat,cp)
%========================================================,rho
% Depending on the value of structure_exists, either obtains data from
% structure or puts data in structure. Also makes lists of the
% characteristics, error and material names.
%-------------------------------------------------------
%Input:     - dataset_name Name of the structure that the data is stored in
%           - dataset Structure in which data is stored
%           - gas name of gas (e.g. 'carbon_Dioxide") 
%           - material The field name material (e.g. ZIF_da_8)
%           - full_name The full, orginal material name (e.g. ZIF-8)
%           - version Name DOI
%           - count_dois Count used to fill the lists correctly
%           - KH_save Value Henry's coeffecient
%          - ns_equ_save Value equilibrium loading
%          - q0_save Value amount adsorbed at feed concentration
%          - m0_save Value slope isotherm at feed concentration
%          - gamma_save Non-linearity
%          - wc_save Working capacity
%          - dH_save Heat of adsorption
%          - dH_list, wc_list,gamma_list, m0_list, q0_list,
%               ns_equ_list,KH_list Lists of the characteristics mentioned above
%          - R2 Value of coefficient of determination
%          - list_of_R2,
%          - structure_exists If 1: data is taken from structure. 
%                             If 0: structure is filled with data. 
%          - material_list List of all materials including DOI
%          - p Parameters from fitting
%          - model Name of model used for fitting
% Output:  -  dH_list, wc_list,gamma_list, m0_list, q0_list,
%               ns_equ_list,KH_list Lists of the characteristics mentioned
%               above to which values have been added
%          - dataset Structure in which characteristics are stored / from
%            which characteristics were obtained.
%          - list_of_R2 list of all values of coefficient of determination 
%          - material_list List of material names including DOI 
%          - p List of parameters used for fitting
%          - model Name of model used for fitting 
%          - full_name The full name of the material. 
%========================================================

dataset_name = string(dataset_name);
if structure_exists == 0
    dataset.(dataset_name).(gas).(material).(version).KH_save       = KH_save;
    dataset.(dataset_name).(gas).(material).(version).ns_equ_save   = ns_equ_save;
    dataset.(dataset_name).(gas).(material).(version).q0_save       = q0_save;
    dataset.(dataset_name).(gas).(material).(version).m0_save       = m0_save;
    dataset.(dataset_name).(gas).(material).(version).gamma_save    = gamma_save;
    dataset.(dataset_name).(gas).(material).(version).wc_save       = wc_save;
    dataset.(dataset_name).(gas).(material).(version).dH_save       = dH_save;
    dataset.(dataset_name).(gas).(material).(version).error         = R2;
    dataset.(dataset_name).(gas).(material).(version).p             = p; % [1/Pa]
    dataset.(dataset_name).(gas).(material).(version).model         = model;
    dataset.(dataset_name).(gas).(material).(version).Full_name     = full_name;
    % added new:
    if isempty(rhoMat)
        rhoMat = 1130; % kg/m3
    end
    dataset.(dataset_name).(gas).(material).(version).rhoMat           = rhoMat;
    dataset.(dataset_name).(gas).(material).(version).T0           = T0;
    
    if isempty(rho)
        rho = 1130*(1-0.35); % kg/m3
    end
    dataset.(dataset_name).(gas).(material).(version).rho           = rho;
    
    if isempty(cp)
        cp = 1070; % J/kg/K
    end
    dataset.(dataset_name).(gas).(material).(version).cp           = cp;
    
    
    
else
    KH_save     = dataset.(dataset_name).(gas).(material).(version).KH_save;
    ns_equ_save = dataset.(dataset_name).(gas).(material).(version).ns_equ_save  ;
    q0_save     = dataset.(dataset_name).(gas).(material).(version).q0_save  ;
    m0_save     = dataset.(dataset_name).(gas).(material).(version).m0_save;
    gamma_save  = dataset.(dataset_name).(gas).(material).(version).gamma_save ;
    wc_save     = dataset.(dataset_name).(gas).(material).(version).wc_save ;
    dH_save     = dataset.(dataset_name).(gas).(material).(version).dH_save  ; 
    R2          = dataset.(dataset_name).(gas).(material).(version).error ;
    p           = dataset.(dataset_name).(gas).(material).(version).p;
    model       = dataset.(dataset_name).(gas).(material).(version).model;
    full_name   = dataset.(dataset_name).(gas).(material).(version).Full_name;
    % added new:
    if isempty(dataset.(dataset_name).(gas).(material).(version).rhoMat)
        rhoMat = 1130; % kg/m3
        dataset.(dataset_name).(gas).(material).(version).rho = rhoMat;
    else
        rhoMat     = dataset.(dataset_name).(gas).(material).(version).rhoMat;
    end
    T0          = dataset.(dataset_name).(gas).(material).(version).T0;
    
    if isempty(dataset.(dataset_name).(gas).(material).(version).rho)
        rho = 1130*(1-0.35); % kg/m3
        dataset.(dataset_name).(gas).(material).(version).rho = rho;
    else
        rho     = dataset.(dataset_name).(gas).(material).(version).rho;
    end
    
    if isempty(dataset.(dataset_name).(gas).(material).(version).cp)
        cp = 1070; % J/kg/K
        dataset.(dataset_name).(gas).(material).(version).rho = cp;
    else
        cp     = dataset.(dataset_name).(gas).(material).(version).cp;
    end
end

dH_list (count_dois)        = dH_save;
wc_list (count_dois)        = wc_save;
gamma_list (count_dois)     = gamma_save;
m0_list (count_dois)        = m0_save;
q0_list (count_dois)        = q0_save;
ns_equ_list (count_dois)    = ns_equ_save;
KH_list (count_dois)        = KH_save;
list_of_R2 (count_dois)     = R2;

vers_name = string(convert_string_name_to_normal ( version, "doi"));
first_cat = strcat(full_name, "  ") ;
second_cat = strcat(first_cat , vers_name );
material_list  = [material_list, string(second_cat)];
end