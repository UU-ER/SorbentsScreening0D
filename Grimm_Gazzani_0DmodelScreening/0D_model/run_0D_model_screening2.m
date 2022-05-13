function out = run_0D_model_screening2(material_name,doi_name,dataset,pCO2_1,Tdes,name_analysis)
%=====================================
% run the 0D model and calculate the process performance parameters
% this file is created to use the screening data as input
% if no screening data is available, use main_0D.m as main file
    
    isothermModel = dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).model;
    p = dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).p;
    
    handles = dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name);
    if isfield(handles,'T0')
        T0 = dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).T0;
    else
        sprintf('no temperature assigned!')
    end

    % check if water isotherm is available. If not, the same water isotherm
    % for all materials is added
    if isfield(handles,'modelH2O')
        own_H2O = true;
        pH2O = dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).modelH2O.pH2O;
    else
        own_H2O = false;
    end
    
    % assign the isotherm parameters and physical properties   
    sorbent_data_screen2;
    
    %% DATA INPUT
    Tdes = Tdes+273; % K
    pvac = 0.15*1e5*1e-6; % MPa
    V_feed = 8.9*1e-06; % m3/s
    yCO2 = pCO2_1*1e-5;
    yH2O = 0.0134;
    feed_concentration = [yCO2,1-yCO2-yH2O,yH2O];
    
    % full saturation or saturation level predicted by a neural network
    full_saturation = false;
    
    if isothermModel == 'langm_freund'
        isothermModel = 'langfr';
    end
    
    if exist('dataModel','var')
    else
        sprintf('not there')
    end
    
    % out = [flag, E_spec, Q_spec, Pr, purity, recovery]
    out = main_0D_model(dataModel,isothermModel,Tdes,pvac,V_feed,feed_concentration,full_saturation); 
    
end