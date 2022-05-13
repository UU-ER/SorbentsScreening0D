%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% Utrecht University, Netherlands 
% Copernicus Institute of Sustainable Development, Section of Energy & Resources 
%
% Year: 	2022
% MATLAB: 	R2021a
% Authors:	Alexa Grimm (AG) (group of Dr. Matteo Gazzani) 
% 
% Purpose: 
% The 0D model simulates a Temperature Vacuum Swing Adsorption (TVSA) cycle 
% for CO2 capture from dilute streams, proposed in <Title of paper, doi>. 
% This model is designed to screen large database, e.g. the NIST/APRA-E database,
% and rank the resulting materials. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% =====================================
% run the 0D model and calculate the process performance parameters

%% INPUT DATA
% here, different materials can be added
    currentSorbent = 1; % 1: exemplary
    isothermModel = "toth_cp";
    % get data (example)
    sorbent_data_screen;
    
    Tdes = 373; % desoprtion temperature, K
    pvac = 0.27/10; % vacuum pressure, MPa
    V_feed = 8.9*1e-06; % m3/s
    yCO2 = 4e-4;
    yH2O = 0.0134;
    feed_concentration = [yCO2,1-yCO2-yH2O,yH2O];

    % full saturation or saturation level calculated with neural network
    full_saturation = false;

%% RUN MODEL
out = main_0D_model(dataModel,isothermModel,Tdes,pvac,V_feed,feed_concentration,full_saturation); 






