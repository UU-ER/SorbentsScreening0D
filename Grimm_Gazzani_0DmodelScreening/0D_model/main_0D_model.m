function out = main_0D_model(dataModel,isothermModel,Tdes,pvac,V_feed,feed_concentration,full_saturation)
%=====================================
% run the 0D model and calculate the process performance parameters

dataModel.currentSorbent = isothermModel; 

% GENERAL DATA
dataModel.general.gasconstant = 8.314; % universal gas constant, J/mol K
dataModel.general.adiabaticConstant = 1.4;
dataModel.general.pumpEfficiency = 0.72;
dataModel.general.airBlowerEfficiency = 0.72;
dataModel.general.densityAir = 1.23; % kg/m3
dataModel.general.MMCO2 = 44; % Molar mass CO2, g/mol
dataModel.general.MMN2 = 14; % g/mol
dataModel.general.MMH2O = 18; % g/mol

% PROCESS SPECIFIC DATA
dataModel.process.Tamb = 293; % ambient temperature 
dataModel.process.Tdes = Tdes; % desoprtion temperature, needs to be higher than 
dataModel.process.pamb = 1.001*1e5*1e-6; % total pressure, MPa
dataModel.process.pvac = pvac; % vacuum pressure, MPa

dataModel.process.adsorbentMass = 1; % Mass of the adsorbent, kg
dataModel.process.voidFraction = (dataModel.sorbent.MaterialDensity- dataModel.sorbent.Density)/dataModel.sorbent.MaterialDensity; % column void fraction
dataModel.process.heatResistance = 16.8; % global heat transfer resistance, (J/(m2 K))
dataModel.process.heightFrame = 0.05; % height of sorbent containing frame, (m)
dataModel.process.area = dataModel.process.adsorbentMass/...
    dataModel.sorbent.Density/dataModel.process.heightFrame/2; % column cross section
dataModel.process.noSteps = 100; % number of steps per unit
dataModel.process.airVelocity = V_feed/(pi*(0.005/2)^2); % air velocity 
dataModel.process.Vol = dataModel.process.adsorbentMass/(dataModel.sorbent.Density); % volume column 

% FEED STREAM: N2, H2O, CO2
    dataModel.feed.yCO2 = feed_concentration(1);
    dataModel.feed.yN2 = feed_concentration(2);
    dataModel.feed.yH2O = feed_concentration(3);
    
% Saturation degree and starting concentration BD
    if full_saturation
        yCO2_sat = dataModel.feed.yCO2;
    else
        load('netSat','netNN_sat');
        rho = dataModel.sorbent.Density;
        sat = netNN_sat([rho;V_feed;Tdes;pvac*10]);
        if sat > 1
            sat=1;
        end
        dataModel.degreeSaturation = sat;
        % calculate CO2 concentration at saturation level
        yCO2_sat = calcSat(dataModel);
    end
    dataModel.startingCondBD = [yCO2_sat,(1-yCO2_sat-dataModel.feed.yH2O),dataModel.feed.yH2O];

%% PRESSURE AND TEMPEARTURE VECTOR
[pressureVector, timeBD, TimeBDProfilePlot] = pressureProfile(dataModel,dataModel.process.pvac);
dataModel.pressureVector = pressureVector; % MPa
dataModel.process.noStepsBD = round(timeBD+2);
dataModel.process.TimeBDProfilePlot = TimeBDProfilePlot; % 

% old (linear temperature profile):
[THeatProfileStep, timeHeating, TimeHeatProfilePlot] = temperatureProfile(dataModel,dataModel.process.Tdes);
dataModel.THeatProfileStep = THeatProfileStep;
dataModel.process.noStepsHeating = round(timeHeating+1);
dataModel.process.TimeHeatProfilePlot = TimeHeatProfilePlot; % 

dataModel.process.coolingTime = 350;
[TCoolProfileStep, TimeCoolProfilePlot, pressureVectorCool] = coolingProfile_1D(dataModel,dataModel.process.Tdes); % 
dataModel.TCoolProfileStepPlot = TCoolProfileStep;
dataModel.TCoolProfileStep = TCoolProfileStep;
dataModel.PressureCoolProfileStep = pressureVectorCool;
dataModel.process.TimeCoolProfilePlot = TimeCoolProfilePlot; % 


%% ADSORPTION FEED 
dataModel.general.MMFeed = dataModel.feed.yCO2*dataModel.general.MMCO2+dataModel.feed.yN2*dataModel.general.MMN2+dataModel.feed.yH2O*dataModel.general.MMH2O; % g/mol

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% RUN THE MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%% %%
try
    %------------------------ blowdown and evacuation ------------------------%
    % outputBlowEvac = [yCO2_save, yN2_save, yH2O_save, Nout_save_sum, Nout_save_CO2, Nout_save_N2, Nout_save_H2O, E_total, T_save, Q_sum_save];
    outputBlowEvac = simulateBlowdownEvacuation(dataModel);

    %-------------------------------- cooling --------------------------------%
    % outputCool = [yCO2_save, yN2_save, yH2O_save, Nin_save, Nin_save_CO2(2:end), Nin_save_N2(2:end), Nin_save_H2O(2:end),T_save, Q_sum_save, Ntotal_CO2, Ntotal_N2, Ntotal_H2O];
    outputCooling = simulateCooling(dataModel,outputBlowEvac);

    %------------------------------ adsorption -------------------------------%
    %  outAdsorption = [yCO2_save, yN2_save, yH2O_save, Nwaste_save_sum, Nwaste_save_CO2, Nwaste_save_N2, Nwaste_save_H2O, E_total (J), N_feed];
    outputAdsorption = simulateAdsorption(dataModel,outputCooling);

    %% results   
    mfeed = sum(outputAdsorption(:,9))*dataModel.general.MMFeed/1000; % kg
    t_ads = mfeed/dataModel.process.airVelocity/dataModel.process.area/dataModel.general.densityAir; % s

    QTotal = outputBlowEvac(end,10); % kJ
    ETotal = (outputBlowEvac(end,8)+sum(outputAdsorption(:,8)))/1000; % kJ
    N_CO2_product = outputBlowEvac(end,5)-outputBlowEvac(dataModel.process.noStepsBD,5); % mol

    %% performance
    % electricity
        E_spec = ETotal/(N_CO2_product*dataModel.general.MMCO2/1000)/1000; % MJ/kg
    % heat
        Q_spec = QTotal/(N_CO2_product*dataModel.general.MMCO2/1000)/1000; % MJ/kg
    % productivity
        Pr = N_CO2_product*dataModel.general.MMCO2/1000/(dataModel.process.adsorbentMass/dataModel.sorbent.Density)/((dataModel.process.coolingTime+t_ads+dataModel.process.noStepsBD+dataModel.process.noStepsHeating)/3600); % kg_CO2/m3_sorbent/h
    % purity
        purity = N_CO2_product/(outputBlowEvac(end-1,4)); % Nout_save_CO2_evac/(Nout_save_total_evac-Nout_save_H2O_evac)
        purity_dry = N_CO2_product/(outputBlowEvac(end-1,4)-outputBlowEvac(end-1,7)); % Nout_save_CO2_evac/(Nout_save_total_evac-Nout_save_H2O_evac)
    %recovery
        recovery = N_CO2_product/(dataModel.feed.yCO2*sum(outputAdsorption(:,9)) + outputCooling(end,5));
    % cyclic working capacity
        WC = N_CO2_product/dataModel.process.adsorbentMass; % mol_CO2_recovered/kg_sorbent

    sprintf(' Specific Q: %0.2f MJ/kg,\n specific E: %0.2f MJ/kg,\n productivity: %0.2f kg/m3/h_{ads},\n purity: %0.2f,\n purity dry: %0.2f,\n capture rate: %0.2f, \n specific working capacity: %0.2f molCO2/kgSorbent.', Q_spec, E_spec, Pr, purity, purity_dry, recovery, WC)
    
    flag = 1;
    
    out = [flag, E_spec, Q_spec, Pr, purity_dry, recovery];

catch
    flag = 0;
    sprintf('Running 0D model did not work.')
    out = [flag, 0, 0, 0, 0, 0];
end
end

