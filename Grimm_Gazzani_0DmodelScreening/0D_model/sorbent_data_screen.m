
% =====================================
% exemplary sorbent data as input
i = 1;

switch isothermModel
    case 'toth_cp'
            dataModel.sorbent(i).name = 'exemplary';
        % CO2: as an example exemplary isotherm
            dataModel.sorbent(i).CO2Isotherm.model = isothermModel;
            dataModel.sorbent(i).CO2Isotherm.T0 = 293; % K
        %     param: Xi_c,dH_c[J/mol],alpha_c,Xi_p,dH_p[J/mol],alpha_p,ns0_c,ns0_c[mol/kg],b0_c[1/MPa],t0_c,ns0_p[mol/kg],b0_p[1/MPa],t0_p
            dataModel.sorbent(i).CO2Isotherm.param(1) = 2.17189502654334; % Xi_c
            dataModel.sorbent(i).CO2Isotherm.param(2) = 64086.2167837738; % dH_c[J/mol]
            dataModel.sorbent(i).CO2Isotherm.param(3) = 0.178250162092376; % alpha_c
            dataModel.sorbent(i).CO2Isotherm.param(4) = 0.107279283494825; % Xi_p
            dataModel.sorbent(i).CO2Isotherm.param(5) = 15489.6860407618; % dH_p[J/mol]
            dataModel.sorbent(i).CO2Isotherm.param(6) = 0; % alpha_p
            dataModel.sorbent(i).CO2Isotherm.param(7) = 2.47567728276134; % ns0_c[mol/kg]
            dataModel.sorbent(i).CO2Isotherm.param(8) = 2.81170386473106d6; % b0_c[1/MPa]
            dataModel.sorbent(i).CO2Isotherm.param(9) = 0.275138295967063; % t0_c
            dataModel.sorbent(i).CO2Isotherm.param(10) = 5.59698351121952; % ns0_p[mol/kg]
            dataModel.sorbent(i).CO2Isotherm.param(11) = 3.32264994071570; % b0_p[1/MPa]
            dataModel.sorbent(i).CO2Isotherm.param(12) = 0.308303702506619; % t0_p

        % H2O: as an example CW isotherm (own fitting)
            dataModel.sorbent(i).H2OIsotherm.model = 'GAB';
            dataModel.sorbent(i).H2OIsotherm.dH_ads_water = -49000; % Heat of adsorption [J/mol]
            
            own_H2O = false;
            if own_H2O 
                dataModel.sorbent(i).H2OIsotherm.param(1) = pH2O(1); % CG0
                dataModel.sorbent(i).H2OIsotherm.param(2) = pH2O(2); % HC[J/mol]
                dataModel.sorbent(i).H2OIsotherm.param(3) = pH2O(3); % K0
                dataModel.sorbent(i).H2OIsotherm.param(4) = pH2O(4); % HK[J/mol]
                dataModel.sorbent(i).H2OIsotherm.param(5) = pH2O(5); % Cm0[mol/kg]
                dataModel.sorbent(i).H2OIsotherm.param(6) = pH2O(6); % beta[K]
            else
                dataModel.sorbent(i).H2OIsotherm.param(1) = 6.86; % CG0
                dataModel.sorbent(i).H2OIsotherm.param(2) = -5088; % HC[J/mol]
                dataModel.sorbent(i).H2OIsotherm.param(3) = 2.27; % K0
                dataModel.sorbent(i).H2OIsotherm.param(4) = -3443; % HK[J/mol]
                dataModel.sorbent(i).H2OIsotherm.param(5) = 0.0208; % Cm0[mol/kg]
                dataModel.sorbent(i).H2OIsotherm.param(6) = 1797; % beta[K]
            end
        dataModel.sorbent(i).Density = 375.4; % density of the bed, kg/m3
        dataModel.sorbent(i).MaterialDensity = 1590; % density material sorbent, kg/m3 
        dataModel.sorbent(i).cp = 1514.2; % [J/K-kg],specific heat capacity solid (Park),,,
    case 's_shaped'
    case 'toth'

end

