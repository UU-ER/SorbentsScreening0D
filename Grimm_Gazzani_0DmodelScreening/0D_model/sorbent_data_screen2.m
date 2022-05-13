
%% sorbent data
i = 1;

switch isothermModel
    case 'toth_cp'
            dataModel.sorbent(i).name = material_name;
        % CO2: as an example exemplary isotherm
            dataModel.sorbent(i).CO2Isotherm.model = isothermModel;
            dataModel.sorbent(i).CO2Isotherm.dH_ads_CO2 = 0; % Heat of adsorption [J/mol]
            dataModel.sorbent(i).CO2Isotherm.T0 = T0; % K
        %     param: Xi_c,dH_c[J/mol],alpha_c,Xi_p,dH_p[J/mol],alpha_p,ns0_c,ns0_c[mol/kg],b0_c[1/MPa],t0_c,ns0_p[mol/kg],b0_p[1/MPa],t0_p
            dataModel.sorbent(i).CO2Isotherm.param(1) = p(1); % Xi_c
            dataModel.sorbent(i).CO2Isotherm.param(2) = p(2)*1000; % dH_c[J/mol]
            dataModel.sorbent(i).CO2Isotherm.param(3) = p(3); % alpha_c
            dataModel.sorbent(i).CO2Isotherm.param(4) = p(4); % Xi_p
            dataModel.sorbent(i).CO2Isotherm.param(5) = p(5)*1000; % dH_p[J/mol]
            dataModel.sorbent(i).CO2Isotherm.param(6) = p(6); % alpha_p
            dataModel.sorbent(i).CO2Isotherm.param(7) = p(7); % ns0_c[mol/kg]
            dataModel.sorbent(i).CO2Isotherm.param(8) = p(8); % b0_c[1/MPa]
            dataModel.sorbent(i).CO2Isotherm.param(9) = p(9); % t0_c
            dataModel.sorbent(i).CO2Isotherm.param(10) = p(10); % ns0_p[mol/kg]
            dataModel.sorbent(i).CO2Isotherm.param(11) = p(11); % b0_p[1/MPa]
            dataModel.sorbent(i).CO2Isotherm.param(12) = p(12); % t0_p

        % H2O: as an example CW isotherm (own fitting)
            dataModel.sorbent(i).H2OIsotherm.model = 'GAB';
            dataModel.sorbent(i).H2OIsotherm.dH_ads_water = -49000; % Heat of adsorption [J/mol]
            
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

        % sorbent specific data    --> assumption
        if ~isempty(dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).rho)
            dataModel.sorbent(i).Density = dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).rho; % density of the bed, kg/m3
        else
            dataModel.sorbent(i).Density = 375.4; % density of the bed, kg/m3
        end
        
        if ~isempty(dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).rhoMat)
            dataModel.sorbent(i).MaterialDensity = dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).rhoMat; % density of the bed, kg/m3
        else
            dataModel.sorbent(i).MaterialDensity = 1590; % density material sorbent, kg/m3 
        end
        
        if ~isempty(dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).cp)
            dataModel.sorbent(i).cp = dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).cp; % density of the bed, kg/m3
        else
            dataModel.sorbent(i).cp = 1514.2; % [J/K-kg],specific heat capacity solid (Park),,,
        end
%             dataModel.sorbent(i).MaterialDensity = 1590; % density material sorbent, kg/m3    
%             dataModel.sorbent(i).cp = 1514.2; % [J/K-kg],specific heat capacity solid (Park),,,

    case 's_shaped'
            dataModel.sorbent(i).name = material_name;
        % CO2: as an example exemplary isotherm
        % q_L0(mol/kg), b_L0(1/Pa), dU_L(kJ/mol), q_U0(kJ/mol),b_U0(1/Pa), dU_U(kJ/mol), b_H0(mol/kg/Pa), dU_H(kJ/mol), Xi_1(-), Xi_2(-), pstep_0(Pa), dH_step(kJ/mol), gamma(-)
            dataModel.sorbent(i).CO2Isotherm.model = isothermModel;
            dataModel.sorbent(i).CO2Isotherm.dH_ads_CO2 = 0; % Heat of adsorption [J/mol]\
            dataModel.sorbent(i).CO2Isotherm.T0 = T0; % K
            dataModel.sorbent(i).CO2Isotherm.param(1) = p(1); % q_L0, mol/kg, a1
            dataModel.sorbent(i).CO2Isotherm.param(2) = p(2); % b_L0, 1/MPa, b0 
            dataModel.sorbent(i).CO2Isotherm.param(3) = p(2)*1000; % dU_L, J/mol, b1 
            dataModel.sorbent(i).CO2Isotherm.param(4) = p(4); % q_U0, J/mol, c1 
            dataModel.sorbent(i).CO2Isotherm.param(5) = p(5); % b_U0, 1/MPa d0
            dataModel.sorbent(i).CO2Isotherm.param(6) = p(6)*1000; % dU_U, J/mol d1
            dataModel.sorbent(i).CO2Isotherm.param(7) = p(7); % b_H0, mol/kg/Pa  bb0
            dataModel.sorbent(i).CO2Isotherm.param(8) = p(8)*1000; % dU_H, J/mol bb1
            dataModel.sorbent(i).CO2Isotherm.param(9) = p(9); % xi1
            dataModel.sorbent(i).CO2Isotherm.param(10) = p(10); % xi_2, 1/K xi2
            dataModel.sorbent(i).CO2Isotherm.param(11) = p(11); % p_step0, MPa ps0
            dataModel.sorbent(i).CO2Isotherm.param(12) = p(12)*1000; % dH_step, J/mol Hst
            dataModel.sorbent(i).CO2Isotherm.param(13) = p(13); % gam

        % H2O: as an example CW isotherm (own fitting)
            dataModel.sorbent(i).H2OIsotherm.model = 'GAB';
            dataModel.sorbent(i).H2OIsotherm.dH_ads_water = -49000; % Heat of adsorption [J/mol]

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
            
        if ~isempty(dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).rho)
            dataModel.sorbent(i).Density = dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).rho; % density of the bed, kg/m3
        else
            dataModel.sorbent(i).Density = 375.4; % density of the bed, kg/m3
        end
        
        if ~isempty(dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).rhoMat)
            dataModel.sorbent(i).MaterialDensity = dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).rhoMat; % density of the bed, kg/m3
        else
            dataModel.sorbent(i).MaterialDensity = 1590; % density material sorbent, kg/m3 
        end
        
        if ~isempty(dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).cp)
            dataModel.sorbent(i).cp = dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).cp; % density of the bed, kg/m3
        else
            dataModel.sorbent(i).cp = 1514.2; % [J/K-kg],specific heat capacity solid (Park),,,
        end
        
    case 'DSL'
        % n1 (mol/kg), b0 (m3/mol), Hb (J/mol), n2 (mol/kg), d0 (m3/mol), Hd (J/mol)
             dataModel.sorbent(i).name = material_name; % e.g. dataModel from Farooq
        % CO2: as an example exemplary isotherm
            dataModel.sorbent(i).CO2Isotherm.model = isothermModel;
            dataModel.sorbent(i).CO2Isotherm.dH_ads_CO2 = 0; % Heat of adsorption [J/mol]\
            dataModel.sorbent(i).CO2Isotherm.param(1) = p(1); % n1 (mol/kg)
            dataModel.sorbent(i).CO2Isotherm.param(2) = p(2); % b0 (m3/Mmol)
            dataModel.sorbent(i).CO2Isotherm.param(3) = p(3); % Hb (J/mol) 
            dataModel.sorbent(i).CO2Isotherm.param(4) = p(4); % n2 (mol/kg) 
            dataModel.sorbent(i).CO2Isotherm.param(5) = p(5); % d0 (m3/Mmol)
            dataModel.sorbent(i).CO2Isotherm.param(6) = p(6); % Hd (J/mol)

        % H2O: as an example CW isotherm (own fitting)
            dataModel.sorbent(i).H2OIsotherm.model = 'GAB';
            dataModel.sorbent(i).H2OIsotherm.dH_ads_water = -49000; % Heat of adsorption [J/mol]
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

        % sorbent specific dataModel    --> change!!!!!!!!!!!!!!
            dataModel.sorbent(i).MaterialDensity = dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).rhoMat;
            dataModel.sorbent(i).Density = dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).rhoMat*0.35; % density sorbent, kg/m3 (particle)
            dataModel.sorbent(i).cp = 1514.2; % could also be added from Farooq data                    
            
    case 'langm_freund'
            dataModel.sorbent(i).name = material_name; % e.g. dataModel from Farooq
        % CO2: as an example exemplary isotherm
            dataModel.sorbent(i).CO2Isotherm.T0 = T0; % K
            dataModel.sorbent(i).CO2Isotherm.model = isothermModel;
            dataModel.sorbent(i).CO2Isotherm.dH_ads_CO2 = 0; % Heat of adsorption [J/mol]\
            dataModel.sorbent(i).CO2Isotherm.param(1) = p(1); % ns0 (mol/kg)
            dataModel.sorbent(i).CO2Isotherm.param(2) = p(2); % Ci
            dataModel.sorbent(i).CO2Isotherm.param(3) = p(3); % t0
            dataModel.sorbent(i).CO2Isotherm.param(4) = p(4); % alpha
            dataModel.sorbent(i).CO2Isotherm.param(5) = p(5); % b0 (1/MPa)
            dataModel.sorbent(i).CO2Isotherm.param(6) = p(6)*1000; % dH (J/mol)

        % H2O: as an example CW isotherm (own fitting)
            dataModel.sorbent(i).H2OIsotherm.model = 'GAB';
            dataModel.sorbent(i).H2OIsotherm.dH_ads_water = -49000; % Heat of adsorption [J/mol]
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
            
        if ~isempty(dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).rho)
            dataModel.sorbent(i).Density = dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).rho; % density of the bed, kg/m3
        else
            dataModel.sorbent(i).Density = 375.4; % density of the bed, kg/m3
        end
        
        if ~isempty(dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).rhoMat)
            dataModel.sorbent(i).MaterialDensity = dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).rhoMat; % density of the bed, kg/m3
        else
            dataModel.sorbent(i).MaterialDensity = 1590; % density material sorbent, kg/m3 
        end
        
        if ~isempty(dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).cp)
            dataModel.sorbent(i).cp = dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).cp; % density of the bed, kg/m3
        else
            dataModel.sorbent(i).cp = 1514.2; % [J/K-kg],specific heat capacity solid (Park),,,
        end
            
            
            
    case 'toth'
        dataModel.sorbent(i).name = material_name; % e.g. dataModel from Farooq
        % CO2: as an example exemplary isotherm
            dataModel.sorbent(i).CO2Isotherm.model = 'toth';
            dataModel.sorbent(i).CO2Isotherm.dH_ads_CO2 = 0; % Heat of adsorption [J/mol]
            dataModel.sorbent(i).CO2Isotherm.T0 = T0; % K
            dataModel.sorbent(i).CO2Isotherm.param(1) = p(1); % Xi
            dataModel.sorbent(i).CO2Isotherm.param(2) = p(2)*1000; % dH[J/mol]
            dataModel.sorbent(i).CO2Isotherm.param(3) = p(3); % alpha
            dataModel.sorbent(i).CO2Isotherm.param(4) = p(4); % ns0[mol/kg]
            dataModel.sorbent(i).CO2Isotherm.param(5) = p(5); % b0[1/MPa]
            dataModel.sorbent(i).CO2Isotherm.param(6) = p(6); % t0   
            
        % H2O: as an example CW isotherm (own fitting)
            dataModel.sorbent(i).H2OIsotherm.model = 'GAB';
            dataModel.sorbent(i).H2OIsotherm.dH_ads_water = -49000; % Heat of adsorption [J/mol]
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
            
        if ~isempty(dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).rho)
            dataModel.sorbent(i).Density = dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).rho; % density of the bed, kg/m3
        else
            dataModel.sorbent(i).Density = 375.4; % density of the bed, kg/m3
        end
        
        if ~isempty(dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).rhoMat)
            dataModel.sorbent(i).MaterialDensity = dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).rhoMat; % density of the bed, kg/m3
        else
            dataModel.sorbent(i).MaterialDensity = 1590; % density material sorbent, kg/m3 
        end
        
        if ~isempty(dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).cp)
            dataModel.sorbent(i).cp = dataset.(name_analysis).Carbon_Dioxide.(material_name).(doi_name).cp; % density of the bed, kg/m3
        else
            dataModel.sorbent(i).cp = 1514.2; % [J/K-kg],specific heat capacity solid (Park),,,
        end

end

