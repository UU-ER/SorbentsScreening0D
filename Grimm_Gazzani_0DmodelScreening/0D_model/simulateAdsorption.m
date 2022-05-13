
%% simulate adsorption step

% two steps:
% 1. saturate with H2O
% 2. saturate with CO2 (start with 100% saturation)

function outAdsorption = simulateAdsorption(dataModel,outputCool)
    
    %% GET ISOTHERM PARAMETERS
    CO2IsothermModel = dataModel.currentSorbent;
    i = 1;

    switch CO2IsothermModel
        case 'toth_cp'
			% CO2 isotherm
			T0 = dataModel.sorbent(i).CO2Isotherm.T0(1);
			Xi_c = dataModel.sorbent(i).CO2Isotherm.param(1);
			dH_c = dataModel.sorbent(i).CO2Isotherm.param(2); % [J/mol]
			alpha_c = dataModel.sorbent(i).CO2Isotherm.param(3);
			Xi_p = dataModel.sorbent(i).CO2Isotherm.param(4);
			dH_p = dataModel.sorbent(i).CO2Isotherm.param(5); % [J/mol]
			alpha_p = dataModel.sorbent(i).CO2Isotherm.param(6);
			ns0_c = dataModel.sorbent(i).CO2Isotherm.param(7); % [mol/kg]
			b0_c = dataModel.sorbent(i).CO2Isotherm.param(8)*10^(6); % [1/MPa]
			t0_c = dataModel.sorbent(i).CO2Isotherm.param(9);
			ns0_p = dataModel.sorbent(i).CO2Isotherm.param(10); % [mol/kg]
			b0_p = dataModel.sorbent(i).CO2Isotherm.param(11)*10^(6); % [1/MPa]
			t0_p = dataModel.sorbent(i).CO2Isotherm.param(12);       
        case 's_shaped'
			T0 = dataModel.sorbent(i).CO2Isotherm.T0; % K
			q_L0 = dataModel.sorbent(i).CO2Isotherm.param(1); % q_L0, mol/kg, a1
			b_L0 = dataModel.sorbent(i).CO2Isotherm.param(2)*10^(6); % b_L0, 1/MPa, b0 
			dU_L = dataModel.sorbent(i).CO2Isotherm.param(3); % dU_L, J/mol, b1 
			q_U0 = dataModel.sorbent(i).CO2Isotherm.param(4); % q_U0, J/mol, c1 
			b_U0 = dataModel.sorbent(i).CO2Isotherm.param(5)*10^(6); % b_U0, 1/MPa d0
			dU_U = dataModel.sorbent(i).CO2Isotherm.param(6); % dU_U, J/mol d1
			b_H0 = dataModel.sorbent(i).CO2Isotherm.param(7); % b_H0, mol/kg/Pa  bb0
			dU_H = dataModel.sorbent(i).CO2Isotherm.param(8); % dU_H, J/mol bb1
			xi_1 = dataModel.sorbent(i).CO2Isotherm.param(9); % xi1
			xi_2 = dataModel.sorbent(i).CO2Isotherm.param(10); % xi_2, 1/K xi2
			p_step0 = dataModel.sorbent(i).CO2Isotherm.param(11)*10^(-6); % p_step0, MPa ps0
			dH_step = dataModel.sorbent(i).CO2Isotherm.param(12); % dH_step, J/mol Hst
			gam = dataModel.sorbent(i).CO2Isotherm.param(13); % gam            
        case 'DSL' % the sorbents by Farooq are already in J/mol -> no need for *1000
			n1 = dataModel.sorbent(i).CO2Isotherm.param(1); % n1 (mol/kg)
			b0 = dataModel.sorbent(i).CO2Isotherm.param(2)/(1e6); % b0 (m3/mol)
			Hb = dataModel.sorbent(i).CO2Isotherm.param(3); % Hb (J/mol) 
			n2 = dataModel.sorbent(i).CO2Isotherm.param(4); % n2 (mol/kg)
			d0 = dataModel.sorbent(i).CO2Isotherm.param(5)/(1e6); % d0 (m3/mol)
			Hd = dataModel.sorbent(i).CO2Isotherm.param(6); % Hd (J/mol)
        case 'DSL2'
			T0 = dataModel.sorbent(i).CO2Isotherm.T0; % K
			n1 = dataModel.sorbent(i).CO2Isotherm.param(1); % n1 (mol/kg)
			b0 = dataModel.sorbent(i).CO2Isotherm.param(2); % b0 (m3/mol)
			Hb = dataModel.sorbent(i).CO2Isotherm.param(3); % Hb (J/mol) 
			n2 = dataModel.sorbent(i).CO2Isotherm.param(4); % n2 (mol/kg) 
			d0 = dataModel.sorbent(i).CO2Isotherm.param(5); % d0 (m3/mol)
			Hd = dataModel.sorbent(i).CO2Isotherm.param(6); % Hd (J/mol)
        case 'toth'
			% CO2 isotherm
			T0 = dataModel.sorbent(i).CO2Isotherm.T0(1);
			Xi = dataModel.sorbent(i).CO2Isotherm.param(1);
			dH = dataModel.sorbent(i).CO2Isotherm.param(2); % [J/mol]
			alpha = dataModel.sorbent(i).CO2Isotherm.param(3);
			ns0 = dataModel.sorbent(i).CO2Isotherm.param(4); % [mol/kg]
			b0 = dataModel.sorbent(i).CO2Isotherm.param(5)*10^(6); % [1/MPa]
			t0 = dataModel.sorbent(i).CO2Isotherm.param(6);
		case 'langfr'
			% CO2 isotherm
			T0 = dataModel.sorbent(i).CO2Isotherm.T0(1);
			ns0 = dataModel.sorbent(i).CO2Isotherm.param(1);
			Xi = dataModel.sorbent(i).CO2Isotherm.param(2);
			t0 = dataModel.sorbent(i).CO2Isotherm.param(3);
			alpha = dataModel.sorbent(i).CO2Isotherm.param(4); % [mol/kg]
			b0 = dataModel.sorbent(i).CO2Isotherm.param(5)*10^(6); % [1/MPa]
			dH = dataModel.sorbent(i).CO2Isotherm.param(6); % [J/mol]
    end       
	
    if dataModel.feed.yH2O > 1e-10
        % H2O isotherm
        CG0 = dataModel.sorbent(i).H2OIsotherm.param(1);
        HC = dataModel.sorbent(i).H2OIsotherm.param(2); % J/mol
        K0 = dataModel.sorbent(i).H2OIsotherm.param(3);
        HK = dataModel.sorbent(i).H2OIsotherm.param(4); % J/mol
        Cm0 = dataModel.sorbent(i).H2OIsotherm.param(5); % mol/kg
        beta = dataModel.sorbent(i).H2OIsotherm.param(6); % K    
    end
    
    %% solver options
    options = optimoptions('lsqnonlin','Display','off');
    options.Algorithm = 'trust-region-reflective'; % levenberg-marquardt, trust-region-reflective
    options.OptimalityTolerance = 1e-16;
    options.FunctionTolerance = 1e-16;
    options.StepTolerance = 1e-16;
    options.MaxFunctionEvaluations = 6000;    
    
    %% RUN ADSORPTION MODEL
   
	Tamb = dataModel.process.Tamb;   
    p_amb = dataModel.process.pamb;
	T = dataModel.process.Tamb;
    yCO2_feed = dataModel.feed.yCO2;
    yN2_feed = dataModel.feed.yN2;
    yH2O_feed = dataModel.feed.yH2O;
    
    yCO2_k = outputCool(end,1);
    yN2_k = outputCool(end,2);
    yH2O_k = outputCool(end,3);
    
%% First adsorption step (saturate with H2O) ------------------------------
    if dataModel.feed.yH2O > 1e-10
        % Initial condition 
        x0 = [1000, yCO2_k, yN2_k, 1000];

        % boundaries
        lbDim = [1, 0, 0, 1];
        ubDim = [5500, 1, 1, 5500];
		
		x0 = (x0-lbDim)./(ubDim-lbDim); % normalized 

        lb = [0, 0, 0, 0];
        ub = [1, 1, 1, 1];
    else

    end   

    %% Run and check for strongly adsorbed component
    % 1.1 Assumption H2O is strongly adsorbed component --------------------
    if dataModel.feed.yH2O > 1e-10
        % calculate composition at given saturation
        yH2O_sat = yH2O_feed;

        funct = @(x) double( adsorption1(x,dataModel.process.adsorbentMass,dataModel.process.voidFraction,dataModel.process.Vol,dataModel.general.gasconstant,yH2O_sat,yCO2_k,yN2_k,yH2O_k,lbDim,ubDim));
        fun = @(x) funct(x);

        [x,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,lb,ub,options);

        % dimensioned
        xopt = x;
        x = xopt.*(ubDim-lbDim)+lbDim;

        N_feed_1 = x(1);
        yCO2_save_1 = x(2);
        yN2_save_1 = x(3); % 
        yH2O_save_1 = yH2O_sat; % 
        Nwaste_save_1 = x(4);
        x_save(:) = residual;

        RH_save = dataModel.process.pamb*yH2O_save_1/vaporPressure(dataModel.process.Tamb)*100;

        if RH_save > 100
            sprintf('Humidity can not be higher than 100%%.')
        end

        % Energy consumption till that given step
        E_total(1) = cmp_air(dataModel.process.pamb,dataModel.general.densityAir,x(1),dataModel.general.MMCO2); % J, air blower
	        
        % check if already saturated with CO2
        yCO2_sat = dataModel.startingCondBD(1);
        if x(2)> yCO2_sat % if yes, then run again fixing yCO2 sat
            % 1.2 Assumption H2O is strongly adsorbed component --------------------
            funct = @(x) double( adsorption11(x,dataModel.process.adsorbentMass,dataModel.process.voidFraction,dataModel.process.Vol,dataModel.general.gasconstant,yCO2_sat,yCO2_k,yN2_k,yH2O_k,lbDim,ubDim));
            fun = @(x) funct(x);

            [x,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,lb,ub,options);
            
            % dimensioned
            xopt = x;
            x = xopt.*(ubDim-lbDim)+lbDim;

            N_feed_1 = x(1);
            yCO2_save_1 = yCO2_sat;
            yN2_save_1 = x(3); % 
            yH2O_save_1 = x(2); % 
            Nwaste_save_1 = x(4);
            x_save(:) = residual;
    %         x_save(:) = fvec;

            RH_save = dataModel.process.pamb*yH2O_save_1/vaporPressure(dataModel.process.Tamb)*100;

            if RH_save > 100
                sprintf('Humidity can not be higher than 100%%.')
            end

            % Energy consumption till that given step
            E_total(1) = cmp_air(dataModel.process.pamb,dataModel.general.densityAir,x(1),dataModel.general.MMCO2); % J, air blower
            
            yCO2_save = [yCO2_save_1]; 
            yN2_save = [yN2_save_1];
            yH2O_save = [yH2O_save_1];

            Nwaste_save_sum(1) = Nwaste_save_1;
            Nfeed_save_sum(1) = N_feed_1;
            Nwaste_save_CO2(1) = Nwaste_save_1*yCO2_save_1;
            Nwaste_save_N2(1) = Nwaste_save_1*yN2_save_1;
            Nwaste_save_H2O(1) = Nwaste_save_1*yH2O_save_1; 
            
        else
            Nwaste_save_sum(1) = Nwaste_save_1;
            Nfeed_save_sum(1) = N_feed_1;
            Nwaste_save_CO2(1) = Nwaste_save_1*yCO2_save_1;
            Nwaste_save_N2(1) = Nwaste_save_1*yN2_save_1;
            Nwaste_save_H2O(1) = Nwaste_save_1*yH2O_save_1; 
        %% Second adsorption step (saturate with CO2) -----------------------------
            % Initial condition 
            % Nin, yN2, Nout
            x0 = [N_feed_1, yN2_k, Nwaste_save_1];
            % boundaries
            lbDim = [0, 0, 0];
            ubDim = [5500, 1, 5500];

            x0 = (x0-lbDim)./(ubDim-lbDim); % normalized 

            lb = [0, 0, 0];
            ub = [1, 1, 1];

            yCO2_k = yCO2_save_1;
            yN2_k = yN2_save_1;
            yH2O_k = yH2O_save_1;

            %% solve equations
            % CO2 concentration at saturation level
            yCO2_sat = dataModel.startingCondBD(1);

            if dataModel.feed.yH2O > 1e-10
                funct = @(x) double( adsorption2(x,dataModel.process.adsorbentMass,dataModel.process.voidFraction,dataModel.process.Vol,dataModel.general.gasconstant,yCO2_sat,yCO2_k,yN2_k,yH2O_k,lbDim,ubDim));
                fun = @(x) funct(x);
            else
                funct = @(x) double( adsorption2_dry(x,dataModel.process.adsorbentMass,dataModel.process.voidFraction,dataModel.process.Vol,dataModel.general.gasconstant,yCO2_sat,yCO2_k,yN2_k,lbDim,ubDim));
                fun = @(x) funct(x);
            end

            [x,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,lb,ub,options);

            % dimensioned
            xopt = x;
            x = xopt.*(ubDim-lbDim)+lbDim;

            N_feed_2 = x(1);
            yCO2_save_2 = yCO2_sat;
            yN2_save_2 = x(2); % 
            yH2O_save_2 = yH2O_sat; % 
            Nwaste_save_2 = x(3);
            x_save2(:) = residual;
   
            Nwaste_save_sum(2) = Nwaste_save_2;
            Nfeed_save_sum(2) = N_feed_2;
            Nwaste_save_CO2(2) = Nwaste_save_2*yCO2_save_2;
            Nwaste_save_N2(2) = Nwaste_save_2*yN2_save_2;
            Nwaste_save_H2O(2) = Nwaste_save_2*yH2O_save_2;  

            RH_save = dataModel.process.pamb*yH2O_save_1/vaporPressure(dataModel.process.Tamb)*100;

            if RH_save > 100
                sprintf('Humidity can not be higher than 100%%.')
            end

            % Energy consumption till that given step
            E_total(2) = cmp_air(dataModel.process.pamb,dataModel.general.densityAir,x(1),dataModel.general.MMCO2); % J, air blower      

            %% Combine both steps
            yCO2_save = [yCO2_save_1; yCO2_save_2]; 
            yN2_save = [yN2_save_1; yN2_save_2];
            yH2O_save = [yH2O_save_1; yH2O_save_2];
        end
    else
        yH2O_sat = 0;
        N_feed_1 = 0;
        yCO2_save_1 = yCO2_k;
        yN2_save_1 = yN2_k; % 
        yH2O_save_1 = 0; % 
        Nwaste_save_1 = 0;

        E_total(1) = 0;
    end 
    
outAdsorption = [yCO2_save  yN2_save  yH2O_save  Nwaste_save_sum'  Nwaste_save_CO2'  Nwaste_save_N2'  Nwaste_save_H2O'  E_total'  Nfeed_save_sum'];    

%% functions
function fvec = adsorption1(x,adsorbentMass,void,V,R,yH2O_sat,yCO2_cool,yN2_cool,yH2O_cool,lb,ub)
    % x(1): N_feed, x(2): yCO2, x(3): yN2, x(4): N_waste

    %% equation 1: material balance CO2
    % CO2 solid and fluid from cooling
        % solid phase CO2, final condition 
        switch CO2IsothermModel
            case 'toth_cp'
                qCO2 = ((ns0_c*exp(Xi_c*(1-T/T0)).*(b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p_amb*yCO2_cool))./((1+((b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p_amb*yCO2_cool)).^(t0_c+alpha_c*(1-T0./T))).^(1./(t0_c+alpha_c*(1-T0./T)))))+((ns0_p*exp(Xi_p*(1-Tamb/T0)).*(b0_p*exp(dH_p/(R*T0)*(T0./Tamb-1))).*(p_amb*yCO2_cool))./((1+((b0_p*exp(dH_p/(R*T0)*(T0./Tamb-1))).*(p_amb*yCO2_cool)).^(t0_p+alpha_p*(1-T0./Tamb))).^(1./(t0_p+alpha_p*(1-T0./Tamb))))); % adsorbed amount of CO2
            case 's_shaped'            
                qCO2 = (q_L0.*(b_L0.*exp(dU_L./(R*T))).*(yCO2_cool*p_amb)./(1+(b_L0.*exp(dU_L./(R*T))).*(yCO2_cool*p_amb))).*(1-((exp((log((yCO2_cool*p_amb))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((yCO2_cool*p_amb)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*T))).*(yCO2_cool*p_amb)./(1+(b_U0.*exp(dU_U./(R*T))).*(yCO2_cool*p_amb))+(b_H0.*exp(dU_H./(R*T))).*(yCO2_cool*p_amb)).*((exp((log((yCO2_cool*p_amb))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((yCO2_cool*p_amb)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam); 
            case 'DSL'
                qCO2 = n1*b0*exp(Hb/R/T)*(yCO2_cool*p_amb)./(1e-6*R*T)./(1+b0*exp(Hb/R/T).*(yCO2_cool*p_amb)./(1e-6*R*T)) + n2*d0*exp(Hd/R/T)*(yCO2_cool*p_amb)./(1e-6*R*T)./(1+d0*exp(Hd/R/T).*(yCO2_cool*p_amb)./(1e-6*R*T));            
            case 'DSL2'
                qCO2 = n1*b0*exp(Hb/R/T)*(yCO2_cool*p_amb)./(1+b0*exp(Hb/R/T).*(yCO2_cool*p_amb)) + n2*d0*exp(Hd/R/T)*(yCO2_cool*p_amb)./(1+d0*exp(Hd/R/T).*(yCO2_cool*p_amb));            
             case 'toth'
                qCO2 = ((ns0*exp(Xi*(1-T/T0)).*(b0*exp(dH/(R*T0)*(T0./T-1))).*(p_amb*yCO2_cool))./((1+((b0*exp(dH/(R*T0)*(T0./T-1))).*(p_amb*yCO2_cool)).^(t0+alpha*(1-T0./T))).^(1./(t0+alpha*(1-T0./T)))));           
            case 'langfr'
                qCO2 = (ns0.*exp(Xi.*(1-T./T0))).*((b0*exp(dH./(R*T0).*(T0./T-1))).*(p_amb*yCO2_cool)).^(1./t0+alpha.*(1-T0./T))./(1+ ((b0*exp(dH./(R*T0).*(T0./T-1))).*(p_amb*yCO2_cool)).^(1./t0+alpha.*(1-T0./T)));
        end        
        NCO2_solid = adsorbentMass*qCO2;
        % fluid phase CO2, final condition  
        NCO2_fluid = yCO2_cool*p_amb*V*void/(R*T)*1e6; % mol
        % solid and fluid CO2, final condition  
        N_CO2_cool = NCO2_solid + NCO2_fluid;

    % solid and fluid at end
        % solid phase CO2, final condition 
        switch CO2IsothermModel
            case 'toth_cp'
                qCO2 = ((ns0_c*exp(Xi_c*(1-T/T0)).*(b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p_amb*(x(2)*(ub(2)-lb(2))+lb(2))))./((1+((b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p_amb*(x(2)*(ub(2)-lb(2))+lb(2)))).^(t0_c+alpha_c*(1-T0./T))).^(1./(t0_c+alpha_c*(1-T0./T)))))+((ns0_p*exp(Xi_p*(1-Tamb/T0)).*(b0_p*exp(dH_p/(R*T0)*(T0./Tamb-1))).*(p_amb*(x(2)*(ub(2)-lb(2))+lb(2))))./((1+((b0_p*exp(dH_p/(R*T0)*(T0./Tamb-1))).*(p_amb*(x(2)*(ub(2)-lb(2))+lb(2)))).^(t0_p+alpha_p*(1-T0./Tamb))).^(1./(t0_p+alpha_p*(1-T0./Tamb))))); % adsorbed amount of CO2
            case 's_shaped'            
                qCO2 = (q_L0.*(b_L0.*exp(dU_L./(R*T))).*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)./(1+(b_L0.*exp(dU_L./(R*T))).*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb))).*(1-((exp((log(((x(2)*(ub(2)-lb(2))+lb(2))*p_amb))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log(((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*T))).*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)./(1+(b_U0.*exp(dU_U./(R*T))).*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb))+(b_H0.*exp(dU_H./(R*T))).*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)).*((exp((log(((x(2)*(ub(2)-lb(2))+lb(2))*p_amb))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log(((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam); 
            case 'DSL'
                qCO2 = n1*b0*exp(Hb/R/T)*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)./(1e-6*R*T)./(1+b0*exp(Hb/R/T).*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)./(1e-6*R*T)) + n2*d0*exp(Hd/R/T)*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)./(1e-6*R*T)./(1+d0*exp(Hd/R/T).*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)./(1e-6*R*T));            
            case 'DSL2'
                qCO2 = n1*b0*exp(Hb/R/T)*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)./(1+b0*exp(Hb/R/T).*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)) + n2*d0*exp(Hd/R/T)*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb)./(1+d0*exp(Hd/R/T).*((x(2)*(ub(2)-lb(2))+lb(2))*p_amb));            
            case 'toth'
                qCO2 = ((ns0*exp(Xi*(1-T/T0)).*(b0*exp(dH/(R*T0)*(T0./T-1))).*(p_amb*(x(2)*(ub(2)-lb(2))+lb(2))))./((1+((b0*exp(dH/(R*T0)*(T0./T-1))).*(p_amb*(x(2)*(ub(2)-lb(2))+lb(2)))).^(t0+alpha*(1-T0./T))).^(1./(t0+alpha*(1-T0./T)))));                       
            case 'langfr'
                qCO2 = (ns0.*exp(Xi.*(1-T./T0))).*((b0*exp(dH./(R*T0).*(T0./T-1))).*(p_amb*(x(2)*(ub(2)-lb(2))+lb(2)))).^(1./t0+alpha.*(1-T0./T))./(1+ ((b0*exp(dH./(R*T0).*(T0./T-1))).*(p_amb*(x(2)*(ub(2)-lb(2))+lb(2)))).^(1./t0+alpha.*(1-T0./T)));
        end        
        NCO2_solid = adsorbentMass*qCO2;
        % fluid phase CO2, final condition  
        NCO2_fluid = (x(2)*(ub(2)-lb(2))+lb(2))*p_amb*V*void/(R*T)*1e6; % mol
        % solid and fluid CO2, final condition  
        N_CO2_end = NCO2_solid + NCO2_fluid;
    
    fvec(1) = N_CO2_cool + yCO2_feed*(x(1)*(ub(1)-lb(1))+lb(1)) - N_CO2_end - yCO2_cool * (x(4)*(ub(4)-lb(4))+lb(4)); % = 0

    %% equation 2: material balance H2O
    % H2O solid and fluid from cooling
            % solid phase H2O, final condition 
            pvap = ((22.064 * exp((647.096/Tamb)*((-7.85951783)*(1-(Tamb/647.096))+(1.84408259)*(1-(Tamb/647.096))^(1.5)+(-11.7866497)*(1-(Tamb/647.096))^3+(22.6807411)*(1-(Tamb/647.096))^(3.5)+(-15.9618719)*(1-(Tamb/647.096))^4+(1.80122502)*(1-(Tamb/647.096))^(7.5))))); % MPa
            qH2O = (Cm0*exp(beta/Tamb))*(CG0*exp(HC/R/Tamb))*(K0*exp(HK/R/Tamb))*((p_amb*yH2O_cool)./pvap)./((1-(K0*exp(HK/R/Tamb))*((p_amb*yH2O_cool)./pvap)).*(1+((CG0*exp(HC/R/Tamb))-1)*(K0*exp(HK/R/Tamb))*((p_amb*yH2O_cool)./pvap)));
            NH2O_solid = adsorbentMass*qH2O;
            % fluid phase H2O, final condition 
            NH2O_fluid = yH2O_cool*p_amb*V*void/(R*Tamb)*1e6; % mol
            % solid and liquid H2O, final condition  
            N_H2O_cool = NH2O_solid + NH2O_fluid;

    % solid and fluid at current point 
            % solid phase H2O, final condition 
            pvap = ((22.064 * exp((647.096/Tamb)*((-7.85951783)*(1-(Tamb/647.096))+(1.84408259)*(1-(Tamb/647.096))^(1.5)+(-11.7866497)*(1-(Tamb/647.096))^3+(22.6807411)*(1-(Tamb/647.096))^(3.5)+(-15.9618719)*(1-(Tamb/647.096))^4+(1.80122502)*(1-(Tamb/647.096))^(7.5))))); % MPa
            qH2O = (Cm0*exp(beta/Tamb))*(CG0*exp(HC/R/Tamb))*(K0*exp(HK/R/Tamb))*((p_amb*yH2O_sat)./pvap)./((1-(K0*exp(HK/R/Tamb))*((p_amb*yH2O_sat)./pvap)).*(1+((CG0*exp(HC/R/Tamb))-1)*(K0*exp(HK/R/Tamb))*((p_amb*yH2O_sat)./pvap)));
            NH2O_solid = adsorbentMass*qH2O;
            % fluid phase H2O, final condition 
            NH2O_fluid = yH2O_sat*p_amb*V*void/(R*Tamb)*1e6; % mol
            % solid and fluid H2O, final condition  
            N_H2O_end = NH2O_solid + NH2O_fluid;

    fvec(2) = N_H2O_cool + yH2O_feed*(x(1)*(ub(1)-lb(1))+lb(1)) - N_H2O_end - yH2O_cool * (x(4)*(ub(4)-lb(4))+lb(4)); % = 0

    %% equation 3: material balance N2
    % N2 solid and fluid from previous step
        % fluid phase N2, final condition 
        NN2_fluid = yN2_cool*p_amb*V*void/(R*Tamb)*1e6; % mol
        N_N2_cool = NN2_fluid;

    % fluid at current step (
        % fluid phase N2, final condition 
        NN2_fluid = (x(3)*(ub(3)-lb(3))+lb(3))*p_amb*V*void/(R*Tamb)*1e6; % mol
        N_N2_end = NN2_fluid;

    fvec(3) = N_N2_cool + yN2_feed*(x(1)*(ub(1)-lb(1))+lb(1)) - N_N2_end - yN2_cool * (x(4)*(ub(4)-lb(4))+lb(4)); % = 0

    %% equation 4: overall concentration balance
   fvec(4) = (x(2)*(ub(2)-lb(2))+lb(2)) + (x(3)*(ub(3)-lb(3))+lb(3)) + yH2O_sat - 1; % = 0

    fvec = real(fvec);
end

function fvec = adsorption11(x,adsorbentMass,void,V,R,yCO2_sat,yCO2_cool,yN2_cool,yH2O_cool,lb,ub)
    % x(1): N_feed, x(2): yCO2, x(3): yN2, x(4): N_waste

    %% equation 1: material balance CO2
    % CO2 solid and fluid from cooling
        % solid phase CO2, final condition 
        switch CO2IsothermModel
            case 'toth_cp'
                qCO2 = ((ns0_c*exp(Xi_c*(1-T/T0)).*(b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p_amb*yCO2_cool))./((1+((b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p_amb*yCO2_cool)).^(t0_c+alpha_c*(1-T0./T))).^(1./(t0_c+alpha_c*(1-T0./T)))))+((ns0_p*exp(Xi_p*(1-Tamb/T0)).*(b0_p*exp(dH_p/(R*T0)*(T0./Tamb-1))).*(p_amb*yCO2_cool))./((1+((b0_p*exp(dH_p/(R*T0)*(T0./Tamb-1))).*(p_amb*yCO2_cool)).^(t0_p+alpha_p*(1-T0./Tamb))).^(1./(t0_p+alpha_p*(1-T0./Tamb))))); % adsorbed amount of CO2
            case 's_shaped'            
                qCO2 = (q_L0.*(b_L0.*exp(dU_L./(R*T))).*(yCO2_cool*p_amb)./(1+(b_L0.*exp(dU_L./(R*T))).*(yCO2_cool*p_amb))).*(1-((exp((log((yCO2_cool*p_amb))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((yCO2_cool*p_amb)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*T))).*(yCO2_cool*p_amb)./(1+(b_U0.*exp(dU_U./(R*T))).*(yCO2_cool*p_amb))+(b_H0.*exp(dU_H./(R*T))).*(yCO2_cool*p_amb)).*((exp((log((yCO2_cool*p_amb))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((yCO2_cool*p_amb)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam); 
            case 'DSL'
                qCO2 = n1*b0*exp(Hb/R/T)*(yCO2_cool*p_amb)./(1e-6*R*T)./(1+b0*exp(Hb/R/T).*(yCO2_cool*p_amb)./(1e-6*R*T)) + n2*d0*exp(Hd/R/T)*(yCO2_cool*p_amb)./(1e-6*R*T)./(1+d0*exp(Hd/R/T).*(yCO2_cool*p_amb)./(1e-6*R*T));            
            case 'DSL2'
                qCO2 = n1*b0*exp(Hb/R/T)*(yCO2_cool*p_amb)./(1+b0*exp(Hb/R/T).*(yCO2_cool*p_amb)) + n2*d0*exp(Hd/R/T)*(yCO2_cool*p_amb)./(1+d0*exp(Hd/R/T).*(yCO2_cool*p_amb));            
             case 'toth'
                qCO2 = ((ns0*exp(Xi*(1-T/T0)).*(b0*exp(dH/(R*T0)*(T0./T-1))).*(p_amb*yCO2_cool))./((1+((b0*exp(dH/(R*T0)*(T0./T-1))).*(p_amb*yCO2_cool)).^(t0+alpha*(1-T0./T))).^(1./(t0+alpha*(1-T0./T)))));           
            case 'langfr'
                qCO2 = (ns0.*exp(Xi.*(1-T./T0))).*((b0*exp(dH./(R*T0).*(T0./T-1))).*(p_amb*yCO2_cool)).^(1./t0+alpha.*(1-T0./T))./(1+ ((b0*exp(dH./(R*T0).*(T0./T-1))).*(p_amb*yCO2_cool)).^(1./t0+alpha.*(1-T0./T)));
        end        
        NCO2_solid = adsorbentMass*qCO2;
        % fluid phase CO2, final condition  
        NCO2_fluid = yCO2_cool*p_amb*V*void/(R*T)*1e6; % mol
        % solid and fluid CO2, final condition  
        N_CO2_cool = NCO2_solid + NCO2_fluid;

    % solid and fluid at end
        % solid phase CO2, final condition 
        switch CO2IsothermModel
            case 'toth_cp'
                qCO2 = ((ns0_c*exp(Xi_c*(1-T/T0)).*(b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p_amb*yCO2_sat))./((1+((b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p_amb*yCO2_sat)).^(t0_c+alpha_c*(1-T0./T))).^(1./(t0_c+alpha_c*(1-T0./T)))))+((ns0_p*exp(Xi_p*(1-Tamb/T0)).*(b0_p*exp(dH_p/(R*T0)*(T0./Tamb-1))).*(p_amb*yCO2_sat))./((1+((b0_p*exp(dH_p/(R*T0)*(T0./Tamb-1))).*(p_amb*yCO2_sat)).^(t0_p+alpha_p*(1-T0./Tamb))).^(1./(t0_p+alpha_p*(1-T0./Tamb))))); % adsorbed amount of CO2
            case 's_shaped'            
                qCO2 = (q_L0.*(b_L0.*exp(dU_L./(R*T))).*(yCO2_sat*p_amb)./(1+(b_L0.*exp(dU_L./(R*T))).*(yCO2_sat*p_amb))).*(1-((exp((log((yCO2_sat*p_amb))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((yCO2_sat*p_amb)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*T))).*(yCO2_sat*p_amb)./(1+(b_U0.*exp(dU_U./(R*T))).*(yCO2_sat*p_amb))+(b_H0.*exp(dU_H./(R*T))).*(yCO2_sat*p_amb)).*((exp((log((yCO2_sat*p_amb))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((yCO2_sat*p_amb)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam); 
            case 'DSL'
                qCO2 = n1*b0*exp(Hb/R/T)*(yCO2_sat*p_amb)./(1e-6*R*T)./(1+b0*exp(Hb/R/T).*(yCO2_sat*p_amb)./(1e-6*R*T)) + n2*d0*exp(Hd/R/T)*(yCO2_sat*p_amb)./(1e-6*R*T)./(1+d0*exp(Hd/R/T).*(yCO2_sat*p_amb)./(1e-6*R*T));            
            case 'DSL2'
                qCO2 = n1*b0*exp(Hb/R/T)*(yCO2_sat*p_amb)./(1+b0*exp(Hb/R/T).*(yCO2_sat*p_amb)) + n2*d0*exp(Hd/R/T)*(yCO2_sat*p_amb)./(1+d0*exp(Hd/R/T).*(yCO2_sat*p_amb));            
            case 'toth'
                qCO2 = ((ns0*exp(Xi*(1-T/T0)).*(b0*exp(dH/(R*T0)*(T0./T-1))).*(p_amb*yCO2_sat))./((1+((b0*exp(dH/(R*T0)*(T0./T-1))).*(p_amb*yCO2_sat)).^(t0+alpha*(1-T0./T))).^(1./(t0+alpha*(1-T0./T)))));                       
            case 'langfr'
                qCO2 = (ns0.*exp(Xi.*(1-T./T0))).*((b0*exp(dH./(R*T0).*(T0./T-1))).*(p_amb*yCO2_sat)).^(1./t0+alpha.*(1-T0./T))./(1+ ((b0*exp(dH./(R*T0).*(T0./T-1))).*(p_amb*yCO2_sat)).^(1./t0+alpha.*(1-T0./T)));
        end        
        NCO2_solid = adsorbentMass*qCO2;
        % fluid phase CO2, final condition  
        NCO2_fluid = yCO2_sat*p_amb*V*void/(R*T)*1e6; % mol
        % solid and fluid CO2, final condition  
        N_CO2_end = NCO2_solid + NCO2_fluid;
    
    fvec(1) = N_CO2_cool + yCO2_feed*(x(1)*(ub(1)-lb(1))+lb(1)) - N_CO2_end - yCO2_cool * (x(4)*(ub(4)-lb(4))+lb(4)); % = 0

    %% equation 2: material balance H2O
    % H2O solid and fluid from cooling
            pvap = ((22.064 * exp((647.096/Tamb)*((-7.85951783)*(1-(Tamb/647.096))+(1.84408259)*(1-(Tamb/647.096))^(1.5)+(-11.7866497)*(1-(Tamb/647.096))^3+(22.6807411)*(1-(Tamb/647.096))^(3.5)+(-15.9618719)*(1-(Tamb/647.096))^4+(1.80122502)*(1-(Tamb/647.096))^(7.5))))); % MPa
            qH2O = (Cm0*exp(beta/Tamb))*(CG0*exp(HC/R/Tamb))*(K0*exp(HK/R/Tamb))*((p_amb*yH2O_cool)./pvap)./((1-(K0*exp(HK/R/Tamb))*((p_amb*yH2O_cool)./pvap)).*(1+((CG0*exp(HC/R/Tamb))-1)*(K0*exp(HK/R/Tamb))*((p_amb*yH2O_cool)./pvap)));
            NH2O_solid = adsorbentMass*qH2O;
            % fluid phase H2O, final condition 
            NH2O_fluid = yH2O_cool*p_amb*V*void/(R*Tamb)*1e6; % mol
            % solid and liquid H2O, final condition  
            N_H2O_cool = NH2O_solid + NH2O_fluid;

    % solid and fluid at current point 
            pvap = ((22.064 * exp((647.096/Tamb)*((-7.85951783)*(1-(Tamb/647.096))+(1.84408259)*(1-(Tamb/647.096))^(1.5)+(-11.7866497)*(1-(Tamb/647.096))^3+(22.6807411)*(1-(Tamb/647.096))^(3.5)+(-15.9618719)*(1-(Tamb/647.096))^4+(1.80122502)*(1-(Tamb/647.096))^(7.5))))); % MPa
            qH2O = (Cm0*exp(beta/Tamb))*(CG0*exp(HC/R/Tamb))*(K0*exp(HK/R/Tamb))*((p_amb*(x(2)*(ub(2)-lb(2))+lb(2)))./pvap)./((1-(K0*exp(HK/R/Tamb))*((p_amb*(x(2)*(ub(2)-lb(2))+lb(2)))./pvap)).*(1+((CG0*exp(HC/R/Tamb))-1)*(K0*exp(HK/R/Tamb))*((p_amb*(x(2)*(ub(2)-lb(2))+lb(2)))./pvap)));
            NH2O_solid = adsorbentMass*qH2O;
            % fluid phase H2O, final condition 
            NH2O_fluid = (x(2)*(ub(2)-lb(2))+lb(2))*p_amb*V*void/(R*Tamb)*1e6; % mol
            % solid and fluid H2O, final condition  
            N_H2O_end = NH2O_solid + NH2O_fluid;

    fvec(2) = N_H2O_cool + yH2O_feed*(x(1)*(ub(1)-lb(1))+lb(1)) - N_H2O_end - yH2O_cool * (x(4)*(ub(4)-lb(4))+lb(4)); % = 0

    %% equation 3: material balance N2
    % N2 solid and fluid from previous step
        NN2_fluid = yN2_cool*p_amb*V*void/(R*Tamb)*1e6; % mol
        N_N2_cool = NN2_fluid;

    % fluid at current step 
        % fluid phase N2, final condition 
        NN2_fluid = (x(3)*(ub(3)-lb(3))+lb(3))*p_amb*V*void/(R*Tamb)*1e6; % mol
        N_N2_end = NN2_fluid;

    fvec(3) = N_N2_cool + yN2_feed*(x(1)*(ub(1)-lb(1))+lb(1)) - N_N2_end - yN2_cool * (x(4)*(ub(4)-lb(4))+lb(4)); % = 0

    %% equation 4: overall concentration balance
   fvec(4) = (x(2)*(ub(2)-lb(2))+lb(2)) + (x(3)*(ub(3)-lb(3))+lb(3)) + yCO2_sat - 1; % = 0

    fvec = real(fvec);
end

function fvec = adsorption2(x,adsorbentMass,void,V,R,yCO2_sat,yCO2_k,yN2_k,yH2O_k,lb,ub)
    % x(1): N_feed, x(3): yN2, x(4): N_waste

    %% equation 1: material balance CO2
    % CO2 solid and fluid from cooling
        % solid phase CO2, final condition 
        switch CO2IsothermModel
            case 'toth_cp'
                qCO2 = ((ns0_c*exp(Xi_c*(1-T/T0)).*(b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p_amb*yCO2_k))./((1+((b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p_amb*yCO2_k)).^(t0_c+alpha_c*(1-T0./T))).^(1./(t0_c+alpha_c*(1-T0./T)))))+((ns0_p*exp(Xi_p*(1-Tamb/T0)).*(b0_p*exp(dH_p/(R*T0)*(T0./Tamb-1))).*(p_amb*yCO2_k))./((1+((b0_p*exp(dH_p/(R*T0)*(T0./Tamb-1))).*(p_amb*yCO2_k)).^(t0_p+alpha_p*(1-T0./Tamb))).^(1./(t0_p+alpha_p*(1-T0./Tamb))))); % adsorbed amount of CO2
            case 's_shaped'            
                qCO2 = (q_L0.*(b_L0.*exp(dU_L./(R*T))).*(yCO2_k*p_amb)./(1+(b_L0.*exp(dU_L./(R*T))).*(yCO2_k*p_amb))).*(1-((exp((log((yCO2_k*p_amb))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((yCO2_k*p_amb)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*T))).*(yCO2_k*p_amb)./(1+(b_U0.*exp(dU_U./(R*T))).*(yCO2_k*p_amb))+(b_H0.*exp(dU_H./(R*T))).*(yCO2_k*p_amb)).*((exp((log((yCO2_k*p_amb))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((yCO2_k*p_amb)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam); 
            case 'DSL'
                qCO2 = n1*b0*exp(Hb/R/T)*(yCO2_k*p_amb)./(1e-6*R*T)./(1+b0*exp(Hb/R/T).*(yCO2_k*p_amb)./(1e-6*R*T)) + n2*d0*exp(Hd/R/T)*(yCO2_k*p_amb)./(1e-6*R*T)./(1+d0*exp(Hd/R/T).*(yCO2_k*p_amb)./(1e-6*R*T));            
            case 'DSL2'
                qCO2 = n1*b0*exp(Hb/R/T)*(yCO2_k*p_amb)./(1+b0*exp(Hb/R/T).*(yCO2_k*p_amb)) + n2*d0*exp(Hd/R/T)*(yCO2_k*p_amb)./(1+d0*exp(Hd/R/T).*(yCO2_k*p_amb));            
             case 'toth'
                qCO2 = ((ns0*exp(Xi*(1-T/T0)).*(b0*exp(dH/(R*T0)*(T0./T-1))).*(p_amb*yCO2_k))./((1+((b0*exp(dH/(R*T0)*(T0./T-1))).*(p_amb*yCO2_k)).^(t0+alpha*(1-T0./T))).^(1./(t0+alpha*(1-T0./T)))));           
            case 'langfr'
                qCO2 = (ns0.*exp(Xi.*(1-T./T0))).*((b0*exp(dH./(R*T0).*(T0./T-1))).*(p_amb*yCO2_k)).^(1./t0+alpha.*(1-T0./T))./(1+ ((b0*exp(dH./(R*T0).*(T0./T-1))).*(p_amb*yCO2_k)).^(1./t0+alpha.*(1-T0./T)));
        end        
        NCO2_solid = adsorbentMass*qCO2;
        % fluid phase CO2, final condition  
        NCO2_fluid = yCO2_k*p_amb*V*void/(R*T)*1e6; % mol
        % solid and fluid CO2, final condition  
        N_CO2_k = NCO2_solid + NCO2_fluid;

    % solid and fluid at end
        % solid phase CO2, final condition 
        switch CO2IsothermModel
            case 'toth_cp'
                qCO2 = ((ns0_c*exp(Xi_c*(1-T/T0)).*(b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p_amb*yCO2_sat))./((1+((b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p_amb*yCO2_sat)).^(t0_c+alpha_c*(1-T0./T))).^(1./(t0_c+alpha_c*(1-T0./T)))))+((ns0_p*exp(Xi_p*(1-Tamb/T0)).*(b0_p*exp(dH_p/(R*T0)*(T0./Tamb-1))).*(p_amb*yCO2_sat))./((1+((b0_p*exp(dH_p/(R*T0)*(T0./Tamb-1))).*(p_amb*yCO2_sat)).^(t0_p+alpha_p*(1-T0./Tamb))).^(1./(t0_p+alpha_p*(1-T0./Tamb))))); % adsorbed amount of CO2
            case 's_shaped'            
                qCO2 = (q_L0.*(b_L0.*exp(dU_L./(R*T))).*(yCO2_sat*p_amb)./(1+(b_L0.*exp(dU_L./(R*T))).*(yCO2_sat*p_amb))).*(1-((exp((log((yCO2_sat*p_amb))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((yCO2_sat*p_amb)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*T))).*(yCO2_sat*p_amb)./(1+(b_U0.*exp(dU_U./(R*T))).*(yCO2_sat*p_amb))+(b_H0.*exp(dU_H./(R*T))).*(yCO2_sat*p_amb)).*((exp((log((yCO2_sat*p_amb))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((yCO2_sat*p_amb)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam); 
            case 'DSL'
                qCO2 = n1*b0*exp(Hb/R/T)*(yCO2_sat*p_amb)./(1e-6*R*T)./(1+b0*exp(Hb/R/T).*(yCO2_sat*p_amb)./(1e-6*R*T)) + n2*d0*exp(Hd/R/T)*(yCO2_sat*p_amb)./(1e-6*R*T)./(1+d0*exp(Hd/R/T).*(yCO2_sat*p_amb)./(1e-6*R*T));            
            case 'DSL2'
                qCO2 = n1*b0*exp(Hb/R/T)*(yCO2_sat*p_amb)./(1+b0*exp(Hb/R/T).*(yCO2_sat*p_amb)) + n2*d0*exp(Hd/R/T)*(yCO2_sat*p_amb)./(1+d0*exp(Hd/R/T).*(yCO2_sat*p_amb));            
            case 'toth'
                qCO2 = ((ns0*exp(Xi*(1-T/T0)).*(b0*exp(dH/(R*T0)*(T0./T-1))).*(p_amb*yCO2_sat))./((1+((b0*exp(dH/(R*T0)*(T0./T-1))).*(p_amb*yCO2_sat)).^(t0+alpha*(1-T0./T))).^(1./(t0+alpha*(1-T0./T)))));                       
            case 'langfr'
                qCO2 = (ns0.*exp(Xi.*(1-T./T0))).*((b0*exp(dH./(R*T0).*(T0./T-1))).*(p_amb*yCO2_sat)).^(1./t0+alpha.*(1-T0./T))./(1+ ((b0*exp(dH./(R*T0).*(T0./T-1))).*(p_amb*yCO2_sat)).^(1./t0+alpha.*(1-T0./T)));
        end        
        NCO2_solid = adsorbentMass*qCO2;
        % fluid phase CO2, final condition  
        NCO2_fluid = yCO2_sat*p_amb*V*void/(R*T)*1e6; % mol
        % solid and fluid CO2, final condition  
        N_CO2_end = NCO2_solid + NCO2_fluid;
    
    fvec(1) = N_CO2_k + yCO2_feed*(x(1)*(ub(1)-lb(1))+lb(1)) - N_CO2_end - yCO2_k * (x(3)*(ub(3)-lb(3))+lb(3)); % = 0

    %% equation 2: material balance H2O
    % H2O solid and fluid from cooling
            % solid phase H2O, final condition 
            pvap = ((22.064 * exp((647.096/Tamb)*((-7.85951783)*(1-(Tamb/647.096))+(1.84408259)*(1-(Tamb/647.096))^(1.5)+(-11.7866497)*(1-(Tamb/647.096))^3+(22.6807411)*(1-(Tamb/647.096))^(3.5)+(-15.9618719)*(1-(Tamb/647.096))^4+(1.80122502)*(1-(Tamb/647.096))^(7.5))))); % MPa
            qH2O = (Cm0*exp(beta/Tamb))*(CG0*exp(HC/R/Tamb))*(K0*exp(HK/R/Tamb))*((p_amb*yH2O_k)./pvap)./((1-(K0*exp(HK/R/Tamb))*((p_amb*yH2O_k)./pvap)).*(1+((CG0*exp(HC/R/Tamb))-1)*(K0*exp(HK/R/Tamb))*((p_amb*yH2O_k)./pvap)));
            NH2O_solid = adsorbentMass*qH2O;
            % fluid phase H2O, final condition 
            NH2O_fluid = yH2O_k*p_amb*V*void/(R*Tamb)*1e6; % mol
            % solid and liquid H2O, final condition  
            N_H2O_k = NH2O_solid + NH2O_fluid;

    % solid and fluid at current point
            % solid phase H2O, final condition 
            pvap = ((22.064 * exp((647.096/Tamb)*((-7.85951783)*(1-(Tamb/647.096))+(1.84408259)*(1-(Tamb/647.096))^(1.5)+(-11.7866497)*(1-(Tamb/647.096))^3+(22.6807411)*(1-(Tamb/647.096))^(3.5)+(-15.9618719)*(1-(Tamb/647.096))^4+(1.80122502)*(1-(Tamb/647.096))^(7.5))))); % MPa
            qH2O = (Cm0*exp(beta/Tamb))*(CG0*exp(HC/R/Tamb))*(K0*exp(HK/R/Tamb))*((p_amb*yH2O_sat)./pvap)./((1-(K0*exp(HK/R/Tamb))*((p_amb*yH2O_sat)./pvap)).*(1+((CG0*exp(HC/R/Tamb))-1)*(K0*exp(HK/R/Tamb))*((p_amb*yH2O_sat)./pvap)));
            NH2O_solid = adsorbentMass*qH2O;
            % fluid phase H2O, final condition 
            NH2O_fluid = yH2O_sat*p_amb*V*void/(R*Tamb)*1e6; % mol
            % solid and fluid H2O, final condition  
            N_H2O_end = NH2O_solid + NH2O_fluid;

    %% equation 3: material balance N2
        NN2_fluid = yN2_k*p_amb*V*void/(R*Tamb)*1e6; % mol
        N_N2_k = NN2_fluid;

    % fluid at current step (with T=Tamb, p=pamb, yCO2=yCO2_feed)
        % fluid phase N2, final condition 
        NN2_fluid = (x(2)*(ub(2)-lb(2))+lb(2))*p_amb*V*void/(R*Tamb)*1e6; % mol
        N_N2_end = NN2_fluid;

    fvec(2) = N_N2_k + yN2_feed*(x(1)*(ub(1)-lb(1))+lb(1)) - N_N2_end - yN2_k * (x(3)*(ub(3)-lb(3))+lb(3)); % = 0

    %% equation 4: overall concentration balance
   fvec(3) = yCO2_sat + (x(2)*(ub(2)-lb(2))+lb(2)) + yH2O_sat - 1; % = 0

    fvec = real(fvec);
end


function fvec = adsorption2_dry(x,adsorbentMass,void,V,R,yCO2_sat,yCO2_k,yN2_k,lb,ub)
    % x(1): N_feed, x(3): yN2, x(4): N_waste
    %% equation 1: material balance CO2
    % CO2 solid and fluid from cooling
        % solid phase CO2, final condition 
        switch CO2IsothermModel
            case 'toth_cp'
                qCO2 = ((ns0_c*exp(Xi_c*(1-T/T0)).*(b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p_amb*yCO2_k))./((1+((b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p_amb*yCO2_k)).^(t0_c+alpha_c*(1-T0./T))).^(1./(t0_c+alpha_c*(1-T0./T)))))+((ns0_p*exp(Xi_p*(1-Tamb/T0)).*(b0_p*exp(dH_p/(R*T0)*(T0./Tamb-1))).*(p_amb*yCO2_k))./((1+((b0_p*exp(dH_p/(R*T0)*(T0./Tamb-1))).*(p_amb*yCO2_k)).^(t0_p+alpha_p*(1-T0./Tamb))).^(1./(t0_p+alpha_p*(1-T0./Tamb))))); % adsorbed amount of CO2
            case 's_shaped'            
                qCO2 = (q_L0.*(b_L0.*exp(dU_L./(R*T))).*(yCO2_k*p_amb)./(1+(b_L0.*exp(dU_L./(R*T))).*(yCO2_k*p_amb))).*(1-((exp((log((yCO2_k*p_amb))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((yCO2_k*p_amb)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*T))).*(yCO2_k*p_amb)./(1+(b_U0.*exp(dU_U./(R*T))).*(yCO2_k*p_amb))+(b_H0.*exp(dU_H./(R*T))).*(yCO2_k*p_amb)).*((exp((log((yCO2_k*p_amb))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((yCO2_k*p_amb)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam); 
            case 'DSL'
                qCO2 = n1*b0*exp(Hb/R/T)*(yCO2_k*p_amb)./(1e-6*R*T)./(1+b0*exp(Hb/R/T).*(yCO2_k*p_amb)./(1e-6*R*T)) + n2*d0*exp(Hd/R/T)*(yCO2_k*p_amb)./(1e-6*R*T)./(1+d0*exp(Hd/R/T).*(yCO2_k*p_amb)./(1e-6*R*T));            
            case 'DSL2'
                qCO2 = n1*b0*exp(Hb/R/T)*(yCO2_k*p_amb)./(1+b0*exp(Hb/R/T).*(yCO2_k*p_amb)) + n2*d0*exp(Hd/R/T)*(yCO2_k*p_amb)./(1+d0*exp(Hd/R/T).*(yCO2_k*p_amb));            
             case 'toth'
                qCO2 = ((ns0*exp(Xi*(1-T/T0)).*(b0*exp(dH/(R*T0)*(T0./T-1))).*(p_amb*yCO2_k))./((1+((b0*exp(dH/(R*T0)*(T0./T-1))).*(p_amb*yCO2_k)).^(t0+alpha*(1-T0./T))).^(1./(t0+alpha*(1-T0./T)))));           
            case 'langfr'
                qCO2 = (ns0.*exp(Xi.*(1-T./T0))).*((b0*exp(dH./(R*T0).*(T0./T-1))).*(p_amb*yCO2_k)).^(1./t0+alpha.*(1-T0./T))./(1+ ((b0*exp(dH./(R*T0).*(T0./T-1))).*(p_amb*yCO2_k)).^(1./t0+alpha.*(1-T0./T)));
        end        
        NCO2_solid = adsorbentMass*qCO2;
        % fluid phase CO2, final condition  
        NCO2_fluid = yCO2_k*p_amb*V*void/(R*T)*1e6; % mol
        % solid and fluid CO2, final condition  
        N_CO2_k = NCO2_solid + NCO2_fluid;

    % solid and fluid at end
        % solid phase CO2, final condition 
        switch CO2IsothermModel
            case 'toth_cp'
                qCO2 = ((ns0_c*exp(Xi_c*(1-T/T0)).*(b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p_amb*yCO2_sat))./((1+((b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p_amb*yCO2_sat)).^(t0_c+alpha_c*(1-T0./T))).^(1./(t0_c+alpha_c*(1-T0./T)))))+((ns0_p*exp(Xi_p*(1-Tamb/T0)).*(b0_p*exp(dH_p/(R*T0)*(T0./Tamb-1))).*(p_amb*yCO2_sat))./((1+((b0_p*exp(dH_p/(R*T0)*(T0./Tamb-1))).*(p_amb*yCO2_sat)).^(t0_p+alpha_p*(1-T0./Tamb))).^(1./(t0_p+alpha_p*(1-T0./Tamb))))); % adsorbed amount of CO2
            case 's_shaped'            
                qCO2 = (q_L0.*(b_L0.*exp(dU_L./(R*T))).*(yCO2_sat*p_amb)./(1+(b_L0.*exp(dU_L./(R*T))).*(yCO2_sat*p_amb))).*(1-((exp((log((yCO2_sat*p_amb))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((yCO2_sat*p_amb)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*T))).*(yCO2_sat*p_amb)./(1+(b_U0.*exp(dU_U./(R*T))).*(yCO2_sat*p_amb))+(b_H0.*exp(dU_H./(R*T))).*(yCO2_sat*p_amb)).*((exp((log((yCO2_sat*p_amb))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((yCO2_sat*p_amb)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam); 
            case 'DSL'
                qCO2 = n1*b0*exp(Hb/R/T)*(yCO2_sat*p_amb)./(1e-6*R*T)./(1+b0*exp(Hb/R/T).*(yCO2_sat*p_amb)./(1e-6*R*T)) + n2*d0*exp(Hd/R/T)*(yCO2_sat*p_amb)./(1e-6*R*T)./(1+d0*exp(Hd/R/T).*(yCO2_sat*p_amb)./(1e-6*R*T));            
            case 'DSL2'
                qCO2 = n1*b0*exp(Hb/R/T)*(yCO2_sat*p_amb)./(1+b0*exp(Hb/R/T).*(yCO2_sat*p_amb)) + n2*d0*exp(Hd/R/T)*(yCO2_sat*p_amb)./(1+d0*exp(Hd/R/T).*(yCO2_sat*p_amb));            
            case 'toth'
                qCO2 = ((ns0*exp(Xi*(1-T/T0)).*(b0*exp(dH/(R*T0)*(T0./T-1))).*(p_amb*yCO2_sat))./((1+((b0*exp(dH/(R*T0)*(T0./T-1))).*(p_amb*yCO2_sat)).^(t0+alpha*(1-T0./T))).^(1./(t0+alpha*(1-T0./T)))));                       
            case 'langfr'
                qCO2 = (ns0.*exp(Xi.*(1-T./T0))).*((b0*exp(dH./(R*T0).*(T0./T-1))).*(p_amb*yCO2_sat)).^(1./t0+alpha.*(1-T0./T))./(1+ ((b0*exp(dH./(R*T0).*(T0./T-1))).*(p_amb*yCO2_sat)).^(1./t0+alpha.*(1-T0./T)));
        end        
        NCO2_solid = adsorbentMass*qCO2;
        % fluid phase CO2, final condition  
        NCO2_fluid = yCO2_sat*p_amb*V*void/(R*T)*1e6; % mol
        % solid and fluid CO2, final condition  
        N_CO2_end = NCO2_solid + NCO2_fluid;
    
    fvec(1) = N_CO2_k + yCO2_feed*(x(1)*(ub(1)-lb(1))+lb(1)) - N_CO2_end - yCO2_k * (x(3)*(ub(3)-lb(3))+lb(3)); % = 0

    %% equation 3: material balance N2
    % unbekannt: N_feed, N_waste  
    % N2 solid and fluid from previous step
        % fluid phase N2, final condition 
        NN2_fluid = yN2_k*p_amb*V*void/(R*Tamb)*1e6; % mol
        N_N2_k = NN2_fluid;

    % fluid at current step
        % fluid phase N2, final condition 
        NN2_fluid = (x(2)*(ub(2)-lb(2))+lb(2))*p_amb*V*void/(R*Tamb)*1e6; % mol
        N_N2_end = NN2_fluid;

    fvec(2) = N_N2_k + yN2_feed*(x(1)*(ub(1)-lb(1))+lb(1)) - N_N2_end - yN2_k * (x(3)*(ub(3)-lb(3))+lb(3)); % = 0

    %% equation 4: overall concentration balance
   fvec(3) = yCO2_sat + (x(2)*(ub(2)-lb(2))+lb(2)) - 1; % = 0

    fvec = real(fvec);
end

end