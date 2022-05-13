
%% COOLING STEP

function outputCool = simulateCooling(dataModel,outputBlowEvac)
    
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
    
    %% define general data
    adsorbentMass = dataModel.process.adsorbentMass;
    R = dataModel.general.gasconstant;
    void = dataModel.process.voidFraction;
    V = dataModel.process.Vol;   
    
    %% RUN COOLING MODEL

    % feed stream
    p = dataModel.process.pamb; %
    yCO2_feed = dataModel.feed.yCO2;
    yH2O_feed = dataModel.feed.yH2O;
    yN2_feed = dataModel.feed.yN2;
	coolingTime = dataModel.process.noSteps;
    TempVector = dataModel.TCoolProfileStep;   

	T_k = dataModel.process.Tdes;   
    p_k = p; %
	pressure_vector = dataModel.PressureCoolProfileStep;	
	
    % Initial condition 
	if dataModel.feed.yH2O > 1e-10
		% (1) yCO2, (2), yN2, y(3) yH2O, (4) Nin
		x0 = [outputBlowEvac(end,1),outputBlowEvac(end,2),outputBlowEvac(end,3),0];
		yCO2_k = x0(1);
		yN2_k = x0(2);
		yH2O_k = x0(3);
		
		% boundaries
		lb = [0,0,0,0];
		ub = [1,1,1,50];
	else
		% (1) yCO2, (2), yN2, (4) Nin
		x0 = [outputBlowEvac(end,1),outputBlowEvac(end,2),0];
		yCO2_k = x0(1);
		yN2_k = x0(2);
		
		% boundaries
		lb = [0,0,0];
		ub = [1,1,50];
	end   

    % create vectors for saving results
    yCO2_save = zeros(coolingTime,1); 
    T_save = zeros(coolingTime,1); 
    yN2_save = zeros(coolingTime,1); 
    yH2O_save = zeros(coolingTime,1);  
    Nin_save = zeros(coolingTime,1); 
    if dataModel.feed.yH2O > 1e-10
		x_save = nan(coolingTime,4);
	else
		x_save = nan(coolingTime,3);
	end
    RH_save = zeros(coolingTime,1); 
    Nin_save_sum = zeros(coolingTime,1); 
    Nin_save_CO2 = zeros(coolingTime,1); 
    Nin_save_N2 = zeros(coolingTime,1); 
    Nin_save_H2O = zeros(coolingTime,1); 
    Q_save_sum = zeros(coolingTime,1);
    Ntotal_CO2 = zeros(coolingTime,1);
    Ntotal_N2 = zeros(coolingTime,1);
    Ntotal_H2O = zeros(coolingTime,1);
	
	% solver options
    options = optimoptions('lsqnonlin','Display','off');
    options.Algorithm = 'trust-region-reflective';
    options.OptimalityTolerance = 1e-9;
    options.FunctionTolerance = 1e-9;
    options.StepTolerance = 1e-9;
    options.MaxFunctionEvaluations = 300;  

    for k = 1:coolingTime

        T = TempVector(k);
		p = pressure_vector(k+1);
		p_k = pressure_vector(k);

		if dataModel.feed.yH2O > 1e-10
			funct = @(x) double( fcn(x));
		else
			funct = @(x) double( fcn_dry(x));
		end
        fun = @(x) funct(x);
		
		[x,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,lb,ub,options);

        yCO2_save(k) = x(1);
        yN2_save(k) = x(2);
		if dataModel.feed.yH2O > 1e-10
			yH2O_save(k) = x(3);
			Nin_save(k) = x(4);
			
			RH_save(k) = p*x(3)/vaporPressure(T)*100;
			if RH_save(k) > 100
				sprintf('WARNING COOL: Humidity is higher than 100%% (%0.6f)!', RH_save(k))
			end
			yCO2_k = x(1);
			yN2_k = x(2);
			yH2O_k = x(3);
		else
			Nin_save(k) = x(3);
			yCO2_k = x(1);
			yN2_k = x(2);
		end
		
        T_save(k) = T;
		x_save(k,:) = [residual];
       
        Nin_save_sum(k+1) = Nin_save_sum(k) + Nin_save(k);
        Nin_save_CO2(k+1) = Nin_save_CO2(k) + Nin_save(k)*x(1);
        Nin_save_N2(k+1) = Nin_save_N2(k) + Nin_save(k)*x(2);
        if dataModel.feed.yH2O > 1e-10
			Nin_save_H2O(k+1) = Nin_save_H2O(k) + Nin_save(k)*x(3);  
		else
			Nin_save_H2O(k+1) = 0; 
		end
        
        [Q, NCO2, NH2O, NN2, qCO2s, qH2O] = Qcool(x,dataModel.process.adsorbentMass,dataModel.sorbent(i).cp,dataModel.general.gasconstant);    
        Q_save_sum(k+1) =  Q_save_sum(k) + Q; %
		
		qCO2_save(k) = qCO2s;

        T_k = T;
        x0 = x; 
        
        Ntotal_CO2(k) = NCO2;
        Ntotal_N2(k) = NN2;
        Ntotal_H2O(k) = NH2O;
     end
    
    delta_CO2 = yCO2_save(end) - dataModel.feed.yCO2; % 
    if delta_CO2 > 0
        sprintf('WARNING: concentration at the end of the cooling step higher than CO2 composition in feed (%0.6f)!', delta_CO2)
    end
    
    %% output
    outputCool = [yCO2_save, yN2_save, yH2O_save, Nin_save, Nin_save_CO2(2:end), Nin_save_N2(2:end), Nin_save_H2O(2:end),T_save, Q_save_sum(2:end), Ntotal_CO2, Ntotal_N2, Ntotal_H2O];

%% functions

function fvec = fcn(x)
	% x(1): yCO2; x(2): yN2; x(3): yH2O; x(4): Nin
	pvap_k = ((22.064 * exp((647.096/T_k)*((-7.85951783)*(1-(T_k/647.096))+(1.84408259)*(1-(T_k/647.096))^(1.5)+(-11.7866497)*(1-(T_k/647.096))^3+(22.6807411)*(1-(T_k/647.096))^(3.5)+(-15.9618719)*(1-(T_k/647.096))^4+(1.80122502)*(1-(T_k/647.096))^(7.5))))); % MPa
	pvap = ((22.064 * exp((647.096/T)*((-7.85951783)*(1-(T/647.096))+(1.84408259)*(1-(T/647.096))^(1.5)+(-11.7866497)*(1-(T/647.096))^3+(22.6807411)*(1-(T/647.096))^(3.5)+(-15.9618719)*(1-(T/647.096))^4+(1.80122502)*(1-(T/647.096))^(7.5))))); % MPa

    %% equation 1: material balance CO2
    % CO2
        % solid phase CO2, initial condition
            switch CO2IsothermModel
                case 'toth_cp'
                    qCO2_k = ((ns0_c*exp(Xi_c*(1-T_k/T0)).*(b0_c*exp(dH_c/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k))./((1+((b0_c*exp(dH_c/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k)).^(t0_c+alpha_c*(1-T0./T_k))).^(1./(t0_c+alpha_c*(1-T0./T_k)))))+((ns0_p*exp(Xi_p*(1-T_k/T0)).*(b0_p*exp(dH_p/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k))./((1+((b0_p*exp(dH_p/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k)).^(t0_p+alpha_p*(1-T0./T_k))).^(1./(t0_p+alpha_p*(1-T0./T_k))))); % adsorbed amount of CO2
                case 's_shaped'
                    qCO2_k = (q_L0.*(b_L0.*exp(dU_L./(R*T_k))).*(yCO2_k*p_k)./(1+(b_L0.*exp(dU_L./(R*T_k))).*(yCO2_k*p_k))).*(1-((exp((log((yCO2_k*p_k))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k))))./(1+exp(((log((yCO2_k*p_k)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*T_k))).*(yCO2_k*p_k)./(1+(b_U0.*exp(dU_U./(R*T_k))).*(yCO2_k*p_k))+(b_H0.*exp(dU_H./(R*T_k))).*(yCO2_k*p_k)).*((exp((log((yCO2_k*p_k))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k))))./(1+exp(((log((yCO2_k*p_k)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k)))))).^gam); 
                case 'DSL'
                    qCO2_k = n1*b0*exp(Hb/R/T_k)*(yCO2_k*p_k)./(1e-6*R*T)./(1+b0*exp(Hb/R/T_k).*(yCO2_k*p_k)./(1e-6*R*T)) + n2*d0*exp(Hd/R/T_k)*(yCO2_k*p_k)./(1e-6*R*T)./(1+d0*exp(Hd/R/T_k).*(yCO2_k*p_k)./(1e-6*R*T));  
				case 'DSL2'
                    qCO2_k = n1*b0*exp(Hb/R/T_k)*(yCO2_k*p_k)./(1+b0*exp(Hb/R/T_k).*(yCO2_k*p_k)) + n2*d0*exp(Hd/R/T_k)*(yCO2_k*p_k)./(1+d0*exp(Hd/R/T_k).*(yCO2_k*p_k));            
                case 'toth'
                    qCO2_k = ((ns0*exp(Xi*(1-T_k/T0)).*(b0*exp(dH/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k))./((1+((b0*exp(dH/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k)).^(t0+alpha*(1-T0./T_k))).^(1./(t0+alpha*(1-T0./T_k)))));            
				case 'langfr'
					qCO2_k = (ns0.*exp(Xi.*(1-T_k./T0))).*((b0*exp(dH./(R*T0).*(T0./T_k-1))).*(p_k*yCO2_k)).^(1./t0+alpha.*(1-T0./T_k))./(1+ ((b0*exp(dH./(R*T0).*(T0./T_k-1))).*(p_k*yCO2_k)).^(1./t0+alpha.*(1-T0./T_k)));
			end
        NCO2_solid_k = adsorbentMass*qCO2_k;
        % fluid phase CO2, initial condition
        NCO2_fluid_k = yCO2_k*p_k*V*void/(R*T_k)*1e6; % mol
        % solid and liquid CO2, initial condition    
        NCO2_k = NCO2_solid_k + NCO2_fluid_k;

        % solid phase CO2, final condition 
            switch CO2IsothermModel
                case 'toth_cp'
                    qCO2 = ((ns0_c*exp(Xi_c*(1-T/T0)).*(b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p*x(1)))./((1+((b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p*x(1))).^(t0_c+alpha_c*(1-T0./T))).^(1./(t0_c+alpha_c*(1-T0./T)))))+((ns0_p*exp(Xi_p*(1-T/T0)).*(b0_p*exp(dH_p/(R*T0)*(T0./T-1))).*(p*x(1)))./((1+((b0_p*exp(dH_p/(R*T0)*(T0./T-1))).*(p*x(1))).^(t0_p+alpha_p*(1-T0./T))).^(1./(t0_p+alpha_p*(1-T0./T))))); % adsorbed amount of CO2
                case 's_shaped'         
                    qCO2 = (q_L0.*(b_L0.*exp(dU_L./(R*T))).*(x(1)*p)./(1+(b_L0.*exp(dU_L./(R*T))).*(x(1)*p))).*(1-((exp((log((x(1)*p))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((x(1)*p)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*T))).*(x(1)*p)./(1+(b_U0.*exp(dU_U./(R*T))).*(x(1)*p))+(b_H0.*exp(dU_H./(R*T))).*(x(1)*p)).*((exp((log((x(1)*p))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((x(1)*p)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam); % adsorbed amount of CO2                    
                case 'DSL'
                    qCO2 = n1*b0*exp(Hb/R/T)*(x(1)*p)./(1e-6*R*T)./(1+b0*exp(Hb/R/T).*(x(1)*p)./(1e-6*R*T)) + n2*d0*exp(Hd/R/T)*(x(1)*p)./(1e-6*R*T)./(1+d0*exp(Hd/R/T).*(x(1)*p)./(1e-6*R*T)); 
				case 'DSL2'
                    qCO2 = n1*b0*exp(Hb/R/T)*(x(1)*p)./(1+b0*exp(Hb/R/T).*(x(1)*p)) + n2*d0*exp(Hd/R/T)*(x(1)*p)./(1+d0*exp(Hd/R/T).*(x(1)*p));            
                case 'toth'
                    qCO2 = ((ns0*exp(Xi*(1-T/T0)).*(b0*exp(dH/(R*T0)*(T0./T-1))).*(x(1)*p))./((1+((b0*exp(dH/(R*T0)*(T0./T-1))).*(x(1)*p)).^(t0+alpha*(1-T0./T))).^(1./(t0+alpha*(1-T0./T)))));			
				case 'langfr'
					qCO2 = (ns0.*exp(Xi.*(1-T./T0))).*((b0*exp(dH./(R*T0).*(T0./T-1))).*(x(1)*p)).^(1./t0+alpha.*(1-T0./T))./(1+ ((b0*exp(dH./(R*T0).*(T0./T-1))).*(x(1)*p)).^(1./t0+alpha.*(1-T0./T)));
			end
        NCO2_solid = adsorbentMass*qCO2;
        % fluid phase CO2, final condition  
        NCO2_fluid = x(1)*p*V*void/(R*T)*1e6; % mol
        % solid and liquid CO2, final condition  
        NCO2 = NCO2_solid + NCO2_fluid;

        % material balance CO2
        fvec(1) = - NCO2 + yCO2_feed*x(4) + NCO2_k; % == 0

    %% equation 2: material balance H2O
    % H2O
        % solid phase H2O, initial condition
        qH2O_k = (Cm0*exp(beta/T_k))*(CG0*exp(HC/R/T_k))*(K0*exp(HK/R/T_k))*((p_k*yH2O_k)./pvap_k)./((1-(K0*exp(HK/R/T_k))*((p_k*yH2O_k)./pvap_k)).*(1+((CG0*exp(HC/R/T_k))-1)*(K0*exp(HK/R/T_k))*((p_k*yH2O_k)./pvap_k)));
        NH2O_solid_k = adsorbentMass*qH2O_k;  
        % fluid phase H2O, initial condition
        NH2O_fluid_k = yH2O_k*p_k*V*void/(R*T_k)*1e6; % mol
        % solid and liquid H2O, initial condition   
        NH2O_k = NH2O_solid_k + NH2O_fluid_k;

        % solid phase H2O, final condition 
        qH2O = (Cm0*exp(beta/T))*(CG0*exp(HC/R/T))*(K0*exp(HK/R/T))*((p*x(3))./pvap)./((1-(K0*exp(HK/R/T))*((p*x(3))./pvap)).*(1+((CG0*exp(HC/R/T))-1)*(K0*exp(HK/R/T))*((p*x(3))./pvap)));
        NH2O_solid = adsorbentMass*qH2O;
        % fluid phase H2O, final condition 
        NH2O_fluid = x(3)*p*V*void/(R*T)*1e6; % mol
        % solid and liquid H2O, final condition  
        NH2O = NH2O_solid + NH2O_fluid;

        % material balance H2O
        fvec(2) = - NH2O + yH2O_feed*x(4) + NH2O_k; % == 0

    %% equation 3 material balance N2
    % N2
        % fluid phase N2, initial condition
        NN2_k = yN2_k*p_k*V*void/(R*T_k)*1e6; % mol
        % fluid phase N2, final condition
        NN2 = x(2)*p*V*void/(R*T)*1e6; % mol

		fvec(3) =  - (NN2 + NH2O + NCO2) + x(4) + (NH2O_k + NCO2_k + NN2_k); % == 0

    %% equation 4: overall material balance
        fvec(4) = 1-x(1)-x(2)-x(3); % == 0
		
		fvec = real(fvec);
    end

    function [Q, NCO2, NH2O, NN2, qCO2, qH2O] = Qcool(x,adsorbentMass,cp,R)
	
       %% equation 5: energy balance
        % CO2
            % solid phase CO2, initial condition
            switch CO2IsothermModel
                case 'toth_cp'
                    qCO2_k = ((ns0_c*exp(Xi_c*(1-T_k/T0)).*(b0_c*exp(dH_c/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k))./((1+((b0_c*exp(dH_c/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k)).^(t0_c+alpha_c*(1-T0./T_k))).^(1./(t0_c+alpha_c*(1-T0./T_k)))))+((ns0_p*exp(Xi_p*(1-T_k/T0)).*(b0_p*exp(dH_p/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k))./((1+((b0_p*exp(dH_p/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k)).^(t0_p+alpha_p*(1-T0./T_k))).^(1./(t0_p+alpha_p*(1-T0./T_k))))); % adsorbed amount of CO2
                case 's_shaped'   
                    qCO2_k = (q_L0.*(b_L0.*exp(dU_L./(R*T_k))).*(yCO2_k*p_k)./(1+(b_L0.*exp(dU_L./(R*T_k))).*(yCO2_k*p_k))).*(1-((exp((log((yCO2_k*p_k))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k))))./(1+exp(((log((yCO2_k*p_k)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*T_k))).*(yCO2_k*p_k)./(1+(b_U0.*exp(dU_U./(R*T_k))).*(yCO2_k*p_k))+(b_H0.*exp(dU_H./(R*T_k))).*(yCO2_k*p_k)).*((exp((log((yCO2_k*p_k))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k))))./(1+exp(((log((yCO2_k*p_k)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k)))))).^gam); 
                case 'DSL'
                    qCO2_k = n1*b0*exp(Hb/R/T_k)*(yCO2_k*p_k)./(1e-6*R*T_k)./(1+b0*exp(Hb/R/T_k).*(yCO2_k*p_k)./(1e-6*R*T_k)) + n2*d0*exp(Hd/R/T_k)*(yCO2_k*p_k)./(1e-6*R*T_k)./(1+d0*exp(Hd/R/T_k).*(yCO2_k*p_k)./(1e-6*R*T_k));
				case 'DSL2'
                    qCO2_k = n1*b0*exp(Hb/R/T_k)*(yCO2_k*p_k)./(1+b0*exp(Hb/R/T_k).*(yCO2_k*p_k)) + n2*d0*exp(Hd/R/T_k)*(yCO2_k*p_k)./(1+d0*exp(Hd/R/T_k).*(yCO2_k*p_k));
                case 'toth'
                    qCO2_k = ((ns0*exp(Xi*(1-T_k/T0)).*(b0*exp(dH/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k))./((1+((b0*exp(dH/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k)).^(t0+alpha*(1-T0./T_k))).^(1./(t0+alpha*(1-T0./T_k)))));            
				case 'langfr'
					qCO2_k = (ns0.*exp(Xi.*(1-T_k./T0))).*((b0*exp(dH./(R*T0).*(T0./T_k-1))).*(p_k*yCO2_k)).^(1./t0+alpha.*(1-T0./T_k))./(1+ ((b0*exp(dH./(R*T0).*(T0./T_k-1))).*(p_k*yCO2_k)).^(1./t0+alpha.*(1-T0./T_k)));
			end         
            NCO2_solid_k = adsorbentMass*qCO2_k;

            % solid phase CO2, final condition 
            switch CO2IsothermModel
                case 'toth_cp'
                    qCO2 = ((ns0_c*exp(Xi_c*(1-T/T0)).*(b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p*x(1)))./((1+((b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p*x(1))).^(t0_c+alpha_c*(1-T0./T))).^(1./(t0_c+alpha_c*(1-T0./T)))))+((ns0_p*exp(Xi_p*(1-T/T0)).*(b0_p*exp(dH_p/(R*T0)*(T0./T-1))).*(p*x(1)))./((1+((b0_p*exp(dH_p/(R*T0)*(T0./T-1))).*(p*x(1))).^(t0_p+alpha_p*(1-T0./T))).^(1./(t0_p+alpha_p*(1-T0./T))))); % adsorbed amount of CO2                    
                case 's_shaped'   
                    qCO2 = (q_L0.*(b_L0.*exp(dU_L./(R*T))).*(x(1)*p)./(1+(b_L0.*exp(dU_L./(R*T))).*(x(1)*p))).*(1-((exp((log((x(1)*p))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((x(1)*p)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*T))).*(x(1)*p)./(1+(b_U0.*exp(dU_U./(R*T))).*(x(1)*p))+(b_H0.*exp(dU_H./(R*T))).*(x(1)*p)).*((exp((log((x(1)*p))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((x(1)*p)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam); % adsorbed amount of CO2                    
               case 'DSL'
                    qCO2 = n1*b0*exp(Hb/R/T)*(x(1)*p)./(1e-6*R*T)./(1+b0*exp(Hb/R/T).*(x(1)*p)./(1e-6*R*T)) + n2*d0*exp(Hd/R/T)*(x(1)*p)./(1e-6*R*T)./(1+d0*exp(Hd/R/T).*(x(1)*p)./(1e-6*R*T));
				case 'DSL2'
                    qCO2 = n1*b0*exp(Hb/R/T)*(x(1)*p)./(1+b0*exp(Hb/R/T).*(x(1)*p)) + n2*d0*exp(Hd/R/T)*(x(1)*p)./(1+d0*exp(Hd/R/T).*(x(1)*p));            
                case 'toth'
                    qCO2 = ((ns0*exp(Xi*(1-T/T0)).*(b0*exp(dH/(R*T0)*(T0./T-1))).*(x(1)*p))./((1+((b0*exp(dH/(R*T0)*(T0./T-1))).*(x(1)*p)).^(t0+alpha*(1-T0./T))).^(1./(t0+alpha*(1-T0./T)))));			            
				case 'langfr'
					qCO2 = (ns0.*exp(Xi.*(1-T./T0))).*((b0*exp(dH./(R*T0).*(T0./T-1))).*(x(1)*p)).^(1./t0+alpha.*(1-T0./T))./(1+ ((b0*exp(dH./(R*T0).*(T0./T-1))).*(x(1)*p)).^(1./t0+alpha.*(1-T0./T)));
			end
            NCO2_solid = adsorbentMass*qCO2;
            % fluid phase CO2, final condition  
            NCO2_fluid = x(1)*p*V*void/(R*T)*1e6; % mol
            % solid and liquid CO2, final condition  
            NCO2 = NCO2_solid + NCO2_fluid;


        % H2O
		if dataModel.feed.yH2O > 1e-10
			pvap_k = ((22.064 * exp((647.096/T_k)*((-7.85951783)*(1-(T_k/647.096))+(1.84408259)*(1-(T_k/647.096))^(1.5)+(-11.7866497)*(1-(T_k/647.096))^3+(22.6807411)*(1-(T_k/647.096))^(3.5)+(-15.9618719)*(1-(T_k/647.096))^4+(1.80122502)*(1-(T_k/647.096))^(7.5))))); % MPa
			pvap = ((22.064 * exp((647.096/T)*((-7.85951783)*(1-(T/647.096))+(1.84408259)*(1-(T/647.096))^(1.5)+(-11.7866497)*(1-(T/647.096))^3+(22.6807411)*(1-(T/647.096))^(3.5)+(-15.9618719)*(1-(T/647.096))^4+(1.80122502)*(1-(T/647.096))^(7.5))))); % MPa

            % solid phase H2O, initial condition
            qH2O_k = (Cm0*exp(beta/T_k))*(CG0*exp(HC/R/T_k))*(K0*exp(HK/R/T_k))*((p_k*yH2O_k)./pvap_k)./((1-(K0*exp(HK/R/T_k))*((p_k*yH2O_k)./pvap_k)).*(1+((CG0*exp(HC/R/T_k))-1)*(K0*exp(HK/R/T_k))*((p_k*yH2O_k)./pvap_k)));
            NH2O_solid_k = adsorbentMass*qH2O_k;  

            % solid phase H2O, final condition 
            qH2O = (Cm0*exp(beta/T))*(CG0*exp(HC/R/T))*(K0*exp(HK/R/T))*((p*x(3))./pvap)./((1-(K0*exp(HK/R/T))*((p*x(3))./pvap)).*(1+((CG0*exp(HC/R/T))-1)*(K0*exp(HK/R/T))*((p*x(3))./pvap)));
            NH2O_solid = adsorbentMass*qH2O;
            
            % fluid phase H2O, final condition 
            NH2O_fluid = x(3)*p*V*void/(R*T)*1e6; % mol
            % solid and liquid H2O, final condition  
            NH2O = NH2O_solid + NH2O_fluid;
		else
			qH2O = 0;
			NH2O = 0;
		end            
        % N2
            % fluid phase N2, final condition
            NN2 = x(2)*p*V*void/(R*T)*1e6; % mol

        %% equation 5: energy balance
        % CO2
            switch CO2IsothermModel
                case 'toth_cp'
					dqCO2_dTequ_k = - b0_c.*ns0_c.*(p_k.*yCO2_k).*exp((dH_c.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi_c.*(T_k./T0 - 1)).*(((T0.*alpha_c.*log(b0_c.*(p_k.*yCO2_k).*exp((dH_c.*(T0./T_k - 1))./(R.*T0))).*(b0_c.*(p_k.*yCO2_k).*exp((dH_c.*(T0./T_k - 1))./(R.*T0))).^(t0_c - alpha_c.*(T0./T_k - 1)))./T_k.^2 - (b0_c.*dH_c.*(p_k.*yCO2_k).*exp((dH_c.*(T0./T_k - 1))./(R.*T0)).*(t0_c - alpha_c.*(T0./T_k - 1)).*(b0_c.*(p_k.*yCO2_k).*exp((dH_c.*(T0./T_k - 1))./(R.*T0))).^(t0_c - alpha_c.*(T0./T_k - 1) - 1))./(R.*T_k.^2))./((t0_c - alpha_c.*(T0./T_k - 1)).*((b0_c.*(p_k.*yCO2_k).*exp((dH_c.*(T0./T_k - 1))./(R.*T0))).^(t0_c - alpha_c.*(T0./T_k - 1)) + 1).^(1./(t0_c - alpha_c.*(T0./T_k - 1)) + 1)) - (T0.*alpha_c.*log((b0_c.*(p_k.*yCO2_k).*exp((dH_c.*(T0./T_k - 1))./(R.*T0))).^(t0_c - alpha_c.*(T0./T_k - 1)) + 1))./(T_k.^2.*(t0_c - alpha_c.*(T0./T_k - 1)).^2.*((b0_c.*(p_k.*yCO2_k).*exp((dH_c.*(T0./T_k - 1))./(R.*T0))).^(t0_c - alpha_c.*(T0./T_k - 1)) + 1).^(1./(t0_c - alpha_c.*(T0./T_k - 1))))) - b0_p.*ns0_p.*(p_k.*yCO2_k).*exp((dH_p.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi_p.*(T_k./T0 - 1)).*(((T0.*alpha_p.*log(b0_p.*(p_k.*yCO2_k).*exp((dH_p.*(T0./T_k - 1))./(R.*T0))).*(b0_p.*(p_k.*yCO2_k).*exp((dH_p.*(T0./T_k - 1))./(R.*T0))).^(t0_p - alpha_p.*(T0./T_k - 1)))./T_k.^2 - (b0_p.*dH_p.*(p_k.*yCO2_k).*exp((dH_p.*(T0./T_k - 1))./(R.*T0)).*(t0_p - alpha_p.*(T0./T_k - 1)).*(b0_p.*(p_k.*yCO2_k).*exp((dH_p.*(T0./T_k - 1))./(R.*T0))).^(t0_p - alpha_p.*(T0./T_k - 1) - 1))./(R.*T_k.^2))./((t0_p - alpha_p.*(T0./T_k - 1)).*((b0_p.*(p_k.*yCO2_k).*exp((dH_p.*(T0./T_k - 1))./(R.*T0))).^(t0_p - alpha_p.*(T0./T_k - 1)) + 1).^(1./(t0_p - alpha_p.*(T0./T_k - 1)) + 1)) - (T0.*alpha_p.*log((b0_p.*(p_k.*yCO2_k).*exp((dH_p.*(T0./T_k - 1))./(R.*T0))).^(t0_p - alpha_p.*(T0./T_k - 1)) + 1))./(T_k.^2.*(t0_p - alpha_p.*(T0./T_k - 1)).^2.*((b0_p.*(p_k.*yCO2_k).*exp((dH_p.*(T0./T_k - 1))./(R.*T0))).^(t0_p - alpha_p.*(T0./T_k - 1)) + 1).^(1./(t0_p - alpha_p.*(T0./T_k - 1))))) - (Xi_c.*b0_c.*ns0_c.*(p_k.*yCO2_k).*exp((dH_c.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi_c.*(T_k./T0 - 1)))./(T0.*((b0_c.*(p_k.*yCO2_k).*exp((dH_c.*(T0./T_k - 1))./(R.*T0))).^(t0_c - alpha_c.*(T0./T_k - 1)) + 1).^(1./(t0_c - alpha_c.*(T0./T_k - 1)))) - (Xi_p.*b0_p.*ns0_p.*(p_k.*yCO2_k).*exp((dH_p.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi_p.*(T_k./T0 - 1)))./(T0.*((b0_p.*(p_k.*yCO2_k).*exp((dH_p.*(T0./T_k - 1))./(R.*T0))).^(t0_p - alpha_p.*(T0./T_k - 1)) + 1).^(1./(t0_p - alpha_p.*(T0./T_k - 1)))) - (b0_c.*dH_c.*ns0_c.*(p_k.*yCO2_k).*exp((dH_c.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi_c.*(T_k./T0 - 1)))./(R.*T_k.^2.*((b0_c.*(p_k.*yCO2_k).*exp((dH_c.*(T0./T_k - 1))./(R.*T0))).^(t0_c - alpha_c.*(T0./T_k - 1)) + 1).^(1./(t0_c - alpha_c.*(T0./T_k - 1)))) - (b0_p.*dH_p.*ns0_p.*(p_k.*yCO2_k).*exp((dH_p.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi_p.*(T_k./T0 - 1)))./(R.*T_k.^2.*((b0_p.*(p_k.*yCO2_k).*exp((dH_p.*(T0./T_k - 1))./(R.*T0))).^(t0_p - alpha_p.*(T0./T_k - 1)) + 1).^(1./(t0_p - alpha_p.*(T0./T_k - 1))));
					dqCO2_dpCO2_k = (b0_c.*ns0_c.*exp((dH_c.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi_c.*(T_k./T0 - 1)))./((b0_c.*(p_k.*yCO2_k).*exp((dH_c.*(T0./T_k - 1))./(R.*T0))).^(t0_c - alpha_c.*(T0./T_k - 1)) + 1).^(1./(t0_c - alpha_c.*(T0./T_k - 1))) + (b0_p.*ns0_p.*exp((dH_p.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi_p.*(T_k./T0 - 1)))./((b0_p.*(p_k.*yCO2_k).*exp((dH_p.*(T0./T_k - 1))./(R.*T0))).^(t0_p - alpha_p.*(T0./T_k - 1)) + 1).^(1./(t0_p - alpha_p.*(T0./T_k - 1))) - (b0_c.^2.*ns0_c.*(p_k.*yCO2_k).*exp((2.*dH_c.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi_c.*(T_k./T0 - 1)).*(b0_c.*(p_k.*yCO2_k).*exp((dH_c.*(T0./T_k - 1))./(R.*T0))).^(t0_c - alpha_c.*(T0./T_k - 1) - 1))./((b0_c.*(p_k.*yCO2_k).*exp((dH_c.*(T0./T_k - 1))./(R.*T0))).^(t0_c - alpha_c.*(T0./T_k - 1)) + 1).^(1./(t0_c - alpha_c.*(T0./T_k - 1)) + 1) - (b0_p.^2.*ns0_p.*(p_k.*yCO2_k).*exp((2.*dH_p.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi_p.*(T_k./T0 - 1)).*(b0_p.*(p_k.*yCO2_k).*exp((dH_p.*(T0./T_k - 1))./(R.*T0))).^(t0_p - alpha_p.*(T0./T_k - 1) - 1))./((b0_p.*(p_k.*yCO2_k).*exp((dH_p.*(T0./T_k - 1))./(R.*T0))).^(t0_p - alpha_p.*(T0./T_k - 1)) + 1).^(1./(t0_p - alpha_p.*(T0./T_k - 1)) + 1);
				case 's_shaped'  
					dqCO2_dTequ_k = gam.*((exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1).*((dH_step.*exp(-xi_2.*(1./T0 - 1./T_k)))./(R.*T_k.^2.*xi_1) + (xi_2.*exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./(T_k.^2.*xi_1)))./(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1) - (exp(-(2.*exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1).*((dH_step.*exp(-xi_2.*(1./T0 - 1./T_k)))./(R.*T_k.^2.*xi_1) + (xi_2.*exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./(T_k.^2.*xi_1)))./(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1).^2).*(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1)./(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1)).^(gam - 1).*(b_H0.*(p_k.*yCO2_k).*exp(dU_H./(R.*T_k)) + (b_U0.*(p_k.*yCO2_k).*q_U0.*exp(dU_U./(R.*T_k)))./(b_U0.*(p_k.*yCO2_k).*exp(dU_U./(R.*T_k)) + 1)) - (exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1)./(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1)).^gam.*((b_H0.*dU_H.*(p_k.*yCO2_k).*exp(dU_H./(R.*T_k)))./(R.*T_k.^2) - (b_U0.^2.*dU_U.*(p_k.*yCO2_k).^2.*q_U0.*exp((2.*dU_U)./(R.*T_k)))./(R.*T_k.^2.*(b_U0.*(p_k.*yCO2_k).*exp(dU_U./(R.*T_k)) + 1).^2) + (b_U0.*dU_U.*(p_k.*yCO2_k).*q_U0.*exp(dU_U./(R.*T_k)))./(R.*T_k.^2.*(b_U0.*(p_k.*yCO2_k).*exp(dU_U./(R.*T_k)) + 1))) - (b_L0.*gam.*(p_k.*yCO2_k).*q_L0.*exp(dU_L./(R.*T_k)).*((exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1).*((dH_step.*exp(-xi_2.*(1./T0 - 1./T_k)))./(R.*T_k.^2.*xi_1) + (xi_2.*exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./(T_k.^2.*xi_1)))./(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1) - (exp(-(2.*exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1).*((dH_step.*exp(-xi_2.*(1./T0 - 1./T_k)))./(R.*T_k.^2.*xi_1) + (xi_2.*exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./(T_k.^2.*xi_1)))./(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1).^2).*(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1)./(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1)).^(gam - 1))./(b_L0.*(p_k.*yCO2_k).*exp(dU_L./(R.*T_k)) + 1) - (b_L0.^2.*dU_L.*(p_k.*yCO2_k).^2.*q_L0.*exp((2.*dU_L)./(R.*T_k)).*((exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1)./(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1)).^gam - 1))./(R.*T_k.^2.*(b_L0.*(p_k.*yCO2_k).*exp(dU_L./(R.*T_k)) + 1).^2) + (b_L0.*dU_L.*(p_k.*yCO2_k).*q_L0.*exp(dU_L./(R.*T_k)).*((exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1)./(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1)).^gam - 1))./(R.*T_k.^2.*(b_L0.*(p_k.*yCO2_k).*exp(dU_L./(R.*T_k)) + 1));
					dqCO2_dpCO2_k =  (exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1)./(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1)).^gam.*(b_H0.*exp(dU_H./(R.*T_k)) + (b_U0.*q_U0.*exp(dU_U./(R.*T_k)))./(b_U0.*(p_k.*yCO2_k).*exp(dU_U./(R.*T_k)) + 1) - (b_U0.^2.*(p_k.*yCO2_k).*q_U0.*exp((2.*dU_U)./(R.*T_k)))./(b_U0.*(p_k.*yCO2_k).*exp(dU_U./(R.*T_k)) + 1).^2) + gam.*((exp(-xi_2.*(1./T0 - 1./T_k)).*exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1))./((p_k.*yCO2_k).*xi_1.*(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1)) - (exp(-xi_2.*(1./T0 - 1./T_k)).*exp(-(2.*exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1))./((p_k.*yCO2_k).*xi_1.*(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1).^2)).*(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1)./(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1)).^(gam - 1).*(b_H0.*(p_k.*yCO2_k).*exp(dU_H./(R.*T_k)) + (b_U0.*(p_k.*yCO2_k).*q_U0.*exp(dU_U./(R.*T_k)))./(b_U0.*(p_k.*yCO2_k).*exp(dU_U./(R.*T_k)) + 1)) - (b_L0.*q_L0.*exp(dU_L./(R.*T_k)).*((exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1)./(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1)).^gam - 1))./(b_L0.*(p_k.*yCO2_k).*exp(dU_L./(R.*T_k)) + 1) + (b_L0.^2.*(p_k.*yCO2_k).*q_L0.*exp((2.*dU_L)./(R.*T_k)).*((exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1)./(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1)).^gam - 1))./(b_L0.*(p_k.*yCO2_k).*exp(dU_L./(R.*T_k)) + 1).^2 - (b_L0.*gam.*(p_k.*yCO2_k).*q_L0.*exp(dU_L./(R.*T_k)).*((exp(-xi_2.*(1./T0 - 1./T_k)).*exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1))./((p_k.*yCO2_k).*xi_1.*(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1)) - (exp(-xi_2.*(1./T0 - 1./T_k)).*exp(-(2.*exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1))./((p_k.*yCO2_k).*xi_1.*(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1).^2)).*(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1)./(exp(-(exp(-xi_2.*(1./T0 - 1./T_k)).*(log(p_step0.*exp(-(dH_step.*(1./T0 - 1./T_k))./R)) - log((p_k.*yCO2_k))))./xi_1) + 1)).^(gam - 1))./(b_L0.*(p_k.*yCO2_k).*exp(dU_L./(R.*T_k)) + 1);			
				case 'DSL'
					dqCO2_dTequ_k = (1000000.*b0.*n1.*(yCO2_k.*p_k).*yCO2_k.*exp(Hb./(R.*T_k)).*((1000000.*b0.*(yCO2_k.*p_k).*yCO2_k.*exp(Hb./(R.*T_k)))./(R.*T_k.^2) + (1000000.*Hb.*b0.*(yCO2_k.*p_k).*yCO2_k.*exp(Hb./(R.*T_k)))./(R.^2.*T_k.^3)))./(R.*T_k.*((1000000.*b0.*(yCO2_k.*p_k).*yCO2_k.*exp(Hb./(R.*T_k)))./(R.*T_k) + 1).^2) - (1000000.*d0.*n2.*(yCO2_k.*p_k).*yCO2_k.*exp(Hd./(R.*T_k)))./(R.*T_k.^2.*((1000000.*d0.*(yCO2_k.*p_k).*yCO2_k.*exp(Hd./(R.*T_k)))./(R.*T_k) + 1)) - (1000000.*Hb.*b0.*n1.*(yCO2_k.*p_k).*yCO2_k.*exp(Hb./(R.*T_k)))./(R.^2.*T_k.^3.*((1000000.*b0.*(yCO2_k.*p_k).*yCO2_k.*exp(Hb./(R.*T_k)))./(R.*T_k) + 1)) - (1000000.*Hd.*d0.*n2.*(yCO2_k.*p_k).*yCO2_k.*exp(Hd./(R.*T_k)))./(R.^2.*T_k.^3.*((1000000.*d0.*(yCO2_k.*p_k).*yCO2_k.*exp(Hd./(R.*T_k)))./(R.*T_k) + 1)) - (1000000.*b0.*n1.*(yCO2_k.*p_k).*yCO2_k.*exp(Hb./(R.*T_k)))./(R.*T_k.^2.*((1000000.*b0.*(yCO2_k.*p_k).*yCO2_k.*exp(Hb./(R.*T_k)))./(R.*T_k) + 1)) + (1000000.*d0.*n2.*(yCO2_k.*p_k).*yCO2_k.*exp(Hd./(R.*T_k)).*((1000000.*d0.*(yCO2_k.*p_k).*yCO2_k.*exp(Hd./(R.*T_k)))./(R.*T_k.^2) + (1000000.*Hd.*d0.*(yCO2_k.*p_k).*yCO2_k.*exp(Hd./(R.*T_k)))./(R.^2.*T_k.^3)))./(R.*T_k.*((1000000.*d0.*(yCO2_k.*p_k).*yCO2_k.*exp(Hd./(R.*T_k)))./(R.*T_k) + 1).^2);
					dqCO2_dpCO2_k = (1000000.*b0.*n1.*yCO2_k.*exp(Hb./(R.*T_k)))./(R.*T_k.*((1000000.*b0.*(yCO2_k.*p_k).*yCO2_k.*exp(Hb./(R.*T_k)))./(R.*T_k) + 1)) + (1000000.*d0.*n2.*yCO2_k.*exp(Hd./(R.*T_k)))./(R.*T_k.*((1000000.*d0.*(yCO2_k.*p_k).*yCO2_k.*exp(Hd./(R.*T_k)))./(R.*T_k) + 1)) - (1000000000000.*b0.^2.*n1.*(yCO2_k.*p_k).*yCO2_k.^2.*exp((2.*Hb)./(R.*T_k)))./(R.^2.*T_k.^2.*((1000000.*b0.*(yCO2_k.*p_k).*yCO2_k.*exp(Hb./(R.*T_k)))./(R.*T_k) + 1).^2) - (1000000000000.*d0.^2.*n2.*(yCO2_k.*p_k).*yCO2_k.^2.*exp((2.*Hd)./(R.*T_k)))./(R.^2.*T_k.^2.*((1000000.*d0.*(yCO2_k.*p_k).*yCO2_k.*exp(Hd./(R.*T_k)))./(R.*T_k) + 1).^2);			
				case 'DSL2'
					dqCO2_dTequ_k = (Hb.*b0.^2.*n1.*(p_k.*yCO2_k).^2.*exp((2.*Hb)./(R.*T_k)))./(R.*T_k.^2.*(b0.*(p_k.*yCO2_k).*exp(Hb./(R.*T_k)) + 1).^2) - (Hd.*d0.*n2.*(p_k.*yCO2_k).*exp(Hd./(R.*T_k)))./(R.*T_k.^2.*(d0.*(p_k.*yCO2_k).*exp(Hd./(R.*T_k)) + 1)) - (Hb.*b0.*n1.*(p_k.*yCO2_k).*exp(Hb./(R.*T_k)))./(R.*T_k.^2.*(b0.*(p_k.*yCO2_k).*exp(Hb./(R.*T_k)) + 1)) + (Hd.*d0.^2.*n2.*(p_k.*yCO2_k).^2.*exp((2.*Hd)./(R.*T_k)))./(R.*T_k.^2.*(d0.*(p_k.*yCO2_k).*exp(Hd./(R.*T_k)) + 1).^2);
					dqCO2_dpCO2_k = (b0.*n1.*exp(Hb./(R.*T_k)))./(b0.*(p_k.*yCO2_k).*exp(Hb./(R.*T_k)) + 1) + (d0.*n2.*exp(Hd./(R.*T_k)))./(d0.*(p_k.*yCO2_k).*exp(Hd./(R.*T_k)) + 1) - (b0.^2.*n1.*(p_k.*yCO2_k).*exp((2.*Hb)./(R.*T_k)))./(b0.*(p_k.*yCO2_k).*exp(Hb./(R.*T_k)) + 1).^2 - (d0.^2.*n2.*(p_k.*yCO2_k).*exp((2.*Hd)./(R.*T_k)))./(d0.*(p_k.*yCO2_k).*exp(Hd./(R.*T_k)) + 1).^2;
				case 'toth'
                    dqCO2_dTequ_k = - b0.*ns0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi.*(T_k./T0 - 1)).*(((T0.*alpha.*log(b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).*(b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(t0 - alpha.*(T0./T_k - 1)))./T_k.^2 - (b0.*dH.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0)).*(t0 - alpha.*(T0./T_k - 1)).*(b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(t0 - alpha.*(T0./T_k - 1) - 1))./(R.*T_k.^2))./((t0 - alpha.*(T0./T_k - 1)).*((b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(t0 - alpha.*(T0./T_k - 1)) + 1).^(1./(t0 - alpha.*(T0./T_k - 1)) + 1)) - (T0.*alpha.*log((b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(t0 - alpha.*(T0./T_k - 1)) + 1))./(T_k.^2.*(t0 - alpha.*(T0./T_k - 1)).^2.*((b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(t0 - alpha.*(T0./T_k - 1)) + 1).^(1./(t0 - alpha.*(T0./T_k - 1))))) - (Xi.*b0.*ns0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi.*(T_k./T0 - 1)))./(T0.*((b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(t0 - alpha.*(T0./T_k - 1)) + 1).^(1./(t0 - alpha.*(T0./T_k - 1)))) - (b0.*dH.*ns0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi.*(T_k./T0 - 1)))./(R.*T_k.^2.*((b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(t0 - alpha.*(T0./T_k - 1)) + 1).^(1./(t0 - alpha.*(T0./T_k - 1))));
                    dqCO2_dpCO2_k = (b0.*ns0.*exp((dH.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi.*(T_k./T0 - 1)))./((b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(t0 - alpha.*(T0./T_k - 1)) + 1).^(1./(t0 - alpha.*(T0./T_k - 1))) - (b0.^2.*ns0.*(p_k.*yCO2_k).*exp((2.*dH.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi.*(T_k./T0 - 1)).*(b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(t0 - alpha.*(T0./T_k - 1) - 1))./((b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(t0 - alpha.*(T0./T_k - 1)) + 1).^(1./(t0 - alpha.*(T0./T_k - 1)) + 1);
                case 'langfr'
					dqCO2_dTequ_k = (ns0.*exp(-Xi.*(T_k./T0 - 1)).*((T0.*alpha.*log(b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).*(b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(1./t0 - alpha.*(T0./T_k - 1)))./T_k.^2 + (b0.*dH.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0)).*(alpha.*(T0./T_k - 1) - 1./t0).*(b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(1./t0 - alpha.*(T0./T_k - 1) - 1))./(R.*T_k.^2)))./((b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(1./t0 - alpha.*(T0./T_k - 1)) + 1) - (ns0.*exp(-Xi.*(T_k./T0 - 1)).*((T0.*alpha.*log(b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).*(b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(1./t0 - alpha.*(T0./T_k - 1)))./T_k.^2 + (b0.*dH.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0)).*(alpha.*(T0./T_k - 1) - 1./t0).*(b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(1./t0 - alpha.*(T0./T_k - 1) - 1))./(R.*T_k.^2)).*(b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(1./t0 - alpha.*(T0./T_k - 1)))./((b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(1./t0 - alpha.*(T0./T_k - 1)) + 1).^2 - (Xi.*ns0.*exp(-Xi.*(T_k./T0 - 1)).*(b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(1./t0 - alpha.*(T0./T_k - 1)))./(T0.*((b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(1./t0 - alpha.*(T0./T_k - 1)) + 1));
					dqCO2_dpCO2_k = (b0.*ns0.*exp((dH.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi.*(T_k./T0 - 1)).*(alpha.*(T0./T_k - 1) - 1./t0).*(b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(1./t0 - alpha.*(T0./T_k - 1)).*(b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(1./t0 - alpha.*(T0./T_k - 1) - 1))./((b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(1./t0 - alpha.*(T0./T_k - 1)) + 1).^2 - (b0.*ns0.*exp((dH.*(T0./T_k - 1))./(R.*T0)).*exp(-Xi.*(T_k./T0 - 1)).*(alpha.*(T0./T_k - 1) - 1./t0).*(b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(1./t0 - alpha.*(T0./T_k - 1) - 1))./((b0.*(p_k.*yCO2_k).*exp((dH.*(T0./T_k - 1))./(R.*T0))).^(1./t0 - alpha.*(T0./T_k - 1)) + 1);		
			end      
            % dH
            dHCO2 = R*T_k^2/(p_k*yCO2_k) *(-dqCO2_dTequ_k)/(dqCO2_dpCO2_k);
            % CO2 desorbed
            dNCO2_solid = NCO2_solid - NCO2_solid_k;

        % H2O
		if dataModel.feed.yH2O > 1e-10
            % dq/dT, inital condition
            dqH2O_dT_k = (CG0*Cm0*K0*(p_k*yH2O_k)*exp(HC/(R*T_k))*exp(HK/(R*T_k))*exp(-(80887*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k))*exp(beta/T_k)*((K0*(p_k*yH2O_k)*exp(HK/(R*T_k))*exp(-(80887*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k))*((80887*((3114378624436657125*(1 - (125*T_k)/80887)^(1/2))/728565326118234619904 - (2488235212353120375*((125*T_k)/80887 - 1)^2)/45535332882389663744 + (1123216880327743625*((125*T_k)/80887 - 1)^3)/11383833220597415936 + (2793026719395026625*(1 - (125*T_k)/80887)^(5/2))/22767666441194831872 + (7604996558327263125*(1 - (125*T_k)/80887)^(13/2))/364282663059117309952 - 553064399539058875/45535332882389663744))/(125*T_k) + (80887*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k^2)))/22064 - (HK*K0*(p_k*yH2O_k)*exp(HK/(R*T_k))*exp(-(80887*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k)))/(22064*R*T_k^2)))/(22064*((K0*(p_k*yH2O_k)*exp(HK/(R*T_k))*exp(-(80887*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k))*(CG0*exp(HC/(R*T_k)) - 1))/22064 + 1)*((K0*(p_k*yH2O_k)*exp(-((553064399539058875*T_k)/562949953421312 + (536709017657618260727*((125*T_k)/80887 - 1)^3)/562949953421312 - (726829150392561588763*((125*T_k)/80887 - 1)^4)/562949953421312 + (671767316786154359653*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (258193774001949164133*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (328077523527155910609*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 357885760684126841777/562949953421312)/(125*T_k))*exp(HK/(R*T_k)))/22064 - 1)^2) - (CG0*Cm0*K0*(p_k*yH2O_k)*exp(HC/(R*T_k))*exp(HK/(R*T_k))*exp(-(80887*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k))*exp(beta/T_k)*((80887*((3114378624436657125*(1 - (125*T_k)/80887)^(1/2))/728565326118234619904 - (2488235212353120375*((125*T_k)/80887 - 1)^2)/45535332882389663744 + (1123216880327743625*((125*T_k)/80887 - 1)^3)/11383833220597415936 + (2793026719395026625*(1 - (125*T_k)/80887)^(5/2))/22767666441194831872 + (7604996558327263125*(1 - (125*T_k)/80887)^(13/2))/364282663059117309952 - 553064399539058875/45535332882389663744))/(125*T_k) + (80887*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k^2)))/(22064*((K0*(p_k*yH2O_k)*exp(HK/(R*T_k))*exp(-(80887*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k))*(CG0*exp(HC/(R*T_k)) - 1))/22064 + 1)*((K0*(p_k*yH2O_k)*exp(HK/(R*T_k))*exp(-(80887*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k)))/22064 - 1)) - (CG0*Cm0*K0*(p_k*yH2O_k)*exp(HC/(R*T_k))*exp(HK/(R*T_k))*exp(-(80887*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k))*exp(beta/T_k)*((HK*K0*(p_k*yH2O_k)*exp(HK/(R*T_k))*exp(-(80887*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k))*(CG0*exp(HC/(R*T_k)) - 1))/(22064*R*T_k^2) - (K0*(p_k*yH2O_k)*exp(HK/(R*T_k))*exp(-(80887*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k))*(CG0*exp(HC/(R*T_k)) - 1)*((80887*((3114378624436657125*(1 - (125*T_k)/80887)^(1/2))/728565326118234619904 - (2488235212353120375*((125*T_k)/80887 - 1)^2)/45535332882389663744 + (1123216880327743625*((125*T_k)/80887 - 1)^3)/11383833220597415936 + (2793026719395026625*(1 - (125*T_k)/80887)^(5/2))/22767666441194831872 + (7604996558327263125*(1 - (125*T_k)/80887)^(13/2))/364282663059117309952 - 553064399539058875/45535332882389663744))/(125*T_k) + (80887*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k^2)))/22064 + (CG0*HC*K0*(p_k*yH2O_k)*exp(HC/(R*T_k))*exp(HK/(R*T_k))*exp(-(80887*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k)))/(22064*R*T_k^2)))/(22064*((K0*(p_k*yH2O_k)*exp(-((553064399539058875*T_k)/562949953421312 + (536709017657618260727*((125*T_k)/80887 - 1)^3)/562949953421312 - (726829150392561588763*((125*T_k)/80887 - 1)^4)/562949953421312 + (671767316786154359653*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (258193774001949164133*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (328077523527155910609*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 357885760684126841777/562949953421312)/(125*T_k))*exp(HK/(R*T_k))*(CG0*exp(HC/(R*T_k)) - 1))/22064 + 1)^2*((K0*(p_k*yH2O_k)*exp(HK/(R*T_k))*exp(-(80887*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k)))/22064 - 1)) + (CG0*Cm0*K0*beta*(p_k*yH2O_k)*exp(HC/(R*T_k))*exp(HK/(R*T_k))*exp(-(80887*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k))*exp(beta/T_k))/(22064*T_k^2*((K0*(p_k*yH2O_k)*exp(HK/(R*T_k))*exp(-(80887*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k))*(CG0*exp(HC/(R*T_k)) - 1))/22064 + 1)*((K0*(p_k*yH2O_k)*exp(HK/(R*T_k))*exp(-(80887*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k)))/22064 - 1)) + (CG0*Cm0*HC*K0*(p_k*yH2O_k)*exp(HC/(R*T_k))*exp(HK/(R*T_k))*exp(-(80887*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k))*exp(beta/T_k))/(22064*R*T_k^2*((K0*(p_k*yH2O_k)*exp(HK/(R*T_k))*exp(-(80887*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k))*(CG0*exp(HC/(R*T_k)) - 1))/22064 + 1)*((K0*(p_k*yH2O_k)*exp(HK/(R*T_k))*exp(-(80887*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k)))/22064 - 1)) + (CG0*Cm0*HK*K0*(p_k*yH2O_k)*exp(HC/(R*T_k))*exp(HK/(R*T_k))*exp(-(80887*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k))*exp(beta/T_k))/(22064*R*T_k^2*((K0*(p_k*yH2O_k)*exp(HK/(R*T_k))*exp(-(80887*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k))*(CG0*exp(HC/(R*T_k)) - 1))/22064 + 1)*((K0*(p_k*yH2O_k)*exp(HK/(R*T_k))*exp(-(80887*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k)))/22064 - 1));
            % dq/dp, inital condition
            dqH2O_dpH2O_k = (CG0*Cm0*K0^2*(p_k*yH2O_k)*exp(HC/(R*T_k))*exp((2*HK)/(R*T_k))*exp(-(161774*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k))*exp(beta/T_k))/(486820096*((K0*(p_k*yH2O_k)*exp(HK/(R*T_k))*exp(-(80887*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k))*(CG0*exp(HC/(R*T_k)) - 1))/22064 + 1)*((K0*(p_k*yH2O_k)*exp(-((553064399539058875*T_k)/562949953421312 + (536709017657618260727*((125*T_k)/80887 - 1)^3)/562949953421312 - (726829150392561588763*((125*T_k)/80887 - 1)^4)/562949953421312 + (671767316786154359653*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (258193774001949164133*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (328077523527155910609*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 357885760684126841777/562949953421312)/(125*T_k))*exp(HK/(R*T_k)))/22064 - 1)^2) - (CG0*Cm0*K0*exp(HC/(R*T_k))*exp(HK/(R*T_k))*exp(-(80887*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k))*exp(beta/T_k))/(22064*((K0*(p_k*yH2O_k)*exp(HK/(R*T_k))*exp(-(80887*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k))*(CG0*exp(HC/(R*T_k)) - 1))/22064 + 1)*((K0*(p_k*yH2O_k)*exp(HK/(R*T_k))*exp(-(80887*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k)))/22064 - 1)) + (CG0*Cm0*K0^2*(p_k*yH2O_k)*exp(HC/(R*T_k))*exp((2*HK)/(R*T_k))*exp(-(161774*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k))*exp(beta/T_k)*(CG0*exp(HC/(R*T_k)) - 1))/(486820096*((K0*(p_k*yH2O_k)*exp(-((553064399539058875*T_k)/562949953421312 + (536709017657618260727*((125*T_k)/80887 - 1)^3)/562949953421312 - (726829150392561588763*((125*T_k)/80887 - 1)^4)/562949953421312 + (671767316786154359653*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (258193774001949164133*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (328077523527155910609*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 357885760684126841777/562949953421312)/(125*T_k))*exp(HK/(R*T_k))*(CG0*exp(HC/(R*T_k)) - 1))/22064 + 1)^2*((K0*(p_k*yH2O_k)*exp(HK/(R*T_k))*exp(-(80887*((553064399539058875*T_k)/45535332882389663744 + (6635293899608321*((125*T_k)/80887 - 1)^3)/562949953421312 - (8985735042621949*((125*T_k)/80887 - 1)^4)/562949953421312 + (8305009665164419*(1 - (125*T_k)/80887)^(3/2))/4503599627370496 + (3192030536451459*(1 - (125*T_k)/80887)^(7/2))/140737488355328 + (4055998164441207*(1 - (125*T_k)/80887)^(15/2))/2251799813685248 - 4424515196312471/562949953421312))/(125*T_k)))/22064 - 1));
            % dH
            dHH2O = R*T_k^2/(p_k*yH2O_k) *(-dqH2O_dT_k)/(dqH2O_dpH2O_k);
            % H2O desorbed
            dNH2O_solid = NH2O_solid - NH2O_solid_k;
		else
			dHH2O = 0;
			dNH2O_solid = 0;
		end

        % overall  
        Q = (adsorbentMass*cp*(T-T_k) - dHCO2*dNCO2_solid - dHH2O*dNH2O_solid)/1000; %kJ
        
        out = [Q, NCO2, NH2O, NN2, qCO2, qH2O];
    end
	
	function fvec = fcn_dry(x)
	% x(1): yCO2; x(2): yN2; x(3): Nin
    
    %% equation 1: material balance CO2
    % CO2
        % solid phase CO2, initial condition
            switch CO2IsothermModel
                case 'toth_cp'
                    qCO2_k = ((ns0_c*exp(Xi_c*(1-T_k/T0)).*(b0_c*exp(dH_c/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k))./((1+((b0_c*exp(dH_c/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k)).^(t0_c+alpha_c*(1-T0./T_k))).^(1./(t0_c+alpha_c*(1-T0./T_k)))))+((ns0_p*exp(Xi_p*(1-T_k/T0)).*(b0_p*exp(dH_p/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k))./((1+((b0_p*exp(dH_p/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k)).^(t0_p+alpha_p*(1-T0./T_k))).^(1./(t0_p+alpha_p*(1-T0./T_k))))); % adsorbed amount of CO2
                case 's_shaped'
                    qCO2_k = (q_L0.*(b_L0.*exp(dU_L./(R*T_k))).*(yCO2_k*p_k)./(1+(b_L0.*exp(dU_L./(R*T_k))).*(yCO2_k*p_k))).*(1-((exp((log((yCO2_k*p_k))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k))))./(1+exp(((log((yCO2_k*p_k)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*T_k))).*(yCO2_k*p_k)./(1+(b_U0.*exp(dU_U./(R*T_k))).*(yCO2_k*p_k))+(b_H0.*exp(dU_H./(R*T_k))).*(yCO2_k*p_k)).*((exp((log((yCO2_k*p_k))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k))))./(1+exp(((log((yCO2_k*p_k)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T_k))))./(xi_1.*exp(xi_2.*(1./T0-1./T_k)))))).^gam); 
                case 'DSL'
                    qCO2_k = n1*b0*exp(Hb/R/T_k)*(yCO2_k*p_k)./(1e-6*R*T)./(1+b0*exp(Hb/R/T_k).*(yCO2_k*p_k)./(1e-6*R*T)) + n2*d0*exp(Hd/R/T_k)*(yCO2_k*p_k)./(1e-6*R*T)./(1+d0*exp(Hd/R/T_k).*(yCO2_k*p_k)./(1e-6*R*T));  
				case 'DSL2'
                    qCO2_k = n1*b0*exp(Hb/R/T_k)*(yCO2_k*p_k)./(1+b0*exp(Hb/R/T_k).*(yCO2_k*p_k)) + n2*d0*exp(Hd/R/T_k)*(yCO2_k*p_k)./(1+d0*exp(Hd/R/T_k).*(yCO2_k*p_k));            
                case 'toth'
                    qCO2_k = ((ns0*exp(Xi*(1-T_k/T0)).*(b0*exp(dH/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k))./((1+((b0*exp(dH/(R*T0)*(T0./T_k-1))).*(p_k*yCO2_k)).^(t0+alpha*(1-T0./T_k))).^(1./(t0+alpha*(1-T0./T_k)))));            
				case 'langfr'
					qCO2_k = (ns0.*exp(Xi.*(1-T_k./T0))).*((b0*exp(dH./(R*T0).*(T0./T_k-1))).*(p_k*yCO2_k)).^(1./t0+alpha.*(1-T0./T_k))./(1+ ((b0*exp(dH./(R*T0).*(T0./T_k-1))).*(p_k*yCO2_k)).^(1./t0+alpha.*(1-T0./T_k)));
			end
        NCO2_solid_k = adsorbentMass*qCO2_k;
        % fluid phase CO2, initial condition
        NCO2_fluid_k = yCO2_k*p_k*V*void/(R*T_k)*1e6; % mol
        % solid and liquid CO2, initial condition    
        NCO2_k = NCO2_solid_k + NCO2_fluid_k;

        % solid phase CO2, final condition 
            switch CO2IsothermModel
                case 'toth_cp'
                    qCO2 = ((ns0_c*exp(Xi_c*(1-T/T0)).*(b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p*x(1)))./((1+((b0_c*exp(dH_c/(R*T0)*(T0./T-1))).*(p*x(1))).^(t0_c+alpha_c*(1-T0./T))).^(1./(t0_c+alpha_c*(1-T0./T)))))+((ns0_p*exp(Xi_p*(1-T/T0)).*(b0_p*exp(dH_p/(R*T0)*(T0./T-1))).*(p*x(1)))./((1+((b0_p*exp(dH_p/(R*T0)*(T0./T-1))).*(p*x(1))).^(t0_p+alpha_p*(1-T0./T))).^(1./(t0_p+alpha_p*(1-T0./T))))); % adsorbed amount of CO2
                case 's_shaped'         
                    qCO2 = (q_L0.*(b_L0.*exp(dU_L./(R*T))).*(x(1)*p)./(1+(b_L0.*exp(dU_L./(R*T))).*(x(1)*p))).*(1-((exp((log((x(1)*p))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((x(1)*p)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam))+(q_U0.*(b_U0.*exp(dU_U./(R*T))).*(x(1)*p)./(1+(b_U0.*exp(dU_U./(R*T))).*(x(1)*p))+(b_H0.*exp(dU_H./(R*T))).*(x(1)*p)).*((exp((log((x(1)*p))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T))))./(1+exp(((log((x(1)*p)))-log(p_step0.*exp(-dH_step./R.*(1/T0-1/T))))./(xi_1.*exp(xi_2.*(1./T0-1./T)))))).^gam); % adsorbed amount of CO2                    
                case 'DSL'
                    qCO2 = n1*b0*exp(Hb/R/T)*(x(1)*p)./(1e-6*R*T)./(1+b0*exp(Hb/R/T).*(x(1)*p)./(1e-6*R*T)) + n2*d0*exp(Hd/R/T)*(x(1)*p)./(1e-6*R*T)./(1+d0*exp(Hd/R/T).*(x(1)*p)./(1e-6*R*T)); 
				case 'DSL2'
                    qCO2 = n1*b0*exp(Hb/R/T)*(x(1)*p)./(1+b0*exp(Hb/R/T).*(x(1)*p)) + n2*d0*exp(Hd/R/T)*(x(1)*p)./(1+d0*exp(Hd/R/T).*(x(1)*p));            
                case 'toth'
                    qCO2 = ((ns0*exp(Xi*(1-T/T0)).*(b0*exp(dH/(R*T0)*(T0./T-1))).*(x(1)*p))./((1+((b0*exp(dH/(R*T0)*(T0./T-1))).*(x(1)*p)).^(t0+alpha*(1-T0./T))).^(1./(t0+alpha*(1-T0./T)))));			
				case 'langfr'
					qCO2 = (ns0.*exp(Xi.*(1-T./T0))).*((b0*exp(dH./(R*T0).*(T0./T-1))).*(x(1)*p)).^(1./t0+alpha.*(1-T0./T))./(1+ ((b0*exp(dH./(R*T0).*(T0./T-1))).*(x(1)*p)).^(1./t0+alpha.*(1-T0./T)));
			end
        NCO2_solid = adsorbentMass*qCO2;
        % fluid phase CO2, final condition  
        NCO2_fluid = x(1)*p*V*void/(R*T)*1e6; % mol
        % solid and liquid CO2, final condition  
        NCO2 = NCO2_solid + NCO2_fluid;

        % material balance CO2
        fvec(1) = - NCO2 + yCO2_feed*x(3) + NCO2_k; % == 0

    %% equation 2 material balance N2
    % N2
        % fluid phase N2, initial condition
        NN2_k = yN2_k*p_k*V*void/(R*T_k)*1e6; % mol
        % fluid phase N2, final condition
        NN2 = x(2)*p*V*void/(R*T)*1e6; % mol

        % material balance N2 
		fvec(2) =  - (NN2 + NCO2) + x(3) + (NCO2_k + NN2_k); % == 0

    %% equation 4: overall material balance
        fvec(3) = 1-x(1)-x(2); % == 0
		
		fvec = real(fvec);
    end

end