function [KH_save, ns_equ_save, q0_save, m0_save, gamma_save, wc_save, dH_save,T0] = characteristics_4(p, model,pCO2_1,rhoMat,T0,Tdes)
%========================================================
% Calculating the different characterisitics according to the isotherm
% model used.
%-------------------------------------------------------
%Input:     - p Parameters belonging to the model
%           - model Model used to fit the data, characteristics will be
%               computed accordingly.
%           - pCO2_1 Partial pressure of adsorption of CO2
% Output:  - KH_save Value Henry's coeffecient
%          - ns_equ_save Value equilibrium loading
%          - q0_save Value amount adsorbed at feed concentration
%          - m0_save Value slope isotherm at feed concentration
%          - gamma_save Non-linearity
%          - wc_save Working capacity
%          - dH_save Heat of adsorption
%========================================================

%% input from fitting
    no = 1;
    pCO2_2 = 1e4; % partial pressure desorption, Pa
    
    T1 = 298;
    T2 = Tdes+273;

    % different isotherm model
    switch model
        case 'langm_freund'
            % parameters langmuir freundlich
            % 1: n0 [mmol/g], 2: Xi [-], 3: t0 [-], 4: alpha [-], 5:b0 [1/Pa], 6:Q [kJ/mol]
            if isempty(T0)
                T0 = 293;
            end
        case 'toth_cp'
            % Xi_c, dH_c, alpha_c, Xi_p, dH_p [kJ/mol] ,alpha_p, ns0_c, b0_c, t0_c, ns0_p, b0_p [1/Pa],    t0_p
            if isempty(T0)
                T0 = 293;
            end
            % Xi_c, dH_c, alpha_c, Xi_p, dH_p, alpha_p, ns0_c, b0_c, t0_c, ns0_p, b0_p, t0_p: [-], [kJ/mol], [-], [-], [kJ/mol], [-]
        case 's_shaped'
            if isempty(T0)
                T0 = 313.5;
            end
            % q_L0(mol/kg), b_L0(1/Pa), dU_L(kJ/mol), q_U0(kJ/mol),b_U0(1/Pa), dU_U(kJ/mol), b_H0(mol/kg/Pa), dU_H(kJ/mol), Xi_1(-), Xi_2(-), pstep_0(Pa), dH_step(kJ/mol), gamma(-)
        case 'DSL'
            T0 = 0; 
            % n1 (mol/m3), b0 (m3/mol), Hb (J/mol), n2 (mol/m3), d0 (m3/mol), Hd (J/mol)
        case 'toth'
            if isempty(T0)
                T0 = 293;
            end
            % Xi, dH, alpha, ns0, b0, t0[-], [kJ/mol], [-]
    end

%% calculate characteristics
switch model
    case 'langm_freund'
        % 1: n0 [mmol/g], 2: Xi [-], 3: t0 [-], 4: alpha [-], 5:b0 [1/Pa], 6:Q [kJ/mol]
            R = 8.314/1000; % (kJ/mol/K)
        % Slope at low pressure
            T = T1;
            xm = 0.05;
            b = p(5)*exp(p(6)/(R*T0)*(T0/T-1));
            ns = p(1)*exp(p(2)*(1-T/T0));
            t = (1/p(3)+p(4)*(1-T0/T));
            q0 = ns*(b*xm)^(t)/(1+ (b*xm)^(t));
            K_H = q0/xm;
            KH_save(:,no) = K_H;
        % equilibrium loading
            ns_equ = ns;
            ns_equ_save(:,no) = ns_equ;
        % CO2 loading at feed concentration
            xm = pCO2_1; % Pa
            t = (1/p(3)+p(4)*(1-T0/T)); % = 1/t
            q0 = ns*(b*xm)^(t)/(1+ (b*xm)^(t));
            q0_save(:,no) = q0;
        % Local slope at feed concentration dq_0/dxm
        xm = pCO2_1; % Pa
        T = T1;
            n0 = p(1);
            Xi = p(2);
            t0 = p(3);
            alpha = p(4);
            b0 = p(5);
            Q = p(6);
            m0 =       (b0*n0*exp((Q*(T0/T - 1))/(R*T0))*exp(-Xi*(T/T0 - 1))*(alpha*(T0/T - 1) - 1/t0)*(b0*xm*exp((Q*(T0/T - 1))/(R*T0)))^(1/t0 - alpha*(T0/T - 1))*(b0*xm*exp((Q*(T0/T - 1))/(R*T0)))^(1/t0 - alpha*(T0/T - 1) - 1))/((b0*xm*exp((Q*(T0/T - 1))/(R*T0)))^(1/t0 - alpha*(T0/T - 1)) + 1)^2 - (b0*n0*exp((Q*(T0/T - 1))/(R*T0))*exp(-Xi*(T/T0 - 1))*(alpha*(T0/T - 1) - 1/t0)*(b0*xm*exp((Q*(T0/T - 1))/(R*T0)))^(1/t0 - alpha*(T0/T - 1) - 1))/((b0*xm*exp((Q*(T0/T - 1))/(R*T0)))^(1/t0 - alpha*(T0/T - 1)) + 1);
            m0_save(:,no) = m0;
        % nonlinearity
            k=0;
            for r = 1:40 
                k=k+1;
                xm=r;
                m(k) = (b0*n0*exp((Q*(T0/T - 1))/(R*T0))*exp(-Xi*(T/T0 - 1))*(alpha*(T0/T - 1) - 1/t0)*(b0*xm*exp((Q*(T0/T - 1))/(R*T0)))^(1/t0 - alpha*(T0/T - 1))*(b0*xm*exp((Q*(T0/T - 1))/(R*T0)))^(1/t0 - alpha*(T0/T - 1) - 1))/((b0*xm*exp((Q*(T0/T - 1))/(R*T0)))^(1/t0 - alpha*(T0/T - 1)) + 1)^2 - (b0*n0*exp((Q*(T0/T - 1))/(R*T0))*exp(-Xi*(T/T0 - 1))*(alpha*(T0/T - 1) - 1/t0)*(b0*xm*exp((Q*(T0/T - 1))/(R*T0)))^(1/t0 - alpha*(T0/T - 1) - 1))/((b0*xm*exp((Q*(T0/T - 1))/(R*T0)))^(1/t0 - alpha*(T0/T - 1)) + 1);
            end
            m_rmse = sqrt(sum((m-mean(m)).^2)/size(m,2));
            gamma_save(:,no) = m_rmse; 
        % working capacity CO2
            T = T2;
            xm = pCO2_2; % Pa
            b = p(5)*exp(p(6)/(R*T0)*(T0/T-1));
            ns = p(1)*exp(p(2)*(1-T/T0));
            t = (1/p(3)+p(4)*(1-T0/T)); % = 1/t
            q1 = ns*(b*xm)^(t)/(1+ (b*xm)^(t));
            wc = q0-q1;
            wc_save(:,no) = wc; 
        % heat of adsorption
            dH = p(6); % kJ/mol
            dH_save(:,no) = dH;
    case 'toth_cp'
            R = 8.314/1000; % (kJ/mol/K)
            % Xi_c, dH_c (kJ/mol), alpha_c, Xi_p, dH_p (kJ/mol), alpha_p, ns0_c (mol/kg), b0_c (1/Pa), t0_c, ns0_p (mol/kg), b0_p (1/Pa), t0_p
        % Slope at low pressure
            xm = 0.05;
            q0 = ((p(7)*exp(p(1)*(1-T/T0)).*(p(8)*exp(p(2)/(R*T0)*(T0./T-1))).*xm)./((1+((p(8)*exp(p(2)/(R*T0)*(T0./T-1))).*xm).^(p(9)+p(3)*(1-T0./T))).^(1./(p(9)+p(3)*(1-T0./T)))))+((p(10)*exp(p(4)*(1-T/T0)).*(p(11)*exp(p(5)/(R*T0)*(T0./T-1))).*xm)./((1+((p(11)*exp(p(5)/(R*T0)*(T0./T-1))).*xm).^(p(12)+p(6)*(1-T0./T))).^(1./(p(12)+p(6)*(1-T0./T))))); % adsorbed amount of CO2
            K_H = q0/xm; % 
            KH_save(:,no) = K_H;
        % equilibrium loading
            ns_equ = nsc+nsp;
            ns_equ_save(:,no) = ns_equ;
        % CO2 loading at feed concentration
            T = T1; % K
            xm = pCO2_1; % Pa
            q0 = ((p(7)*exp(p(1)*(1-T/T0)).*(p(8)*exp(p(2)/(R*T0)*(T0./T-1))).*xm)./((1+((p(8)*exp(p(2)/(R*T0)*(T0./T-1))).*xm).^(p(9)+p(3)*(1-T0./T))).^(1./(p(9)+p(3)*(1-T0./T)))))+((p(10)*exp(p(4)*(1-T/T0)).*(p(11)*exp(p(5)/(R*T0)*(T0./T-1))).*xm)./((1+((p(11)*exp(p(5)/(R*T0)*(T0./T-1))).*xm).^(p(12)+p(6)*(1-T0./T))).^(1./(p(12)+p(6)*(1-T0./T))))); % adsorbed amount of CO2
            q0_save(:,no) = q0;
        % Local slope at feed concentration dq_0/dxm
            xm = pCO2_1; % Pa
            t_c = p(9)+p(3)*(1-T0/T1); 
            t_p = p(12)+p(6)*(1-T0/T1); 
            m0 = b_c*nsc*(1+(b_c*xm)^t_c)^(-((1+t_c)/t_c))+b_p*nsp*(1+(b_p*xm)^t_p)^(-((1+t_p)/t_p));
            m0_save(:,no) = m0;
        % nonlinearity
            T = T1;
            t_c = p(9)+p(3)*(1-T0/T); 
            t_p = p(12)+p(6)*(1-T0/T); 
            nsc = p(7)*exp(p(1)*(1-(T/T0)));
            k=0;
            for r = 1:40 
                k=k+1;
                xm=r;
                m(k) = b_c*nsc*(1+(b_c*xm)^t_c)^(-((1+t_c)/t_c))+b_p*nsp*(1+(b_p*xm)^t_p)^(-((1+t_p)/t_p));
            end
            m_rmse = sqrt(sum((m-mean(m)).^2)/size(m,2));
            gamma_save(:,no) = m_rmse;        
        % working capacity CO2
            T = T2; % K
            xm = pCO2_2; % Pa
            q1 = ((p(7)*exp(p(1)*(1-T/T0)).*(p(8)*exp(p(2)/(R*T0)*(T0./T-1))).*xm)./((1+((p(8)*exp(p(2)/(R*T0)*(T0./T-1))).*xm).^(p(9)+p(3)*(1-T0./T))).^(1./(p(9)+p(3)*(1-T0./T)))))+((p(10)*exp(p(4)*(1-T/T0)).*(p(11)*exp(p(5)/(R*T0)*(T0./T-1))).*xm)./((1+((p(11)*exp(p(5)/(R*T0)*(T0./T-1))).*xm).^(p(12)+p(6)*(1-T0./T))).^(1./(p(12)+p(6)*(1-T0./T)))));
            wc = q0-q1;
            wc_save(:,no) = wc;   
        % heat of adsorption
            dH_p = p(5);
            dH_c = p(2);
            dH = dH_c+dH_p;
            dH_save(:,no) = dH; 
        
    case 's_shaped'
        % q_L0(mol/kg), b_L0(1/Pa), dU_L(kJ/mol), q_U0(kJ/mol),b_U0(1/Pa), dU_U(kJ/mol), b_H0(mol/kg/Pa), dU_H(kJ/mol), Xi_1(-), Xi_2(-), pstep_0(Pa), dH_step(kJ/mol), gamma(-)
        R = 8.314/1000; % (kJ/mol/K)
        % Slope at low pressure
            xm = 0.05;
            T = T1; % K
            b_L = p(2)*exp(p(3)/(R*T)); % b
            b_U = p(5)*exp(p(6)/(R*T)); % d
            b_H = p(7)*exp(p(8)/(R*T)); % bb
            q_L = p(1)*b_L*xm./(1+b_L*xm); % low pressure region
            q_U = p(4)*b_U.*xm./(1+b_U.*xm)+b_H.*xm; % high pressure region ->  changed this, beacuse I want a saturation loading for my characterstics
            p_step = p(11)*exp(-p(12)/R*(1/T0-1/T)); % ps
            sigm = p(9)*exp(p(10)*(1/T0-1/T)); % s
            w = (exp((log(xm)-log(p_step))./sigm)./(1+exp((log(xm)-log(p_step))./sigm))).^p(13); % 
            q0 = q_L.*(1-w)+q_U.*w;
            KH_save(:,no) = q0/xm;
        % saturation loading -> no real value
            ns_equ = p(4);
            ns_equ_save(:,no) = ns_equ;
        % CO2 loading at feed concentration
            T = T1; % K
            xm = pCO2_1*1e-06; % MPa
                b_L = p(2)*exp(p(3)/(R*T)); % b
                b_U = p(5)*exp(p(6)/(R*T)); % d
                b_H = p(7)*exp(p(8)/(R*T)); % bb
                q_L = p(1)*b_L*xm./(1+b_L*xm); % low pressure region
                q_U = p(4)*b_U.*xm./(1+b_U.*xm)+b_H.*xm; % high pressure region ->  changed this, beacuse I want a saturation loading for my characterstics
                p_step = p(11)*exp(-p(12)/R*(1/T0-1/T)); % ps
                sigm = p(9)*exp(p(10)*(1/T0-1/T)); % s
                w = (exp((log(xm)-log(p_step))./sigm)./(1+exp((log(xm)-log(p_step))./sigm))).^p(13); % 
               q0 = q_L.*(1-w)+q_U.*w;
            q0_save(:,no) = q0;
        % Local slope at feed concentration dq_0/dxm
            xm = pCO2_1*1e-06; % MPa
                b_L = p(2)*exp(p(3)/(R*T));
                b_U = p(5)*exp(p(6)/(R*T));
                b_H = p(7)*exp(p(8)/(R*T));
                p_step = p(11)*exp(-p(12)/R*(1/T0-1/T));
                sigm = p(9)*exp(p(10)*(1/T0-1/T));
                m0 = (exp((log(xm) - log(p_step)) /sigm)./(exp((log(xm) - log(p_step)) /sigm) + 1)).^p(13)*(b_H + (b_U*p(4))./(b_U*xm + 1) -  (b_U^2*xm*p(4))./(b_U*xm + 1)^2) + p(13)*(exp((log(xm) - log(p_step))./sigm)./(xm*sigm*(exp((log(xm) - log(p_step))./sigm) + 1)) - exp((2*(log(xm) - log(p_step)))./sigm)/(xm*sigm*(exp((log(xm) - log(p_step))./sigm) + 1)^2)).*(b_H*xm + (b_U*xm*p(4))./(b_U*xm + 1)).*(exp((log(xm) - log(p_step))./sigm)./(exp((log(xm) - log(p_step))./sigm) + 1)).^(p(13) - 1) - (b_L*p(1).*((exp((log(xm) - log(p_step))./sigm)./(exp((log(xm) - log(p_step))./sigm) + 1)).^p(13) - 1))./(b_L*xm + 1) + (b_L^2.*xm*p(1).*((exp((log(xm) - log(p_step))./sigm)./(exp((log(xm) - log(p_step))./sigm) + 1)).^p(13) - 1))./(b_L*xm + 1).^2 - (b_L*p(13)*xm*p(1)*(exp((log(xm) - log(p_step))./sigm)/(xm*sigm*(exp((log(xm) - log(p_step))./sigm) + 1)) - exp((2*(log(xm) - log(p_step)))./sigm)/(xm*sigm*(exp((log(xm) - log(p_step))./sigm) + 1)^2)).*(exp((log(xm) - log(p_step))./sigm)/(exp((log(xm) - log(p_step))./sigm) + 1)).^(p(13) - 1))./(b_L*xm + 1);
            m0_save(:,no) = m0/1e5;
        % nonlinearity
            T = T1;
            k=0;
            for r = 1:40 
                k=k+1;
                xm=r;
                m(k) = (exp((log(xm) - log(p_step))./sigm)./(exp((log(xm) - log(p_step))./sigm) + 1)).^p(13)*(b_H + (b_U*p(4))./(b_U*xm + 1) - (b_U^2*xm*p(4))./(b_U*xm + 1)^2) + p(13)*(exp((log(xm) - log(p_step))./sigm)./(xm*sigm*(exp((log(xm) - log(p_step))./sigm) + 1)) - exp((2*(log(xm) - log(p_step)))./sigm)./(xm.*sigm*(exp((log(xm) - log(p_step))./sigm) + 1)^2)).*(b_H*xm + (b_U*xm*p(4))./(b_U*xm + 1)).*(exp((log(xm) - log(p_step))./sigm)./(exp((log(xm) - log(p_step))./sigm) + 1)).^(p(13) - 1) - (b_L*p(1).*((exp((log(xm) - log(p_step))./sigm)./(exp((log(xm) - log(p_step))./sigm) + 1)).^p(13) - 1))./(b_L.*xm + 1) + (b_L^2.*xm.*p(1)*((exp((log(xm) - log(p_step))./sigm)./(exp((log(xm) - log(p_step))./sigm) + 1)).^p(13) - 1))./(b_L*xm + 1)^2 - (b_L*p(13)*xm*p(1).*(exp((log(xm) - log(p_step))./sigm)./(xm*sigm*(exp((log(xm) - log(p_step))./sigm) + 1)) - exp((2*(log(xm) - log(p_step)))./sigm)./(xm.*sigm.*(exp((log(xm) - log(p_step))./sigm) + 1)^2)).*(exp((log(xm) - log(p_step))/sigm)./(exp((log(xm) - log(p_step))./sigm) + 1)).^(p(13) - 1))./(b_L*xm + 1);
            end
            m_rmse = sqrt(sum((m-mean(m)).^2)/size(m,2));
            gamma_save(:,no) = m_rmse;
        % working capacity CO2
            T = T2; % K
            xm = pCO2_2*1e-06; % MPa
                b_L = p(2)*exp(p(3)/(R*T)); % b
                b_U = p(5)*exp(p(6)/(R*T)); % d
                b_H = p(7)*exp(p(8)/(R*T)); % bb
                q_L = p(1)*b_L*xm./(1+b_L*xm); % low pressure region
                q_U = p(4)*b_U.*xm./(1+b_U.*xm)+b_H.*xm; % high pressure region ->  changed this, beacuse I want a saturation loading for my characterstics
                p_step = p(11)*exp(-p(12)/R*(1/T0-1/T)); % ps
                sigm = p(9)*exp(p(10)*(1/T0-1/T)); % s
                w = (exp((log(xm)-log(p_step))./sigm)./(1+exp((log(xm)-log(p_step))./sigm))).^p(13); % 
                q1 = q_L.*(1-w)+q_U.*w;
            wc = q0-q1;
            wc_save(:,no) = wc; 
        % heat of adsorption
            dH = p(3)+p(6)+p(8);
            dH_save(:,no) = dH;
            
    case 'DSL'  
        % qs1 (mol/kg), b0 (m3/mol), dH1 (J/mol), qs2 (mol/kg), d0 (m3/mol), dH2 (J/mol)
            R = 8.314; % (J/mol/K) = (m3*Pa/mol/K)
            % Slope at low pressure
                xm = 0.05;
                T = T1; % K
                q0 = p(1)*p(2)*exp(p(3)/R/T)*xm./(R*T)./(1+p(2)*exp(p(3)/R/T)*xm./(R*T))...
                    + p(4)*p(5)*exp(p(6)/R/T)*xm./(R*T)./(1+p(5)*exp(p(6)/R/T)*xm./(R*T));    
                K_H = q0/xm; % y=mx+c -> c=0, m=y/x
                KH_save(:,no) = K_H;
            % equilibrium loading
                ns_equ = (n1+n2); % (mol/kg)
                ns_equ_save(:,no) = ns_equ;
            % CO2 loading at feed concentration
                xm = pCO2_1; % Pa
                T = T1;
                q0 = p(1)*p(2)*exp(p(3)/R/T)*xm./(R*T)./(1+p(2)*exp(p(3)/R/T)*xm./(R*T))...
                        + p(4)*p(5)*exp(p(6)/R/T)*xm./(R*T)./(1+p(5)*exp(p(6)/R/T)*xm./(R*T));    
                q0_save(:,no) = q0;
            % Local slope at feed concentration dq_0/dxm
            xm = pCO2_1; % Pa
            T = T1;
            yCO2 = 1;
                m0 = (p(2)*p(1)*yCO2*exp(p(3)/(R*T)))./(R*T*((p(2).*xm.*yCO2*exp(p(3)./(R*T)))./(R*T) + 1))...
                    + (p(5)*p(4)*yCO2*exp(p(6)/(R*T)))./(R*T*((p(5).*xm.*yCO2*exp(p(6)/(R*T)))./(R*T) + 1))...
                    - (p(2)^2*p(1)*xm*yCO2^2*exp((2*p(3))./(R*T)))./(R^2*T^2*((p(2).*xm.*yCO2*exp(p(3)./(R*T)))./(R*T) + 1)^2)...
                    - (p(5)^2*p(4)*xm*yCO2^2*exp((2*p(6))/(R*T)))./(R^2*T^2*((p(5).*xm.*yCO2*exp(p(6)./(R*T)))./(R*T) + 1)^2);
                m0_save(:,no) = m0;
            % nonlinearity
                T = T1;
                k=0;
                for r = 1:40 
                    k=k+1;
                    xm=r;
                    m(k) = (p(2)*p(1)*yCO2*exp(p(3)/(R*T)))./(R*T*((p(2).*xm.*yCO2*exp(p(3)./(R*T)))./(R*T) + 1))...
                            + (p(5)*p(4)*yCO2*exp(p(6)/(R*T)))./(R*T*((p(5).*xm.*yCO2*exp(p(6)/(R*T)))./(R*T) + 1))...
                            - (p(2)^2*p(1)*xm*yCO2^2*exp((2*p(3))./(R*T)))./(R^2*T^2*((p(2).*xm.*yCO2*exp(p(3)./(R*T)))./(R*T) + 1)^2)...
                            - (p(5)^2*p(4)*xm*yCO2^2*exp((2*p(6))/(R*T)))./(R^2*T^2*((p(5).*xm.*yCO2*exp(p(6)./(R*T)))./(R*T) + 1)^2);
                end
                m_rmse = sqrt(sum((m-mean(m)).^2)/size(m,2));
                gamma_save(:,no) = m_rmse; 
            % working capacity CO2
                T = T2;
                xm = pCO2_2; % Pa
                q1 =  p(1)*p(2)*exp(p(3)/R/T)*xm./(R*T)./(1+p(2)*exp(p(3)/R/T)*xm./(R*T))...
                                + p(4)*p(5)*exp(p(6)/R/T)*xm./(R*T)./(1+p(5)*exp(p(6)/R/T)*xm./(R*T));  
                wc = q0-q1;
                wc_save(:,no) = wc; 
            % heat of adsorption
                dH = (p(3)/1000+p(6)/1000); % kJ/mol
                dH_save(:,no) = dH;
    case 'toth'
        % Xi, dH(kJ/mol), alpha, ns0(mol/kg), b0 [1/Pa], t0
            R = 8.314/1000; % (kJ/mol/K)
        % Henry's law coefficient
            xm = 0.05;
            T = T1; % K
            q0 = ((p(4)*exp(p(1)*(1-T/T0)).*(p(5)*exp(p(2)/(R*T0)*(T0./T-1))).*xm)./((1+((p(5)*exp(p(2)/(R*T0)*(T0./T-1))).*xm).^(p(6)+p(3)*(1-T0./T))).^(1./(p(6)+p(3)*(1-T0./T))))); % adsorbed amount of CO2
            K_H = q0/xm;
            KH_save(:,no) = K_H;
        % equilibrium loading
            ns_equ = ns;
            ns_equ_save(:,no) = ns_equ;
        % CO2 loading at feed concentration
            T = T1; % K
            xm = pCO2_1; % Pa
            q0 = ((p(4)*exp(p(1)*(1-T/T0)).*(p(5)*exp(p(2)/(R*T0)*(T0./T-1))).*xm)./((1+((p(5)*exp(p(2)/(R*T0)*(T0./T-1))).*xm).^(p(6)+p(3)*(1-T0./T))).^(1./(p(6)+p(3)*(1-T0./T))))); % adsorbed amount of CO2
            q0_save(:,no) = q0;
        % Local slope at feed concentration dq_0/dxm
            xm = pCO2_1; % Pa
            t = p(6)+p(3)*(1-T0/T1); 
            m0 = (b*ns)/((b*xm)^t + 1)^t - (b^2*ns*t^2*xm*(b*xm)^(t - 1))/((b*xm)^t + 1)^(t + 1);
            m0_save(:,no) = m0;
        % nonlinearity
            T = T1;
            k=0;
            for r = 1:40 
                k=k+1;
                xm=r;
                m(k) = (b*ns)/((b*xm)^t + 1)^t - (b^2*ns*t^2*xm*(b*xm)^(t - 1))/((b*xm)^t + 1)^(t + 1);
            end
            m_rmse = sqrt(sum((m-mean(m)).^2)/size(m,2));
            gamma_save(:,no) = m_rmse;        
        % working capacity CO2
            T = T2; % K
            xm = pCO2_2; % Pa
            q1 = ((p(4)*exp(p(1)*(1-T/T0)).*(p(5)*exp(p(2)/(R*T0)*(T0./T-1))).*xm)./((1+((p(5)*exp(p(2)/(R*T0)*(T0./T-1))).*xm).^(p(6)+p(3)*(1-T0./T))).^(1./(p(6)+p(3)*(1-T0./T)))));
            wc = q0-q1;
            wc_save(:,no) = wc;   
        % heat of adsorption
            dH = p(2);
            dH_save(:,no) = dH; 
    end
end

