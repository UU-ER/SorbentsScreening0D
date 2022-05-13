function [p_vap] = vaporPressure(T)
    
    % parameters:
    Tcrit = 647.096;                  % K critical temperature
    Pcrit = 22.064;                   % MPa critical pressure


    % model parameters from Wagner, Pruss 2002 equation 2.5
    a1 =- 7.85951783;
    a2 =  1.84408259;
    a3 =-11.7866497;
    a4 = 22.6807411;
    a5 =-15.9618719;
    a6 =  1.80122502;

    x  = 1-(T/647.096);
    p_vap = Pcrit * exp((Tcrit/T)*(a1*x+a2*x^(1.5)+a3*x^3+a4*x^(3.5)+a5*x^4+a6*x^(7.5))); % MPa
end

% p_vap = Pcrit * exp((Tcrit/T)*(a1*x+a2*x^(1.5)+a3*x^3+a4*x^(3.5)+a5*x^4+a6*x^(7.5))); % MPa