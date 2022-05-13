function E= cmp_W_vac2(pin,pout,T,Nout) % kJ vacuum pump CO2 prod
poleffHH = 0.72; % efficiency vacuum pump p = 1.01 bar
poleffH = 0.7; % efficiency vacuum pump p = 0.8 bar
poleffM = 0.6; % efficiency vaccum pump p = 0.1 bar
poleffL = 0.3; % efficiency vacuum pump p = 0.01 bar
poleff=interp1([0.01 0.1 0.8 1.01],[poleffL poleffM poleffH poleffHH],pin);
R0 = 8.314; % J/mol/K

cp_in = 37/44; % J/mol/K
gamma = cp_in/(cp_in - R0/44);
theta = (gamma-1)/gamma;

E = 1/poleff*1/theta*R0*T*((pout/pin)^(theta)-1)*sum(Nout)/1000; % 

if isnan(E)
    E=0;
else 
end

end