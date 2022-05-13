function [TCoolProfileStep, TimeCoolProfilePlot, pressureVectorCool] = coolingProfile_1D(dataModel,Tdes)

    pvac = dataModel.process.pvac;
    p = [0.174853800090986,-0.0262921350934970];

    timeCooling = dataModel.process.coolingTime-1; % directly from 1D model

    % pressure vector: very fast, takes only 50 seconds 
    xm = linspace(0,50,50);
    pp = [0.000221560995722903,0.00170366982722347,0.141834668956426,0.000322209875008416,0.00168615364119241];
    pressureVector = pvac.*10.*exp(((pp(1))-(pp(4)).*pvac.*10).*xm)+((pp(2))-(pp(5)).*pvac.*10).*exp((pp(3)).*xm); % bar
    pressureVector = [pressureVector./10, ones(1,(timeCooling-30)).*dataModel.process.pamb]; % MPa

    max_steps = dataModel.process.noSteps;

    timeEnd = 0;
    timeEnd2 = 0;
    TCoolProfileStep(1) = Tdes;
    a = (timeCooling-1)/max_steps;
    b = (timeCooling-1)/max_steps; % for time profile
    pressureVectorCool(1) = pvac;
        for m = (1:max_steps)
            pressureVectorCool(m+1) = pressureVector(round(m*a));
            TCoolProfileStep(m+1) = 293 + (p(1)).*Tdes*exp((p(2)).*timeEnd); % 
            timeEnd = timeEnd+a;
            TimeCoolProfilePlot(m) = m*b;
            timeEnd2 = timeEnd2+b;
        end
end