

%% calculate temperature profile and heating time

function [THeatProfileStep, timeHeatingNN, TimeHeatProfilePlot] = temperatureProfile(dataModel,Tdes)

    % 1. calculate time dependent on the specific sorbent (rho) and time
    % 2. make the profile the same length as the time
    
    % calculate time using NN
    load('net_new4','netNN_new');
    ParamNN = netNN_new([dataModel.sorbent.MaterialDensity;...
        dataModel.sorbent.Density;...
        dataModel.sorbent.cp; ...
        dataModel.process.Tdes;dataModel.process.pvac]);
    timeHeatingNN = ParamNN; % 

    % for the calculation, reduce/increase no of steps to 100
    max_steps = dataModel.process.noSteps;
    
    % starting point
        THeatProfileStep = ones(1,max_steps);   
        TimeHeatProfilePlot = ones(1,max_steps);  
    
    % profile
        p = [-687.316091760085,661.535735641533,42.0052174647580,-2791.97322593205,-4.29623797695230,-0.345664545572582,-34162.5562758762,0.00867252521860531]; % from fitting
        a = (timeHeatingNN-1)/(max_steps-1); % 
        b = (2606-1)/(max_steps-1);
        timeEnd = b;
        THeatProfileStep(1) = dataModel.process.Tamb;
        TimeHeatProfilePlot(1) = 0;
            for m = (2:(max_steps-1))
                THeatProfileStep(m) = (p(1)+Tdes-p(4).*Tdes.^p(6)).*atan((timeEnd-p(2).*1e9.*Tdes.^p(5))./p(3))+(Tdes-p(7)).*p(8);
                TimeHeatProfilePlot(m) = m*a;
                timeEnd = timeEnd+b;
            end
        
    % end point
        THeatProfileStep(max_steps) =  dataModel.process.Tdes;
        TimeHeatProfilePlot(max_steps) = timeHeatingNN;
end