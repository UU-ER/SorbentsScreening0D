
function [pressureVector, timeBD, TimeBDProfilePlot] = pressureProfile(dataModel,pvac)

    pvac = pvac*10; % bar
    x = 1;
    y = 0;
    deltaP = abs(pvac-y);
    pp = [0.186277700879350,-0.123808829948893,-7.13805471361795e-05];
    while deltaP>1e-3 && x<10000
        y = (pp(1))./(pvac.^0.85).*exp((pp(2)).*x) +...
            pvac.*exp((pp(3))./(pvac.^0.85).*x);
        x = x+1;
        deltaP = abs(pvac-y);
    end

    timeBD = floor((-10*pvac+20.7)./100*200);

    max_steps = dataModel.process.noSteps;
    timeEnd1 = 0;
    a = (x-1)/(max_steps-2);
    b = (timeBD-1)/(max_steps-1); % for time profile
    timeEnd2 = b;
    pressureVector(1) = dataModel.process.pamb;
    TimeBDProfilePlot(1) = 0;
    for m = (2:(max_steps-1))
        pressureVector(m) = ((pp(1))./(pvac.^0.85).*exp((pp(2)).*timeEnd1) +...
            pvac.*exp((pp(3))./(pvac.^0.85).*timeEnd1))./10; % MPa
        timeEnd1 = timeEnd1+a;
        TimeBDProfilePlot(m) = m*b;
        timeEnd2 = timeEnd2+b;
    end
    pressureVector(max_steps) = dataModel.process.pvac;
    TimeBDProfilePlot(max_steps) = timeBD;

end