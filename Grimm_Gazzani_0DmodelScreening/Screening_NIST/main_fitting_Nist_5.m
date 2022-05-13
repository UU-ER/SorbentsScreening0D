 function [xm,ym,ymSpec,sigmaCO2,p_best,R2_best] = main_fitting_Nist_5 (fitting_input,temperature_name_for_input,every_gas_name, model) 
 %========================================================
% Fits isotherm data to an isotherm model. Adapted from original code by Alexa Grimm.
%-------------------------------------------------------
%Input:   - fitting_input Structure containing the data that will be fitted
%         - temperature_name_for_input Temperature that is fitted
%         - every_gas_name Name of gas 
%         - model Model used for fitting
%Output:  - xm, ym, ymSpec, sigmaCO2 Fitted data
%         - p Parameters 
%         - R2 Coefficient of determination 
%========================================================

%% Fitting isotherms with data from NIST databank
    clearvars objective p ym xm yp p0 A b Aeq beq low up nonlcon options %debug
    %% Define model, error function, optimizer

    % define how you want to calculate the error
    error = 'R2'; % R2

    % choose the optimizer
    optimizer = 'fmin'; % fmin, gs, pattern

    %% Data input experimental points
    R = 8.314/1000; % kJ/K/mol
 
     clearvars xm ym ymSpec sigmaCO2
     for every_temp = [1 : length(temperature_name_for_input)]
        temperature_name_for_input(every_temp) = string(temperature_name_for_input(every_temp)); 
        T_name = string (temperature_name_for_input (every_temp));
        
        T_name2num = strip(T_name,'T');
        T_name2num = strip(T_name2num,'a');
        T_name2num = strip(T_name2num,'A');
        T_name2num = strip(T_name2num,'a');
        
        index_num = strfind(T_name2num,"num");
        if isempty(index_num) == 0
            if index_num>0
                T_name2num = replaceBetween(T_name2num,index_num, length(char(T_name2num)),'');
            end
        end
        
        T(every_temp) = str2num(T_name2num);
       
        if every_temp == 1
            exp_data1 = transpose(fitting_input.(every_gas_name).(T_name)); % column 1 is amount adsorbed, columnn 2 is pressure!
            first_row = exp_data1 (1,:);
            while isempty(exp_data1) ==0
                while (first_row(1) == 0 |first_row(2) == 0) && isempty(exp_data1) ==0
                    exp_data1(1,:) = [];
                    if size(exp_data1,1)>0
                        first_row = exp_data1 (1,:);
                    end
                end
                if first_row(1) ~= 0 && first_row(2) ~= 0
                    break
                end
            end
            [row1, column1] = size(exp_data1);                              %debug
            xm(:,1)   = exp_data1(:,2);                                     % x, pressure (Pa)
            ym(:,1)   = exp_data1(1:(length(xm(:,1))),1);                   % y, q (mmol/g)
            largest_row = row1;
        
        elseif every_temp == 2
            exp_data2 = transpose(fitting_input.(every_gas_name).(T_name));
            first_row = exp_data2 (1,:);
            while isempty(exp_data2)==0
                 while (first_row(1) == 0 |first_row(2) == 0) && isempty(exp_data2) ==0
                    exp_data2(1,:) = [];
                    if size(exp_data2,1)>0
                        first_row = exp_data2 (1,:);
                    end
                 end
                 if first_row(1) ~= 0 && first_row(2) ~= 0
                    break
                 end
            end
           
            [row2, column2] = size(exp_data2);
            if  row2 < largest_row
                exp_data2(row2+1:row1,:) = 0;
            elseif row2 > largest_row
                exp_data1(row1+1:row2,:) = 0;
                largest_row = row2;
            end
            clearvars xm ym
            xm(:,1)  = exp_data1(:,2);                                     % x, pressure (Pa)
            ym(:,1)   = exp_data1(1:(length(xm(:,1))),1);                   % y, q (mmol/g)
            xm(:,2)   = exp_data2(1:(length(xm(:,1))),2); % x, pressure (Pa)
            ym(:,2)   = exp_data2(1:(length(xm(:,1))),1);  % y, q (mmol/g)

       elseif every_temp == 3
            exp_data3 = transpose( fitting_input.(every_gas_name).(T_name));
            first_row = exp_data3 (1,:);
            while isempty(exp_data3)==0
                while (first_row(1) == 0 |first_row(2) == 0 ) && isempty(exp_data3) ==0
               
                    exp_data3(1,:) = [];
                    if size(exp_data3,1)>0
                        first_row = exp_data3 (1,:);
                    end       
                end
                if first_row(1) ~= 0 && first_row(2) ~= 0
                    break
                end
            end
            [row3, column3] = size(exp_data3);
            if  row3 < largest_row
                exp_data3(row3+1:largest_row,:) = 0;
            elseif row3 > largest_row
                exp_data1(row1+1:row3,:) = 0;
                exp_data2(row2+1:row3,:) = 0;
                largest_row = row3;
            end
            clearvars xm ym 
            xm(:,1)  = exp_data1(:,2);                    % x, pressure (Pa)
            ym(:,1)   = exp_data1(1:(length(xm(:,1))),1); % y, q (mmol/g)
            xm(:,2)   = exp_data2(1:(length(xm(:,1))),2); % x, pressure (Pa)
            ym(:,2)   = exp_data2(1:(length(xm(:,1))),1);  % y, q (mmol/g)
            xm(:,3)   = exp_data3(1:(length(xm(:,1))),2); % x, pressure (Pa)
            ym(:,3)   = exp_data3(1:(length(xm(:,1))),1);  % y, q (mmol/g)
            
       elseif every_temp == 4
            exp_data4 = transpose( fitting_input.(every_gas_name).(T_name));
            first_row = exp_data4 (1,:);
            while isempty(exp_data4)==0
                while (first_row(1) == 0 |first_row(2) == 0 ) && isempty(exp_data4) ==0
                
                    exp_data4(1,:) = [];
                    if size(exp_data4,1)>0
                        first_row = exp_data4 (1,:);
                    end  
                end
                if first_row(1) ~= 0 && first_row(2) ~= 0
                    break
                end
            end
            [row4, column4] = size(exp_data4);
            if  row4 < largest_row
                exp_data4(row4+1:largest_row,:) = 0;
            elseif row4 > largest_row
                exp_data1(row1+1:row4,:) = 0;
                exp_data2(row2+1:row4,:) = 0;
                exp_data3(row3+1:row4,:) = 0;
                largest_row = row4;
            end
            clearvars xm ym
            xm(:,1)  = exp_data1(:,2);                                     % x, pressure (Pa)
            ym(:,1)   = exp_data1(1:(length(xm(:,1))),1);                   % y, q (mmol/g)
            xm(:,2)   = exp_data2(1:(length(xm(:,1))),2); % x, pressure (Pa)
            ym(:,2)   = exp_data2(1:(length(xm(:,1))),1);  % y, q (mmol/g)
            xm(:,3)   = exp_data3(1:(length(xm(:,1))),2); % x, pressure (Pa)
            ym(:,3)   = exp_data3(1:(length(xm(:,1))),1);  % y, q (mmol/g)
            xm(:,4)   = exp_data4(1:(length(xm(:,1))),2); % x, pressure (Pa)
            ym(:,4)   = exp_data4(1:(length(xm(:,1))),1);  % y, q (mmol/g)          
        
       elseif every_temp == 5
            exp_data5 = transpose( fitting_input.(every_gas_name).(T_name));
            first_row = exp_data5 (1,:);
            while isempty(exp_data5)==0 
                while (first_row(1) == 0 |first_row(2) == 0  ) && isempty(exp_data5) ==0
                    exp_data5(1,:) = [];
                    if size(exp_data5,1)>0
                        first_row = exp_data5 (1,:);
                    end      
                end
                if first_row(1) ~= 0 && first_row(2) ~= 0
                    break
                end
            end
            [row5, column5] = size(exp_data5); 
            if  row5 < largest_row
                exp_data5(row5+1:largest_row,:) = 0;
            elseif row5 > largest_row
                exp_data1(row1+1:row5,:) = 0;
                exp_data2(row2+1:row5,:) = 0;
                exp_data3(row3+1:row5,:) = 0;
                exp_data4(row4+1:row5,:) = 0;
                largest_row = row5;
            end
            clearvars xm ym
            xm(:,1)  = exp_data1(:,2);                                     % x, pressure (Pa)
            ym(:,1)   = exp_data1(1:(length(xm(:,1))),1);                   % y, q (mmol/g)
            xm(:,2)   = exp_data2(1:(length(xm(:,1))),2); % x, pressure (Pa)
            ym(:,2)   = exp_data2(1:(length(xm(:,1))),1);  % y, q (mmol/g)
            xm(:,3)   = exp_data3(1:(length(xm(:,1))),2); % x, pressure (Pa)
            ym(:,3)   = exp_data3(1:(length(xm(:,1))),1);  % y, q (mmol/g)
            xm(:,4)   = exp_data4(1:(length(xm(:,1))),2); % x, pressure (Pa)
            ym(:,4)   = exp_data4(1:(length(xm(:,1))),1);  % y, q (mmol/g)      
            xm(:,5)   = exp_data5(1:(length(xm(:,1))),2); % x, pressure (Pa)
            ym(:,5)   = exp_data5(1:(length(xm(:,1))),1);  % y, q (mmol/g)     
        
        elseif every_temp == 6
            exp_data6 = transpose( fitting_input.(every_gas_name).(T_name));
            first_row = exp_data6 (1,:);
            while isempty(exp_data6)==0
                while ( first_row(1) == 0 |first_row(2) == 0 ) && isempty(exp_data6) ==0
                    exp_data6(1,:) = [];
                    if size(exp_data6,1)>0
                        first_row = exp_data6 (1,:);
                    end             
                end
                if first_row(1) ~= 0 && first_row(2) ~= 0
                    break
                end
            end
            [row6, column6] = size(exp_data6); 
            if  row6 < largest_row
                exp_data6(row6+1:largest_row,:) = 0;
            elseif row6 > largest_row
                exp_data1(row1+1:row6,:) = 0;
                exp_data2(row2+1:row6,:) = 0;
                exp_data3(row3+1:row6,:) = 0;
                exp_data4(row4+1:row6,:) = 0;
                exp_data5(row5+1:row6,:) = 0;
                largest_row = row6;
            end
            clearvars xm ym
            xm(:,1)  = exp_data1(:,2);                                     % x, pressure (Pa)
            ym(:,1)   = exp_data1(1:(length(xm(:,1))),1);                   % y, q (mmol/g)
            xm(:,2)   = exp_data2(1:(length(xm(:,1))),2); % x, pressure (Pa)
            ym(:,2)   = exp_data2(1:(length(xm(:,1))),1);  % y, q (mmol/g)
            xm(:,3)   = exp_data3(1:(length(xm(:,1))),2); % x, pressure (Pa)
            ym(:,3)   = exp_data3(1:(length(xm(:,1))),1);  % y, q (mmol/g)
            xm(:,4)   = exp_data4(1:(length(xm(:,1))),2); % x, pressure (Pa)
            ym(:,4)   = exp_data4(1:(length(xm(:,1))),1);  % y, q (mmol/g)      
            xm(:,5)   = exp_data5(1:(length(xm(:,1))),2); % x, pressure (Pa)
            ym(:,5)   = exp_data5(1:(length(xm(:,1))),1);  % y, q (mmol/g) 
            xm(:,6)   = exp_data6(1:(length(xm(:,1))),2); % x, pressure (Pa)
            ym(:,6)   = exp_data6(1:(length(xm(:,1))),1);  % y, q (mmol/g) 
        end  
     end

    %% Defining constraints, boundaries, parameters 
    % Inequality constraint -> you don't need that
    A = [];
    b = [];

    % Equality constraint -> you don't need that
    Aeq = [];
    beq = [];
 
    switch model
        case 'langm_freund'
            T0 = 293;
            lb = [0,  0,  0,  0,  0,  -100]; % 1: n0 [mmol/g], 2: Xi [-], 3: t0 [-], 4: alpha [-], 5:b0 [1/Pa], 6:Q [kJ/mol]
            ub = [15, 15, 5, 10, 1, 100]; 

            p0 = [2.3, 0.5, 0.4, 0.5, 0.007, 1]; % starting point
            p0 = (p0-lb)./(ub-lb); % normalized

        case 'toth_cp'
            % extended toth isotherm including chemi- and physisorption
            T0 = 293;
            % Xi_c, dH_c (kJ/mol), alpha_c, Xi_p, dH_p (kJ/mol), alpha_p, ns0_c (mol/kg), b0_c (1/Pa), t0_c, ns0_p (mol/kg), b0_p (1/Pa), t0_p
            lb = [0        10      -2        0           1          0          0.5      1e-12     0        0.5      1e-010      0];
            ub = [10        130     2        10          50         2          15       10        5        18       7e-04       5];

            p0 = [5        60      0.6      0.011       5          0.43        10        0.6        0.1      1        5e-8        0.66 ];

            p0 = (p0-lb)./(ub-lb); % normalized

        case 's_shaped'
            %    q_L0(mol/kg), b_L0(1/Pa), dU_L(kJ/mol), q_U0(kJ/mol),b_U0(1/Pa), dU_U(kJ/mol), b_H0(mol/kg/Pa), dU_H(kJ/mol), Xi_1(-), Xi_2(-), pstep_0(Pa), dH_step(kJ/mol), gamma(-)
            T0 = 313.5;
            lb = [0.01      9*1e-11*1e-05   15          1           1*1e-12*1e-05    30         1*1e-08*1e-05   -50             0       0       1*1e-10*1e05    -120        1];
            ub = [3         0.1*1e-05       65          20           1*1e-04*1e-05    100        1*1e-01*1e-05   50              0.5     1e4     1*1e05          -20         10];

            p0 = [0.1     0.009*1e-05       31          2.4         9*10^(-7)*1e-05   59        5*10^(-4)*1e-05 18              0.124    0   5.0000e-04*1e05    -74.1       4];          
            p0 = (p0-lb)./(ub-lb);           
    end
    
    %% Define function and nonlinear constrain
            %p0 = (p0-lb)./(ub-lb); % normalts

    switch model
        case 'langm_freund'
            yp = langm_freund(xm,R,T0,T,lb,ub);
            nonlcon = [];
        case 'toth_cp'
            yp = toth_cp(xm,R,T0,T,lb,ub);
            nonlcon = [];
        case 's_shaped'
            yp = s_shaped(xm,R,T0,T,lb,ub);
            nonlcon = [];
    end

    %% check which parameters need to be considered for error calculation
    % delte rows with zero for the calculation of the error
    for k = 1:length(T)
        ymTmp = ym(:,k);
        indices1 = find(ymTmp==0);
        ymTmp(indices1,:) = [];
        ymSpec{:,k} = ymTmp;
        ymMean(k) = mean(ymTmp);
    end

    %% Error function
    switch error
        case 'ERRSQ'
            error_function = sum(sum((yp(p) - ym).^2));
        case 'HYBRID'
            error_function = 100/(length(xm1)-length(p0))*sum(sum(((ym-yp(p)).^2)./ym));
        case 'ARE'
            error_function = 100/(length(xm1))*sum(sum(abs((yp(p)-ym)./ym)));
        case 'EABS'
            error_function = sum(sum(abs((yp(p) - ym))));
        case 'MPSD'
            error_function = 100*sqrt(1/(length(xm1)-length(p0))*sum(sum(((ym-yp(p))./ym).^2)));
        case 'R2'
             if length(T)==2  
               objective = @(p) ((sum(sum((ym-yp(p)).^2)))./(sum(((sum((ymSpec{:,1}-ymMean(1)).^2)+sum((ymSpec{:,2}-ymMean(2)).^2))))));
             end

            if length(T)==3  
               objective = @(p) ((sum(sum((ym-yp(p)).^2)))./(sum(((sum((ymSpec{:,1}-ymMean(1)).^2)+sum((ymSpec{:,2}-ymMean(2)).^2)+sum((ymSpec{:,3}-ymMean(3)).^2))))));
            end
            
            if length(T)==4
               objective = @(p) ((sum(sum((ym-yp(p)).^2)))./(sum(((sum((ymSpec{:,1}-ymMean(1)).^2)+sum((ymSpec{:,2}-ymMean(2)).^2)+sum((ymSpec{:,3}-ymMean(3)).^2 + sum((ymSpec{:,4}-ymMean(4)).^2)))))));
            end
            
            if length(T)==5
               objective = @(p) ((sum(sum((ym-yp(p)).^2)))./(sum(((sum((ymSpec{:,1}-ymMean(1)).^2)+sum((ymSpec{:,2}-ymMean(2)).^2)+sum((ymSpec{:,3}-ymMean(3)).^2 + sum((ymSpec{:,4}-ymMean(4)).^2 + sum((ymSpec{:,5}-ymMean(5)).^2))))))));
            end
            
            if length(T)==6
               objective = @(p) ((sum(sum((ym-yp(p)).^2)))./(sum(((sum((ymSpec{:,1}-ymMean(1)).^2)+sum((ymSpec{:,2}-ymMean(2)).^2)+sum((ymSpec{:,3}-ymMean(3)).^2 + sum((ymSpec{:,4}-ymMean(4)).^2 + sum((ymSpec{:,5}-ymMean(5)).^2 + sum((ymSpec{:,6}-ymMean(6)).^2)))))))));
            end

        case 'rS'
            error_function = 1-((6*sum(sum((ym-yp(p)).^2)))./(length(xm1)*(length(xm1)-1).^2));
        case 'sRE'
            error_function = sqrt(sum(sum((((ym-yp(p))./ym)- (100/(length(xm1))*sum(sum(abs((yp(p)-ym)./ym))))).^2))/(length(xm1)-1))  ;
        case 'Xi'
           error_function = sum(sum(((yp(p) - ym).^2)./ym));
    end

    %% Find the optimum
    low = 0*ones((length(p0)),1)';
    up = 1*ones((length(p0)),1)';
    
    R2_check = 0;
    m = 1;
   
    while R2_check <= 0.9 && m < 3
        
        if m == 1
        else
            p0 = p;
        end
            
        switch optimizer
            case 'fmin'
                options = optimoptions('fmincon','Algorithm','sqp'); % sqp or interior-point  ,'PlotFcns',@optimplotfval
                options.MaxFunctionEvaluations = 600;
                 options.ConstraintTolerance = 1e-18;  % 1e-6
                 options.StepTolerance = 1e-18;        % 1e-6

                [p,fval,ef,output,lambda] = fmincon(objective,p0,A,b,Aeq,beq,low,up,nonlcon,options);

            case 'gs'
                %%% global search
                rng default
                gs = GlobalSearch;
                opts = optimoptions('fmincon','Algorithm','sqp'); % sqp or interior-point
                opts.StepTolerance = 1e-28;    
                problem = createOptimProblem('fmincon','objective',objective,'x0',p0,'lb',low,'ub',up,'options',opts,'nonlcon',nonlcon); 
                [p,fg,exitflag,output,solutions] = run(gs,problem);

            case 'pattern'
                %%% patternsearch
                options = optimoptions('fmincon','Algorithm','interior-point'); % sqp or interior-point
                p = patternsearch(objective,p0,A,b,Aeq,beq,low,up,nonlcon,options);    
        end

        popt = p;

        %% dimensioned
        popt = popt.*(ub-lb)+lb;

        % check
        for ii=1:length(p0)
            if popt(ii)==lb(ii)
                disp(['Parameter ' ,num2str(ii), ' at lower bound']);
            elseif popt(ii)==ub(ii)
                disp(['Parameter ' ,num2str(ii), ' at upper bound']);
            end
        end

        %% calculate error
        R2(m) = 1-(sum(sum((ym-yp(p)).^2)))/(sum(sum((ym-mean(ym)).^2)));
        
        if isnan(R2(m))
            R2(m) = 0;
        end
        
        param(m,:) = popt;
        param_norm(m,:) = p;
        
        R2_check = R2(m);
        lb_check(m,:) = lb;
        ub_check(m,:) = ub;
        m = m+1;
    end
    
    %% chose highest R2
    [R2_sorted,idx] = sortrows(R2');
    param_sorted = param(idx,:);
    lb_sorted = lb_check(idx,:);
    ub_sorted = ub_check(idx,:);
    param_norm_sorted = param_norm(idx,:);
    
    p_best = param_sorted(end,:);
    R2_best = R2_sorted(end);
    lb = lb_sorted(end,:);
    ub = ub_sorted(end,:);
    p = param_norm_sorted(end,:);

    %% calculate sigmaCO2

    switch model
        case 'langm_freund'
            sigmaCO2 = langm_freund(xm,R,T0,T,lb,ub);
        case 'toth_cp' 
            sigmaCO2 = toth_cp(xm,R,T0,T,lb,ub);
        case 's_shaped'
            sigmaCO2 = s_shaped(xm,R,T0,T,lb,ub);
    end
    sigmaCO2 = sigmaCO2(p);
    
    plot_isotherms_2 = "no";
    if plot_isotherms_2 == "yes"
        count_isotherms_plotted=0;
         figure()
         hold on
         colour = ['b', 'r','m','y', 'k', 'c'];
         for every_temp = [1:length(temperature_name_for_input)]
              count_isotherms_plotted =  count_isotherms_plotted +1;
              zero_list_to_append = zeros(1,size(xm,2));
              scatter( xm(1:length(ymSpec{:,every_temp}),every_temp), ym(1:length(ymSpec{:,every_temp}),every_temp), colour(every_temp) );
              xm = [zero_list_to_append; xm] ; ym = [zero_list_to_append ;ym ] ; sigmaCO2 = [zero_list_to_append; sigmaCO2];  
              specific_ymSpec = [0;ymSpec{:,every_temp}];
              plot(xm(1:length(specific_ymSpec),every_temp),sigmaCO2(1:length(specific_ymSpec),every_temp),'color', colour(every_temp), 'LineStyle', '-');
              temperature_to_display_pre = string(temperature_name_for_input(every_temp));
              temperature_name_to_display = strip(temperature_to_display_pre, 'left', 'T');
              temperature_name_to_display = strip(temperature_name_to_display, 'left', 'a');
              temperature_name_to_display = strip(temperature_name_to_display, 'left', 'A');
              temperature_name_to_display = strip(temperature_name_to_display, 'left', 'a');
              make_legend (:,every_temp) = [strcat('Data: ',temperature_name_to_display) strcat('Fitting: ', temperature_name_to_display)];
         end
        legend(make_legend, 'Location', 'bestoutside');
        axis([0 inf 0 inf]);
        box on
        xlabel('Pressure (Pa)');
        ylabel('Amount adsorbed (mmol/g)');
        fig = gcf;
        hold off
     end
    
    %% function models

    function yp = langm_freund(xm,R,T0,T,lb,ub)
        yp = @(p) ((p(1)*(ub(1)-lb(1))+lb(1)).*exp((p(2)*(ub(2)-lb(2))+lb(2)).*(1-T./T0))).*(((p(5)*(ub(5)-lb(5))+lb(5))*exp((p(6)*(ub(6)-lb(6))+lb(6))./(R*T0).*(T0./T-1))).*xm).^(1./(p(3)*(ub(3)-lb(3))+lb(3))+(p(4)*(ub(4)-lb(4))+lb(4)).*(1-T0./T))./(1+ (((p(5)*(ub(5)-lb(5))+lb(5))*exp((p(6)*(ub(6)-lb(6))+lb(6))./(R*T0).*(T0./T-1))).*xm).^(1./(p(3).*(ub(3)-lb(3))+lb(3))+(p(4)*(ub(4)-lb(4))+lb(4)).*(1-T0./T)));
    end

    function yp = toth_cp(xm,R,T0,T,lb,ub)
        yp = @(p) (((p(7)*(ub(7)-lb(7))+lb(7))*exp((p(1)*(ub(1)-lb(1))+lb(1))*(1-T/T0)).*((p(8)*(ub(8)-lb(8))+lb(8))*exp((p(2)*(ub(2)-lb(2))+lb(2))/(R*T0)*(T0./T-1))).*xm)./((1+(((p(8)*(ub(8)-lb(8))+lb(8))*exp((p(2)*(ub(2)-lb(2))+lb(2))/(R*T0)*(T0./T-1))).*xm).^((p(9)*(ub(9)-lb(9))+lb(9))+(p(3)*(ub(3)-lb(3))+lb(3))*(1-T0./T))).^(1./((p(9)*(ub(9)-lb(9))+lb(9))+(p(3)*(ub(3)-lb(3))+lb(3))*(1-T0./T)))))+(((p(10)*(ub(10)-lb(10))+lb(10))*exp((p(4)*(ub(4)-lb(4))+lb(4))*(1-T/T0)).*((p(11)*(ub(11)-lb(11))+lb(11))*exp((p(5)*(ub(5)-lb(5))+lb(5))/(R*T0)*(T0./T-1))).*xm)./((1+(((p(11)*(ub(11)-lb(11))+lb(11))*exp((p(5)*(ub(5)-lb(5))+lb(5))/(R*T0)*(T0./T-1))).*xm).^((p(12)*(ub(12)-lb(12))+lb(12))+(p(6)*(ub(6)-lb(6))+lb(6))*(1-T0./T))).^(1./((p(12)*(ub(12)-lb(12))+lb(12))+(p(6)*(ub(6)-lb(6))+lb(6))*(1-T0./T)))));
    end

    function yp = s_shaped(xm,R,T0,T,lb,ub)
        yp = @(p) ((p(1)*(ub(1)-lb(1))+lb(1)).*((p(2)*(ub(2)-lb(2))+lb(2)).*exp((p(3)*(ub(3)-lb(3))+lb(3))./(R.*T))).*xm./(1+((p(2)*(ub(2)-lb(2))+lb(2)).*exp((p(3)*(ub(3)-lb(3))+lb(3))./(R.*T))).*xm)).*(1-((exp((log(xm)-log((p(11)*(ub(11)-lb(11))+lb(11)).*exp(-(p(12)*(ub(12)-lb(12))+lb(12))./R.*(1/T0-1./T))))./((p(9)*(ub(9)-lb(9))+lb(9)).*exp((p(10).*(ub(10)-lb(10))+lb(10)).*(1./T0-1./T))))./(1+exp(((log(xm))-log((p(11)*(ub(11)-lb(11))+lb(11)).*exp(-(p(12)*(ub(12)-lb(12))+lb(12))./R.*(1/T0-1./T))))./((p(9)*(ub(9)-lb(9))+lb(9)).*exp((p(10).*(ub(10)-lb(10))+lb(10)).*(1./T0-1./T)))))).^(p(13)*(ub(13)-lb(13))+lb(13))))+((p(4)*(ub(4)-lb(4))+lb(4)).*((p(5)*(ub(5)-lb(5))+lb(5)).*exp((p(6)*(ub(6)-lb(6))+lb(6))./(R.*T))).*xm./(1+((p(5)*(ub(5)-lb(5))+lb(5)).*exp((p(6)*(ub(6)-lb(6))+lb(6))./(R.*T))).*xm)+((p(7)*(ub(7)-lb(7))+lb(7)).*exp((p(8)*(ub(8)-lb(8))+lb(8))./(R.*T))).*xm).*((exp((log(xm)-log((p(11)*(ub(11)-lb(11))+lb(11)).*exp(-(p(12)*(ub(12)-lb(12))+lb(12))./R.*(1/T0-1./T))))./((p(9)*(ub(9)-lb(9))+lb(9)).*exp((p(10).*(ub(10)-lb(10))+lb(10)).*(1./T0-1./T))))./(1+exp(((log(xm))-log((p(11)*(ub(11)-lb(11))+lb(11)).*exp(-(p(12)*(ub(12)-lb(12))+lb(12))./R.*(1/T0-1./T))))./((p(9)*(ub(9)-lb(9))+lb(9)).*exp((p(10).*(ub(10)-lb(10))+lb(10)).*(1./T0-1./T)))))).^(p(13)*(ub(13)-lb(13))+lb(13)));
    end

 end

 
