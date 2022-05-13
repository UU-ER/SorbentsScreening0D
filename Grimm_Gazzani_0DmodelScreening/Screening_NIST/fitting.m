%=====================================
% Initiates the fitting process, makes graphs if plot_isotherms == "yes"
%--------------------------------------

%--- check if material has more than 6 temperatures ---
if flag_doi_not_unique ==1
    if length(T_list) > 6 %because only 6 temperatures can be fitted
        %only delete out of range temperatures if list pure temperatures is too long
        T_list_unique = convert_T_to_numbers(T_list);
        index_to_del = []; % indexes for which T > 473, T <273
        if length(T_list_unique) > 6
            for every_temp = [1:length(T_list_unique)]
                if T_list_unique(every_temp) > 473 | T_list_unique(every_temp) < 273
                    index_to_del = [index_to_del, every_temp];
                end
            end
            T_list(index_to_del) = [];
            if flag_celsius == 1 
                T_list_in_Celsius(index_to_del)=[];
            end
        end
    end

    % if list longer than 6, some temperatures are deleted
    T_list = delete_temperatures_from_list (T_list);

    count_temperature_for_input = 0;
    clearvars temperature_name_for_input

    %--- make structure of data that will be fitted:fitting_data ---
    for every_temp = [1:length(T_list)]
        temperature_name = T_list(every_temp);
        flag_specific_isotherm_fitted = 0;
        if flag_celsius == 1 %need to convert because in structure, temperatures are saved as in their original form (e.g. Celsius)
            temperature_name_Kelvin = temperature_name;
            temperature_name_Celsius = T_list_in_Celsius(every_temp); %because in structure names are in original form (so Celsius)
            temperature_name = temperature_name_Celsius;
        end   
        if isempty(data_NIST.(material_struc_name).(every_gas_name).(temperature_name))
        else
            if ismember(data_NIST.(material_struc_name).(every_gas_name).(temperature_name).Doi_NIST, material_version_original)==1 
                    count_temperature_for_input = count_temperature_for_input + 1; %counts number of actual isotherms
                    clearvars pressure amount_adsorbed

                    sorted_data_NIST_2.(every_gas_name).(material_struc_name).(temperature_name) = data_NIST.(material_struc_name).(every_gas_name).(temperature_name);
                    sorted_data_NIST_2.(every_gas_name).(material_struc_name).(temperature_name).Full_name = material_full_name; 
                    pressure = sorted_data_NIST_2.(every_gas_name).(material_struc_name).(temperature_name).pressure_NIST;
                    amount_adsorbed = sorted_data_NIST_2.(every_gas_name).(material_struc_name).(temperature_name).amount_adsorbed_NIST; 

                    if first_sorting =="yes"
                        %--- check if units are correct ---
                        unit_pressure = sorted_data_NIST_2.(every_gas_name).(material_struc_name).(temperature_name).Unit_pressure;
                        unit_amount_adsorbed = sorted_data_NIST_2.(every_gas_name).(material_struc_name).(temperature_name).Unit_adsorption;
                        check_units
                        list_unit_amount_adsorbed = [list_unit_amount_adsorbed unit_amount_adsorbed];
                    else
                        list_unit_amount_adsorbed =["mmol/g"]
                    end

                    %The temperatures should be displayed in Kelvin 
                    if flag_celsius == 1
                        %convert Celsius to Kelvin
                         temperature_name = temperature_name_Kelvin;
                    end

                    if first_sorting == "no"
                        flag_specific_isotherm_fitted = 1;
                    end

                    if flag_specific_isotherm_fitted ==1 %only if specific isotherm has right units (see check_units)

                        try
                            fitting_input.(every_gas_name).(temperature_name) = [amount_adsorbed; pressure];
                        catch
                            flag_not_fitted = 1;
                            flag_error_fitting_input =1;
                            fprintf("Error: Problem making structure fitting input: probably amount_adsorbed and pressure not of the same length!")
                        end
                        temperature_name_for_input(count_temperature_for_input) = temperature_name;
                    end
            else
                flag_not_fitted = 1;
                flag_doi_unique =1;
            end
        end
    end
end

%--- Fitting the data ---
if flag_T_present == 1  
    Counter_d_step_3_4 = Counter_d_step_3_4 + 1;
    if flag_not_fitted == 0 
        Counter_d_step_4_5 = Counter_d_step_4_5 + 1; 
        if isempty(list_unit_amount_adsorbed)==0 
            Counter_d_step_5_6 = Counter_d_step_5_6 + 1;
            if all(ismember(list_unit_amount_adsorbed, list_convertible_units))== 1 
                Counter_d_step_6_7 = Counter_d_step_6_7 + 1;
                if isfield(fitting_input,(every_gas_name)) == 1 
                   material_number_name_normal = convert_string_name_to_normal ( material_full_name, "name");
                   every_gas_name_normal = convert_string_name_to_normal ( every_gas_name, "gas");
                   flag_fitted = 0;
                   clearvars R2
                   flag_count_dois_toth_s_shape =0;
                   %--- try to fit with the Toth, if R2 <0.98 try S-shaped ---
                   try 
                       [xm_t,ym_t,ymSpec_t,sigmaCO2_t,p_t,R2] = main_fitting_Nist_5(fitting_input,temperature_name_for_input,every_gas_name,'toth_cp');%''
                       R2_toth =R2;
                       model = 'toth_cp';
                       error_R2_toth = 1 - R2_toth;
                       flag_s_shape = 0;
                       
                       % try to fit with langmuir-freundlich model
                       [xm_t2,ym_t2,ymSpec_t2,sigmaCO2_t2,p_t2,R22] = main_fitting_Nist_5(fitting_input,temperature_name_for_input,every_gas_name,'langm_freund');%''
                       R2_langfr = R22;
                       % check if toth-cp or langmuir-fr is better
                       if R2_toth > R2_langfr
                           R2 = R2_toth;
                           model = 'toth_cp';
                       else 
                           R2 = R2_langfr;
                           model = 'langm_freund';
                           error_R2_toth = 1 - R2_langfr;
                           pause(2)
                       end
                       % Try to fit with s_shaped model
                       if error_R2_toth >= 0.02 
                           [xm_s,ym_s,ymSpec_s,sigmaCO2_s,p_s,R2] = main_fitting_Nist_5(fitting_input,temperature_name_for_input,every_gas_name,'s_shaped');
                           R2_s_shaped = R2;
                           error_R2_s_shaped = 1 - R2_s_shaped;
                           if error_R2_toth >= error_R2_s_shaped % take s_shape
                               model = 's_shaped';
                               flag_s_shape = 1;
                               R2 = R2_s_shaped;sigmaCO2 = sigmaCO2_s; p = p_s;
                               xm = xm_s; ym = ym_s; ymSpec = ymSpec_s;
                               
                               if R2_s_shaped < 0.8
                                   no_R2_low = no_R2_low+1;
                                   save_R2_low(no_R2_low).model = 's_shaped';
                                   save_R2_low(no_R2_low).R2 = error_R2_s_shaped;
                                   save_R2_low(no_R2_low).no = material_number;
                                   save_R2_low(no_R2_low).name = material_number_name_normal;
                                   save_R2_low(no_R2_low).doi = every_non_unique_doi;   
                               end
                                   
                           else % take toth or langmuir-fr
                               if R2_toth > R2_langfr
                                   model = 'toth_cp';
                                   flag_toth_cp = 1;
                                   xm = xm_t; ym = ym_t; ymSpec = ymSpec_t;
                                   R2 = R2_toth; sigmaCO2 = sigmaCO2_t; p = p_t;
                                   
                                   if R2_toth < 0.8
                                       no_R2_low = no_R2_low+1;
                                       save_R2_low(no_R2_low).model = 'toth';
                                       save_R2_low(no_R2_low).R2 = R2_toth;
                                       save_R2_low(no_R2_low).no = material_number;
                                       save_R2_low(no_R2_low).name = material_number_name_normal;
                                       save_R2_low(no_R2_low).doi = every_non_unique_doi;   
                                   end
                               else 
                                   model = 'langm_freund';
                                   flag_langmfr = 1;
                                   xm = xm_t2; ym = ym_t2; ymSpec = ymSpec_t2;
                                   sigmaCO2 = sigmaCO2_t2; p = p_t2;
                                   
                                   if R2_langfr < 0.8
                                       no_R2_low = no_R2_low+1;
                                       save_R2_low(no_R2_low).model = 'langfr';
                                       save_R2_low(no_R2_low).R2 = R2_langfr;
                                       save_R2_low(no_R2_low).no = material_number;
                                       save_R2_low(no_R2_low).name = material_number_name_normal;
                                       save_R2_low(no_R2_low).doi = every_non_unique_doi;   
                                   end
                               end
                           end
                       else % check which one is better
                           if R2_toth > R2_langfr
                               flag_toth_cp = 1;
                               xm = xm_t; ym = ym_t; ymSpec = ymSpec_t;
                               R2 = R2_toth; sigmaCO2 = sigmaCO2_t; p = p_t;
                           else 
                               flag_langmfr = 1;
                               xm = xm_t2; ym = ym_t2; ymSpec = ymSpec_t2;
                               sigmaCO2 = sigmaCO2_t2; p = p_t2;
                           end
                       end

                       if flag_s_shape ==1
                          count_dois_s_shape = count_dois_s_shape +1;
                       else
                           if flag_toth_cp ==1
                                count_dois_toth = count_dois_toth +1;
                           else
                               count_dois_langmfr = count_dois_langmfr +1;
                           end
                       end

                       flag_count_dois_toth_s_shape = 1;

                   catch
                        flag_not_fitted = 1;
                        R2 = 0;
                   end
                   
                    if first_sorting =="yes"                           
                        if  flag_not_fitted ==0 && all(ismember( list_unit_amount_adsorbed,list_convertible_units)) == 1
                            flag_wc = 0;
                            % Checks on fitting error, computes
                            % characteristics:
                            making_fitted_data 
                             if flag_wc == 1
                                count_plots = count_plots +1;
                                list_all_units_amount_adsorbed = [list_all_units_amount_adsorbed list_unit_amount_adsorbed];
                                 if plot_isotherms == "yes"
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
                                    title_name = strcat('Material: ',material_number_name_normal , "  Version: ",material_version_original,'  Gas: ', every_gas_name_normal );
                                    title(title_name);
                                    axis([0 inf 0 inf]);
                                    box on
                                    xlabel('Pressure (Pa)');
                                    ylabel('Amount adsorbed (mmol/g)');
                                    fig = gcf;
                                    hold off
                                 end
                             end
                        end
                    else
                        if  R2 >= 0 && R2 <= 1 && R2 > 0.96  
                         list_all_units_amount_adsorbed = [list_all_units_amount_adsorbed list_unit_amount_adsorbed];
                         flag_wc = 0;
                         % Checks on fitting error, computes
                         % characteristics:
                         making_fitted_data
                            if flag_wc == 1
                                count_plots = count_plots +1;
                                    if plot_isotherms == "yes" 
                                        figure()
                                        hold on
                                        colour = ['b', 'r','m','y', 'k', 'c'];
                                         for every_temp = [1:length(temperature_name_for_input)]
                                              count_isotherms_plotted =  count_isotherms_plotted +1;
                                              zero_list_to_append = zeros(1,size(xm,2));
                                              scatter( xm(1:length(ymSpec{:,every_temp}),every_temp), ym(1:length(ymSpec{:,every_temp}),every_temp), colour(every_temp) );%, 'o', 'LineStyle', 'none' );
                                              xm = [zero_list_to_append; xm] ; ym = [zero_list_to_append ;ym ] ; sigmaCO2 = [zero_list_to_append; sigmaCO2]; 
                                              plot(xm(1:length(ymSpec{:,every_temp}),every_temp),sigmaCO2(1:length(ymSpec{:,every_temp}),every_temp),'color', colour(every_temp), 'LineStyle', '-');
                                              temperature_to_display_pre = string(temperature_name_for_input(every_temp));
                                              temperature_name_to_display = strip(temperature_to_display_pre, 'left', 'T');
                                              temperature_name_to_display = strip(temperature_name_to_display, 'left', 'a');
                                              temperature_name_to_display = strip(temperature_name_to_display, 'left', 'A');
                                              temperature_name_to_display = strip(temperature_name_to_display, 'left', 'a');
                                              make_legend (:,every_temp) = [strcat('Data: ',temperature_name_to_display ) strcat('Fitting: ', temperature_name_to_display)];
                                         end
                                        legend(make_legend, 'Location', 'bestoutside');
                                         material_version_normal = convert_string_name_to_normal ( material_version, "doi"); %material_version_original?
                                        title_name = strcat('Material: ',material_number_name_normal , "  Version: ",material_version_normal,'  Gas: ', every_gas_name_normal );
                                        title(title_name);
                                        axis([0 inf 0 inf]); 
                                        xlabel('Pressure (Pa)');
                                        box on
                                        ylabel('Amount adsorbed (mmol/g)');
                                        fig = gcf;
                                        hold off
                                    end
                            end
                        end
                    end
                end
            end
        end
    end
end


                    
    
    
    
    
    
    
    