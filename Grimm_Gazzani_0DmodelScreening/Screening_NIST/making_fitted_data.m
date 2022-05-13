%=====================================
% Checks on fitting error, computes characteristics.
%--------------------------------------
% Check if R2 is within boundaries
if R2 >= 0 && R2 <= 1
   flag_fitted = 1;
   flag_not_fitted = 0;
   Counter_d_step_8_9_13_14 = Counter_d_step_8_9_13_14 + 1; 
   fitted_data.(every_gas_name).(material_struc_name).(material_version_name).parameters = p;
   rhoMat = []; % kg/m3 --> not known
   rho = [];
   cp = [];
else
   Counter_d_step_8_13 = Counter_d_step_8_13 +1;
   flag_not_fitted = 1;
   rhoMat = []; % kg/m3 --> not known
   rho = [];
   cp = [];
end

if flag_fitted == 1 && flag_not_fitted == 0
   T_list = convert_T_to_numbers(T_list);
   T_list = sort(unique(T_list));
   %--- compute the characteristics ---
  [KH_save, ns_equ_save, q0_save, m0_save, gamma_save, wc_save,dH_save,T0] = characteristics_4(p,model,pCO2_1,rhoMat,[],120) ;
  if wc_save < 10e14
    flag_wc = 1;
    structure_exists = 0;
    %--- save the characteristics ---
    [dH_list, wc_list,gamma_list, m0_list, q0_list, ns_equ_list,KH_list,dataset,list_of_R2,material_names,p,model,material_number_name_normal ,T0,rho,rhoMat,cp] = assign_characteristics ("fitted_data", dataset, every_gas_name, material_struc_name, material_number_name_normal ,material_version_name, Counter_d_step_8_9_13_14, KH_save, ns_equ_save , q0_save, m0_save, gamma_save, wc_save,dH_save,dH_list, wc_list,gamma_list, m0_list, q0_list, ns_equ_list,KH_list,R2,list_of_R2,structure_exists,material_names,p,model,T0,rho,rhoMat,cp);
     %--- make data to plot ---
    for every_temp = [1:length(temperature_name_for_input)]
         values_to_put_in = double(fitting_input.Carbon_Dioxide.(temperature_name_for_input(every_temp))); 
         temp_struc =temperature_name_for_input(every_temp);
         structure_from_NIST.(every_gas_name).(material_struc_name).(material_version_name).(temp_struc).amount_adsorbed_NIST =  values_to_put_in(1,:); 
         structure_from_NIST.(every_gas_name).(material_struc_name).(material_version_name).(temp_struc).pressure_NIST =  values_to_put_in(2,:) ;
     end
  else
      flag_wc = 0;
  end
end