%=====================================
% Checks the units; converts them if possible. Data from a DOI is only
% fitted if all units are convertible. The convertible units are defined in
% fit_nist_3 (list_convertible_units). 
%--------------------------------------

flag_specific_isotherm_fitted = 0;

if isfloat(unit_pressure) == 1
    flag_not_fitted = 1; 
end
if isfloat(unit_amount_adsorbed) == 1
    flag_not_fitted = 1 ; 
end

if flag_not_fitted ==0; 
    if unit_pressure == "bar"
        pressure = pressure * 10^5;
        count_bar = count_bar +1;
    elseif unit_pressure == "Pa" 
        pressure = pressure;
        count_Pa = count_Pa + 1;
    elseif unit_pressure == "Pascal"
        pressure = pressure;
        count_pascal = count_pascal +1;
    else
        pressure = pressure;
        unit_pressure %to spot other units than the ones above
    end

    
    if unit_amount_adsorbed == "mmol/g"
        amount_adsorbed = amount_adsorbed;
         count_mmolg =  count_mmolg + 1;
    elseif unit_amount_adsorbed == "cm3(STP)/g" 
        amount_adsorbed = 10^2 * amount_adsorbed / (8.314472 * 293); %10000 (Pa=J/m^-3) * aa *10^-6 (m^3) / (R (J K-1 mol-1) * T (K))) gives mol/g, *1000 gives mmol/g
        count_cm3stpg = count_cm3stpg + 1;
        
    elseif unit_amount_adsorbed == "wt%"
        amount_adsorbed = amount_adsorbed; %Not in list_convertible_units
        count_wt = count_wt + 1;
   elseif unit_amount_adsorbed == "ml/g"
        amount_adsorbed = amount_adsorbed; %Not in list_convertible_units
        count_mlg = count_mlg + 1;
    elseif unit_amount_adsorbed == "mol/mol"
        amount_adsorbed = amount_adsorbed; %Not in list_convertible_units
        count_molmol = count_molmol + 1;
   elseif unit_amount_adsorbed == "mmol/kg"
       amount_adsorbed = amount_adsorbed*10^-3;
        count_mmolkg = count_mmolkg + 1;
  elseif unit_amount_adsorbed == "molecules/unitcell"
       amount_adsorbed = amount_adsorbed; %Not in list_convertible_units
        count_molecules_unitcell = count_molecules_unitcell + 1;
   elseif unit_amount_adsorbed == "mg/g"
       amount_adsorbed = amount_adsorbed; %Not in list_convertible_units
        count_mgg = count_mgg + 1;
   elseif unit_amount_adsorbed == "g/g" 
       amount_adsorbed = amount_adsorbed;%Not in list_convertible_units
        count_gg = count_gg + 1;
   elseif unit_amount_adsorbed == "cm3(STP)/cm3" 
       amount_adsorbed = amount_adsorbed; %Not in list_convertible_units
        count_cm3STPcm3 = count_cm3STPcm3 + 1;
   elseif unit_amount_adsorbed == "mmol/m2" 
       amount_adsorbed = amount_adsorbed;%Not in list_convertible_units
        count_mmolm2 = count_mmolm2 + 1;
   elseif unit_amount_adsorbed == "mol/g" 
       amount_adsorbed = amount_adsorbed * 10 ^3;
        count_molg = count_molg + 1;
  elseif unit_amount_adsorbed == "mol/m3" 
       amount_adsorbed = amount_adsorbed;%Not in list_convertible_units
        count_molm3 = count_molm3 + 1;
  elseif unit_amount_adsorbed == "g/cm3" 
       amount_adsorbed = amount_adsorbed;%Not in list_convertible_units
        count_gcm3 = count_gcm3 + 1;
  elseif unit_amount_adsorbed == "g/l" 
       amount_adsorbed = amount_adsorbed;%Not in list_convertible_units
       count_gL = count_gL + 1;
  elseif unit_amount_adsorbed == "molecules/pore" 
       amount_adsorbed = amount_adsorbed;%Not in list_convertible_units
       count_moleculespore = count_moleculespore + 1;   
    elseif unit_amount_adsorbed == "molecules/8 unit cells" 
       amount_adsorbed = amount_adsorbed; %Not in list_convertible_units
       count_moleculeseightunitcells = count_moleculeseightunitcells + 1 
    elseif unit_amount_adsorbed == "mmol/cm3"
       amount_adsorbed = amount_adsorbed;%Not in list_convertible_units
       count_mmolcm3 = count_mmolcm3 + 1;      
    elseif unit_amount_adsorbed == "kg/mol"
        for i = [1: length(amount_adsorbed)]
            j = 1 / amount_adsorbed(i); %mol/kg = mmol/g
            new_amount_adsorbed (i) = j;
        end
       amount_adsorbed = new_amount_adsorbed;
       count_kgmol = count_kgmol + 1;    
    elseif unit_amount_adsorbed == "mg/m2"
       amount_adsorbed =  amount_adsorbed; %Not in list_convertible_units
       count_mgm2 = count_mgm2 + 1;  
    else
        unit_amount_adsorbed %to spot other units than the ones above
        amount_adsorbed = amount_adsorbed;
    end
    
    %if the unit is the right unit, save it in list_unit_amount adsorbed
    if ismember(unit_amount_adsorbed, list_convertible_units) ==1
        flag_specific_isotherm_fitted =1; %only if this flag is 1 will the isotherm be fitted, else it will be ignored!
    else
        flag_specific_isotherm_fitted =0;
    end
    
end 