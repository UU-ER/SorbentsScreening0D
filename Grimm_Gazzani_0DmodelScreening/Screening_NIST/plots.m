function plots(R2_plot, characteristics_plot, dH_list, wc_list ,gamma_list, m0_list, q0_list, ns_equ_list, KH_list, list_of_R2,pCO2_1,material_list,specific_materials )
    %=====================================
    %Takes a DOI string and makes the DOI string suitable as a structure field name
    %--------------------------------------
    %Input:  - R2_plot If "R2_yes", makes swamchart of R2 distribution
    %        - characteristics_plot If "characteristics_yes" makes plots of
    %          the characteristics
    %        - dH_list List of all values of enthalpy of sdsorption
    %        - wc_list List of all values working capacity
    %        - gamma_list List of all values of non-linearity
    %        - m0_list List of all values of slope at feed concentration
    %        - q0_list List of all values of samount adsorbed at feed concentration
    %        - ns_equ_list List of all values of equilibrium loading
    %        - KH_list List of all values Henry coefficient
    %        - list_of_R2 List fo all values of R2
    %        - pCO2_1 Partial pessure of adsorption
    %        - material_list List of all materials
    %        - specific_materials List specific materials
    %Output: - Plots of the characteristics against each other.
    %--------------------------------------

%--- R2 distribution plot ---
		R2_plot = "false";
if R2_plot == "R2_yes"
    figure;
    swarmchart(zeros(1, length(list_of_R2)),list_of_R2, 'k' ,'filled'); 
    titlename = strcat("Distribution of the coefficient of determination. Partial pressure: ", string(pCO2_1), " Pa");
    title (titlename, 'FontSize', 12)
    box on
    ylabel ("Coefficient of determination")
    set(gca,'XTick',[])
end

%---characteristic plots ---
		characteristics_plot= "false";
if characteristics_plot == "characteristics_yes"
    %--- Working capacity distribution plot ---
    figure;
    swarmchart(zeros(1, length(wc_list)),wc_list,'k' ,'filled');
    titlename = strcat("Distribution the working capacity, partial pressure: ", string(pCO2_1), " Pa");
    title (titlename)
    box on
    ylabel ("Working capacity")
    set(gca,'XTick',[])

    %--- 3Dscatter plots of wc against dH and gamma ---
    figure;
    clearvars wc_larger_than_0 dH_list_wc_larger wc_smaller_than_0 dH_list_wc_smaller gamma_list_wc_larger gamma_list_wc_smaller
    wc_larger_than_0 = [];  wc_smaller_than_0 =[];  
    dH_list_wc_larger = []; dH_list_wc_smaller =[];
    gamma_list_wc_larger =[];gamma_list_wc_smaller =[];
    for i = [1:length(wc_list)]
        if wc_list(i)>=0
            wc_larger_than_0 =  [wc_larger_than_0 wc_list(i)];
            dH_list_wc_larger = [dH_list_wc_larger dH_list(i)];
            gamma_list_wc_larger = [gamma_list_wc_larger gamma_list(i)];
        else
           wc_smaller_than_0 =  [wc_smaller_than_0 wc_list(i)];
           dH_list_wc_smaller = [dH_list_wc_smaller dH_list(i)];
           gamma_list_wc_smaller =  [gamma_list_wc_smaller gamma_list(i)];
        end
    end
    s = scatter3(dH_list_wc_smaller, gamma_list_wc_smaller,wc_smaller_than_0, 'r' , 'filled');
    hold on
    if length(wc_larger_than_0)>0
        s = scatter3(dH_list_wc_larger, gamma_list_wc_larger  ,wc_larger_than_0, 'g' , 'filled');
    end

    title_string = strcat("Working capacity against the heat of adsorption and the nonlinearity, partial pressure: ", string(pCO2_1), " Pa");
    title(title_string)
    box on
    zlabel('working capacity')
    xlabel('Heat of adsorption (kJ/mol)') %dH
    ylabel('nonlinearity') %gamma

     %--- 3Dscatter plots of wc against q0 and m0 ---
    figure;
    clearvars wc_larger_than_0 m0_list_wc_larger wc_smaller_than_0 m0_list_wc_smaller q0_list_wc_larger q0_list_wc_smaller
    wc_larger_than_0 = [];  wc_smaller_than_0 =[];  
    m0_list_wc_larger = []; m0_list_wc_smaller =[];
    q0_list_wc_larger =[];q0_list_wc_smaller =[];
    for i = [1:length(wc_list)]
        if wc_list(i)>=0
            wc_larger_than_0 =  [wc_larger_than_0 wc_list(i)];
            m0_list_wc_larger = [m0_list_wc_larger m0_list(i)];
            q0_list_wc_larger = [q0_list_wc_larger q0_list(i)];
        else
           wc_smaller_than_0 =  [wc_smaller_than_0 wc_list(i)];
           m0_list_wc_smaller = [m0_list_wc_smaller m0_list(i)];
           q0_list_wc_smaller =  [q0_list_wc_smaller q0_list(i)];
        end
    end
    s = scatter3(m0_list_wc_smaller, q0_list_wc_smaller,wc_smaller_than_0, 'r' , 'filled');
    hold on
    if length(wc_larger_than_0)>0
        s = scatter3(m0_list_wc_larger, q0_list_wc_larger  ,wc_larger_than_0, 'g' , 'filled');
    end
    title_string = strcat("Working capacity against the local slope and CO2 loading at feed concentration, partial pressure: ", string(pCO2_1), " Pa");
    title(title_string)
    box on
    zlabel('working capacity')
    xlabel('Local slope at feed concentration') %m0
    ylabel('CO2 loading at feed concentration') %q0

    %--- 3Dscatter plots of wc against KH and ns_equ ---
    figure;
    clearvars wc_larger_than_0 KH_list_wc_larger wc_smaller_than_0 KH_list_wc_smaller ns_equ_list_wc_larger ns_equ_list_wc_smaller
    wc_larger_than_0 = [];  wc_smaller_than_0 =[];  
    KH_list_wc_larger = []; KH_list_wc_smaller =[];
    ns_equ_list_wc_larger =[];ns_equ_list_wc_smaller =[];
    for i = [1:length(wc_list)]
        if wc_list(i)>=0
            wc_larger_than_0 =  [wc_larger_than_0 wc_list(i)];
            KH_list_wc_larger = [KH_list_wc_larger KH_list(i)];
            ns_equ_list_wc_larger = [ns_equ_list_wc_larger ns_equ_list(i)];
        else
           wc_smaller_than_0 =  [wc_smaller_than_0 wc_list(i)];
           KH_list_wc_smaller = [KH_list_wc_smaller KH_list(i)];
           ns_equ_list_wc_smaller =  [ns_equ_list_wc_smaller ns_equ_list(i)];
        end
    end
    s = scatter3(KH_list_wc_smaller, ns_equ_list_wc_smaller,wc_smaller_than_0, 'r' , 'filled');
    hold on
    if length(wc_larger_than_0)>0
        s = scatter3(KH_list_wc_larger, ns_equ_list_wc_larger  ,wc_larger_than_0, 'g' , 'filled');
    end
    title_string = strcat("Working capacity against Henry's law coefficient and the equilibrium loading, partial pressure: ", string(pCO2_1), " Pa");
    title(title_string)
    box on
    zlabel('working capacity')
    xlabel('Equilibrium loading') %ns_equ
    ylabel("Henry's law coefficient")  %KH
    hold off
    
    %--- Heat of adsorption dH against Henry's coefficient KH ---
    figure;
    clearvars wc_larger_than_0 dH_list_wc_larger wc_smaller_than_0 dH_list_wc_smaller KH_list_wc_larger KH_list_wc_smaller
    wc_larger_than_0 = [];  wc_smaller_than_0 =[];  
    dH_list_wc_larger = []; dH_list_wc_smaller =[];
    KH_list_wc_larger = []; KH_list_wc_smaller =[];
    s = scatter(KH_list, dH_list, 'k' , 'filled');

    title_string = strcat("Heat of adsorption against Henry's coefficient, partial pressure: ", string(pCO2_1), " Pa");
    title(title_string);
    box on
    ylabel('Heat of adsorption (kJ/mol)') %dH
    xlabel("Henry's coefficient (mmol/(g Pa))") %gamma
    
    %--- Nonlinearity gamma against Henry's coefficient KH ---
    figure;
    clearvars wc_larger_than_0 gamma_list_wc_larger wc_smaller_than_0 gamma_list_wc_smaller gamma_list_wc_larger gamma_list_wc_smaller
    wc_larger_than_0 = [];  wc_smaller_than_0 =[];  
    gamma_list_wc_larger = []; gamma_list_wc_smaller =[];
    KH_list_wc_larger = []; KH_list_wc_smaller =[];
    s = scatter(KH_list, gamma_list, 'k' , 'filled');

    title_string = strcat("Non-linearity against Henry's coefficient, partial pressure: ", string(pCO2_1), " Pa");
    title(title_string);
    box on
    ylabel('Non-linearity') %dH
    xlabel("Henry's coefficient (mmol/(g Pa))") %gamma

    %--- Working capacity against the other characteristics ---
    figure;
    t = tiledlayout(3,2,'TileSpacing','Compact');
    title(t,strcat ("Working capacity against various characteristics, partial pressure: ", string(pCO2_1), " Pa"))
    
    %setting up colours
    colour_list_length = length(wc_list);
    colour_list = [];
    if length(specific_materials) ~= 0
        index_list = [];
        for every_special_mat = specific_materials
            for every_material_in_list = [1:length(material_list)]
                if startsWith(string(material_list(every_material_in_list)), every_special_mat) == 1
                    index_list = [index_list,every_material_in_list];     
                end
            end
        end
        all_colours = [1 0 0; 0 1 0; 0 1 1; 1 0 1; 1 1 0; 0 0 1; 0 0.4470 0.7410; 0.8500 0.3250 0.0980 ; 0.9290 0.6940 0.125;0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840];                              
        legend_list = strings(length(wc_list),1);
        for i = [1:length(index_list)]
            colour_list(index_list(i),:) = all_colours(i,:);
            legend_list(index_list(i),1) = material_list(index_list(i));
        end
        blue_index_list = [1: length(wc_list)];
        blue_index_list(index_list) = [];
        for j = blue_index_list
            colour_list(j,:) = [0 0 0];          %black
            legend_list(j,1) = 'Other';
        end
        
        %make sure that non-blue entries are printed last so these points
        %are visible
        [colour_list,I] = sortrows(colour_list); % sort colour list
        dH_list = dH_list(I);wc_list = wc_list(I); % sort the other lists accordingly
        gamma_list = gamma_list(I);m0_list = m0_list(I);q0_list = q0_list(I);
        ns_equ_list = ns_equ_list(I);legend_list = legend_list(I); KH_list = KH_list(I);
        
    else
        for i = [1:colour_list_length]
            colour_list = [colour_list; 0 0 0];  %black
        end
    end
  
    nexttile
    scatter(dH_list,wc_list, 20, colour_list, 'filled');
    titlename = strcat("Heat of adsorption");
    title (titlename)
    box on
    axis([0 inf 0 inf]);
    xlabel('Heat of adsorption (kJ/mol)') %dH
    ylabel ("Working capacity (mmol/g)")

    nexttile
    scatter(gamma_list,wc_list, 20, colour_list, 'filled'); 
    titlename = strcat("Non-linearity");
    title (titlename)
    box on
    axis([0 inf 0 inf]);
    ylabel("Working capacity (mmol/g)") 
    xlabel ("Non-linearity") %gamma

    nexttile
    scatter(KH_list,wc_list, 20, colour_list, 'filled');
    titlename = strcat("Henry's law coefficient");
    title (titlename)
    box on
    axis([0 inf 0 inf]);
    ylabel("Working capacity (mmol/g)") 
    xlabel("Henry's law coefficient (mmol/(g Pa))") %KH_list
    
    
    nexttile
    scatter(m0_list,wc_list, 20,  colour_list, 'filled');
    titlename = strcat("Local slope at feed concentration");
    title (titlename)
    box on
    axis([0 inf 0 inf]);
    ylabel("Working capacity (mmol/g)") 
    xlabel('Local slope at feed concentration (mmol/(g Pa))') %m0

    nexttile
    scatter(q0_list,wc_list, 20,  colour_list, 'filled');
    titlename = strcat("CO2 loading at feed concentration");
    title (titlename)
    box on
    axis([0 inf 0 inf]);
    ylabel("Working capacity (mmol/g)") 
    xlabel('CO2 loading at feed concentration (mmol/g)') %q0
    
    nexttile
    if length(specific_materials) ~= 0
        colour_list = unique(colour_list,'rows', 'stable');     
        gscatter(ns_equ_list,wc_list, legend_list, colour_list);
        legend('location','SouthOutside')
        legend('show');
    else
        scatter(ns_equ_list,wc_list, [], colour_list, 'filled');
    end    
    titlename = strcat("Equilibrium loading");
    title (titlename)
    box on
    axis([0 inf 0 inf]);
    ylabel("Working capacity (mmol/g)") 
    xlabel('Equilibrium loading (mmol/g)') %ns_equ_list
    
end

	Figure_new = false;
if Figure_new
%% new figures (without marking special materials)
    figure;
    t = tiledlayout(3,2,'TileSpacing','Compact');
    title(t,strcat ("Working capacity against over characteristics"))
    
    nexttile
    scatter(dH_list,wc_list, 20, colour_list, 'filled');
    titlename = strcat("Heat of adsorption");
    title (titlename)
    box on
    axis([0 inf 0 inf]);
    xlabel('Heat of adsorption (kJ/mol)') %dH
    ylabel ("Working capacity (mmol/g)")
end    


end
