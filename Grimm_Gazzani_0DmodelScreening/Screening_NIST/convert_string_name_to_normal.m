function string_converted = convert_string_name_to_normal ( string_to_convert, sort)
    %=====================================
    %Takes a string and makes the string suitable as a structure name
    % (that is only contains letters of the alphabet, numbers and
    % underscores)
    %--------------------------------------
    %Input:  string_to_convert (string to be converted)
    %        sort what the string represents: either gas, name or temperature.
    %Output: sname_material written in such a way that it can be used as
    %        structure index
    %--------------------------------------
    sort = string(sort);
    string_transform = string ( string_to_convert);
    numbers = ["0" , "1", "2", "3", "4", "5" , "6", "7", "8", "9"];
    numbername = ["z", "o", "w", "h", "f", "v", "x", "s", "e", "n"];
    letters =     ["a","b", "c", "d", "e", "f", "g", "h" , "i" , "j", "k" , "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"];
    lettername = ["_a_","_b_", "_c_", "_d_", "_e_", "_f_", "_g_", "_h_" , "_i_" , "_j_", "_k_" , "_l_", "_m_", "_n_", "_o_", "_p_", "_q_", "_r_", "_s_", "_t_", "_u_", "_v_", "_w_", "_x_", "_y_", "_z_"];
    
    if sort =="gas"
        string_transform = string(regexprep(string_transform,'_',' '));
    elseif sort =="doi"    
        old_doi = string(string_transform); 
        string_transform = strrep(string_transform,'d','.'); 
        string_transform = strrep(string_transform,"_op_",'(');
        
        for number = [1:10]
            string_transform = strrep(string_transform,numbername(number),numbers(number));
        end   
        
        for letter =[1:length(letters)]
            string_transform = strrep(string_transform,lettername(letter),letters(letter));
        end   
    
        string_transform = strrep(string_transform,'_0_','z');
        string_transform = strrep(string_transform,"_1_",'o');
        string_transform = strrep(string_transform,"_2_",'w');
        string_transform = strrep(string_transform,'_3_','h');
        string_transform = strrep(string_transform,"_4_",'f');
        string_transform = strrep(string_transform,"_5_",'v');
        string_transform = strrep(string_transform,'_6_','x'); 
        string_transform = strrep(string_transform,"_7_",'s');
        string_transform = strrep(string_transform,"_8_",'e');
        string_transform = strrep(string_transform,'_9_','n');  
        
        string_transform = strrep(string_transform,"_cl_",')');    
        string_transform = strrep(string_transform,'_._','d');
        string_transform = strrep(string_transform,"__",'/');
        string_transform = strrep(string_transform,"_",'-');
    else
        string_transform = strrep(string_transform,"_apostrophe_","'");
        string_transform = strrep(string_transform,'_da_','-');
        string_transform = strrep(string_transform,'_d_','.');
        string_transform = strrep(string_transform,'_at_','@');
        string_transform = strrep(string_transform,'_about_','~');
        string_transform = strrep(string_transform,'_pr_','%');
        string_transform = strrep(string_transform,'_x_','*');
        string_transform = strrep(string_transform,'_roof_','^');
        string_transform = strrep(string_transform,'_sqop_','[');
        string_transform = strrep(string_transform,'_sqcl_',']');
        string_transform = strrep(string_transform,'_cuop_','{');
        string_transform = strrep(string_transform,'_cucl_','}');
        string_transform = strrep(string_transform,'_op_','(');
        string_transform = strrep(string_transform,'_cl_',')');
        string_transform = strrep(string_transform,'_slash_','/');
        string_transform = strrep(string_transform,'_plus_','+');
         string_transform = strrep(string_transform,'_cross_','×');
        string_transform = strrep(string_transform,'_minus_','-');
        string_transform = strrep(string_transform,'_equal_','=');
        string_transform = strrep(string_transform,'_oumlaut_','ö');
        string_transform = strrep(string_transform,'_times_','·'); % orignally 3 types but that cannot be converted  ⋅ ··
        string_transform = strrep(string_transform,'_comma_',',');
        string_transform = strrep(string_transform,'_IJ_','Ĳ' );
        string_transform = strrep(string_transform,'_dollar_', '$');
        string_transform = strrep(string_transform,'aAa','');
        string_transform = string(regexprep(string_transform,'_',' '));
    end
    string_converted = string_transform;
end