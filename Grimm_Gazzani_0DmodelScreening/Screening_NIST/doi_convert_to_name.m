function new_doi = doi_convert_to_name(old_doi)
    %=====================================
    %Takes a DOI string and makes the DOI string suitable as a structure field name
    %--------------------------------------
    %Input:  - old_doi String to be converted)
    %Output: - new_doi String written in such a way that it can be used as
    %         field name in a structure
    %--------------------------------------

    old_doi = string(old_doi); 
    string_transform = strrep(old_doi,"/",'__');

    numbers = ["0" , "1", "2", "3", "4", "5" , "6", "7", "8", "9"];
    numbername = ["z", "o", "w", "h", "f", "v", "x", "s", "e", "n"];
    letters =     ["a","b", "c", "d", "e", "f", "g", "h" , "i" , "j", "k" , "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"];
    lettername = ["_a_","_b_", "_c_", "_d_", "_e_", "_f_", "_g_", "_h_" , "_i_" , "_j_", "_k_" , "_l_", "_m_", "_n_", "_o_", "_p_", "_q_", "_r_", "_s_", "_t_", "_u_", "_v_", "_w_", "_x_", "_y_", "_z_"];
    
    for letter =[1:length(letters)]
        string_transform = strrep(string_transform,letters(letter),lettername(letter));
    end
    
    for number = [1:10]
        string_transform = strrep(string_transform,numbers(number),numbername(number));
    end
    
    string_transform = strrep(string_transform,"/",'__');
    string_transform = strrep(string_transform,"-",'_');
    string_transform = strrep(string_transform,"(",'_op_');
    string_transform = strrep(string_transform,")",'_cl_');
    string_transform = strrep(string_transform,'.','d');
    new_doi = string_transform;
end