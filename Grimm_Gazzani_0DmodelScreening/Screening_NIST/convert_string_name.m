function string_converted = convert_string_name (string_to_convert, sort)
    %=====================================
    %Takes a string and makes the string suitable as a structure field name
    %--------------------------------------
    %Input:  - string_to_convert String to be converted)
    %        - sort What the string represents: either gas, name or temperature.
    %Output: -string_converted String written in such a way that it can be used as
    %         fieldname for a structure
    %--------------------------------------
    sort = string(sort);
    string_transform = string ( string_to_convert);
    string_transform = strip(string_transform, '''');
    string_transform = string(regexprep(string_transform,' ','_'));
    string_transform = strrep(string_transform,"'",'_apostrophe_');
    string_transform = strrep(string_transform,'.','_d_');
    string_transform = strrep(string_transform,'-','_da_');
    string_transform = strrep(string_transform,'(','_op_');
    string_transform = strrep(string_transform,')','_cl_');
    string_transform = strrep(string_transform,'@','_at_');
    string_transform = strrep(string_transform,'~','_about_');
    string_transform = strrep(string_transform,'%','_pr_');
    string_transform = strrep(string_transform,'*','_x_');
    string_transform = strrep(string_transform,'^','_roof_');
    string_transform = strrep(string_transform,'[','_sqop_');
    string_transform = strrep(string_transform,']','_sqcl_');
    string_transform = strrep(string_transform,'{','_cuop_');
    string_transform = strrep(string_transform,'}','_cucl_');
    string_transform = strrep(string_transform,'/','_slash_');
    string_transform = strrep(string_transform,'+','_plus_');
    string_transform = strrep(string_transform,'×','_cross_');
    string_transform = strrep(string_transform,'-','_minus_');
	string_transform = strrep(string_transform,'–','_minus_');
    string_transform = strrep(string_transform,'=','_equal_');
    string_transform = strrep(string_transform,'ö','_oumlaut_');
    string_transform = strrep(string_transform,'μ','_greekmu_');
    string_transform = strrep(string_transform,'⋅','_times_');
    string_transform = strrep(string_transform,'·','_times_');
    string_transform = strrep(string_transform,'·','_times_');
    string_transform = strrep(string_transform,',','_comma_');
    string_transform = strrep(string_transform,'Ĳ', '_IJ_');  
    string_transform = strrep(string_transform,'$', '_dollar_');
    numbers = [string(0), string(1) ,string(2) , string(3),string(4) ,string(5) , string(6), string(7), string(8),string(9)];
    if startsWith(string_transform, numbers)  | startsWith(string_transform, '_')
        string_transform = strcat('aAa',string_transform );
    end
    string_converted = string_transform;
end