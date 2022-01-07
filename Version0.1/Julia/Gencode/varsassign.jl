function varsassign(mystr::String, varname::String, n::Int, flg::Int, ustr)

if flg==0
    for i = 1:n
        str1 = varname * string(i);
        if contains(ustr, str1)
            str2 = varname * "[" * string(i-1) * "]";
            mystr = mystr * "\t\tT " * str1 * " = " * str2 * ";\n";
        end
    end
else
    for i = 1:n
        str1 = varname * string(i);
        if contains(ustr, str1)
            str2 = varname * "[" * string(i-1) * "*ng+i" * "]";
            mystr = mystr * "\t\tT " * str1 * " = " * str2 * ";\n";
        end
    end
end

return mystr;

end
