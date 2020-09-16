function varsassign(mystr::String, varname::String, n::Int, flg::Int)

if flg==0
    for i = 1:n
        str1 = varname * string(i);
        str2 = varname * "[" * string(i-1) * "]";
        mystr = mystr * "\t\tT " * str1 * " = " * str2 * ";\n";
    end
else
    for i = 1:n
        str1 = varname * string(i);
        str2 = varname * "[" * string(i-1) * "*ng+i" * "]";
        mystr = mystr * "\t\tT " * str1 * " = " * str2 * ";\n";
    end
end

return mystr;

end
