function getccode(f, varstr)

f = f[:];
n = length(f);
# for i=1:n
#     f[i] = sympy.simplify(f[i]);
# end
ts, fs = sympy.cse(f);

str = "";
m = length(ts);
for i = 1:m
    str1 = string(sympy.ccode(ts[i][1]));
    str2 = string(sympy.ccode(ts[i][2]));
    str = str * "\t\tT " * str1 * " = " * str2 * ";\n";
end

if length(fs)==1
    for i = 1:n
        #str1 = "\t\tf[" * string(i-1) * "*ng+i" * "]";
        #strj = [varstr num2str(j) '*ng+i] = 0.0;'];
        str1 = "\t\t" * varstr * string(i-1) * "*ng+i" * "]";
        str2 = string(sympy.ccode(fs[1][i]));
        str = str * str1 * " = " * str2 * ";\n";
    end
elseif length(fs)>1
    for i = 1:n
        #str1 = "\t\tf[" * string(i-1) * "*ng+i" * "]";
        str1 = "\t\t" * varstr * string(i-1) * "*ng+i" * "]";
        str2 = string(sympy.ccode(fs[i]));
        str = str * str1 * " = " * str2 * ";\n";
    end
else
    error("Symbolic expressions are wrong");
end

return str;

end
