import sympy

def getccode(f, varstr):
    n = len(f)
    
    # Perform common subexpression elimination
    ts, fs = sympy.cse(f)
    
    str_code = ""
    m = len(ts)
    
    for i in range(m):
        str1 = sympy.ccode(ts[i][0])
        str2 = sympy.ccode(ts[i][1])
        str_code += "\t\tdstype {} = {};\n".format(str1, str2)

    for i in range(n):
        str1 = "\t\t{}[{}*ng+i]".format(varstr, i)
        str2 = sympy.ccode(fs[i])
        str_code += "{} = {};\n".format(str1, str2)

    # if len(fs) == 1:
    #     for i in range(n):
    #         str1 = "\t\t{}[{}*ng+i]".format(varstr, i)
    #         str2 = sympy.ccode(fs[0][i])
    #         str_code += "{} = {};\n".format(str1, str2)
    # elif len(fs) > 1:
    #     for i in range(n):
    #         str1 = "\t\t{}[{}*ng+i]".format(varstr, i)
    #         str2 = sympy.ccode(fs[i])
    #         str_code += "{} = {};\n".format(str1, str2)
    # else:
    #     raise ValueError("Symbolic expressions are wrong")
    
    return str_code

# Example usage:
# f = sympy.symbols('f0:10')
# varstr = "f"
# code = getccode(f, varstr)
# print(code)
