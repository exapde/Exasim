import sympy

def sympyassign(mystr, f):

    n = len(f);
    #for i in range(0,n):
    #    f[i] = sympy.simplify(f[i]);

    ts, fs = sympy.cse(f);
    
    m = len(ts);
    for i in range(0,m):
        str1 = sympy.ccode(ts[i][0]);
        str2 = sympy.ccode(ts[i][1]);
        mystr = mystr + "\t\tT " + str1 + " = " + str2 + ";\n";

    for i in range(0,n):
        str1 = "\t\tf[" + str(i) + "*ng+i" + "]";
        str2 = sympy.ccode(fs[i]);
        mystr = mystr + str1 + " = " + str2 + ";\n";

    return mystr;
