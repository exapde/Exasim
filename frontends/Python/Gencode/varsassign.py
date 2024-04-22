def varsassign(mystr, varname, n, flg):

    if flg==0:
        for i in range(0,n):
            str1 = varname + str(i+1);
            str2 = varname + "[" + str(i) + "]";
            mystr = mystr + "\t\tT " + str1 + " = " + str2 + ";\n";
    else:
        for i in range(0,n):
            str1 = varname + str(i+1);
            str2 = varname + "[" + str(i) + "*ng+i" + "]";
            mystr = mystr + "\t\tT " + str1 + " = " + str2 + ";\n";

    return mystr
