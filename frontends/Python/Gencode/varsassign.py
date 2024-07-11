def varsassign(mystr, varname, n, flg, ustr):
    if flg == 0:
        for i in range(1, n+1):
            str1 = varname + str(i)
            if str1 in ustr:
                str2 = varname + "[" + str(i-1) + "]"
                mystr += "\t\tT " + str1 + " = " + str2 + ";\n"
    else:
        for i in range(1, n+1):
            str1 = varname + str(i)
            if str1 in ustr:
                str2 = varname + "[" + str(i-1) + "*ng+i" + "]"
                mystr += "\t\tT " + str1 + " = " + str2 + ";\n"
    
    return mystr

# def varsassign(mystr, varname, n, flg):

#     if flg==0:
#         for i in range(0,n):
#             str1 = varname + str(i+1);
#             str2 = varname + "[" + str(i) + "]";
#             mystr = mystr + "\t\tT " + str1 + " = " + str2 + ";\n";
#     else:
#         for i in range(0,n):
#             str1 = varname + str(i+1);
#             str2 = varname + "[" + str(i) + "*ng+i" + "]";
#             mystr = mystr + "\t\tT " + str1 + " = " + str2 + ";\n";

#     return mystr
