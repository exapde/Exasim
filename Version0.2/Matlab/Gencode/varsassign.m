function mystr = varsassign(mystr, varname, n, flg, ustr)

if flg==0
    for i = 1:n
        str1 = varname + num2str(i);
        if any(contains(ustr, str1)) 
            str2 = varname + "[" + num2str(i-1) + "]";
            mystr = mystr + "\t\tT " + str1 + " = " + str2 + ";\n";
        end
    end
else
    for i = 1:n
        str1 = varname + num2str(i);
        if any(contains(ustr, str1)) 
            str2 = varname + "[" + num2str(i-1) + "*ng+i" + "]";
            mystr = mystr + "\t\tT " + str1 + " = " + str2 + ";\n";
        end
    end
end

end


% function mystr = varsassign(mystr, varname, n, flg)
% 
% if flg==0
%     for i = 1:n
%         str1 = varname + num2str(i);
%         str2 = varname + "[" + num2str(i-1) + "]";
%         mystr = mystr + "\t\tT " + str1 + " = " + str2 + ";\n";
%     end
% else
%     for i = 1:n
%         str1 = varname + num2str(i);
%         str2 = varname + "[" + num2str(i-1) + "*ng+i" + "]";
%         mystr = mystr + "\t\tT " + str1 + " = " + str2 + ";\n";
%     end
% end
% 
% end
% 
% % fid = fopen('test.cpp', 'w');
% % fprintf(fid, char(mystr));
% % fclose(fid);
% 
