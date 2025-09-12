function mystr = symsassign2(mystr, f, udg, wdg, uhg)

str1 = getccode(f, 'f[');
str1 = "\t\t{\n" + str1 + "\t\t}\n";
mystr = mystr + str1;

if nargin>=3 && isempty(udg) == 0
  f_udg  = jacobian(f(:),udg);
  if isempty(f_udg) == 0
    str2 = getccode(f_udg, 'f_udg[');
    str2 = "\t\t{\n" + str2 + "\t\t}\n";
    mystr = mystr + str2;
  end    
end

if nargin>=4 && isempty(wdg) == 0
  f_wdg  = jacobian(f(:),wdg);
  if isempty(f_wdg) == 0
    str3 = getccode(f_wdg, 'f_wdg[');
    str3 = "\t\t{\n" + str3 + "\t\t}\n";
    mystr = mystr + str3;
  end
end

if nargin>=5 && isempty(uhg) == 0
  f_uhg  = jacobian(f(:),uhg);
  if isempty(f_uhg) == 0
    str4 = getccode(f_uhg, 'f_uhg[');
    str4 = "\t\t{\n" + str4 + "\t\t}\n";
    mystr = mystr + str4;
  end
end

end


