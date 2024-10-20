function string2cmd(str::String)

ind = findall(" ", str);
n = length(ind);
cmdstr = Array{String,1}(undef,n+1);
cmdstr[1] = str[1:(ind[1][1]-1)];
for i = 2:n
    cmdstr[i] = str[(ind[i-1][1]+1):(ind[i][1]-1)];
end
i = n+1;
cmdstr[i] = str[(ind[i-1][1]+1):end];

return Cmd(cmdstr)

end
