using DelimitedFiles

function writemesh(p, t, filename, mode)

nd,np = size(p);
nve,ne = size(t);
tmp = [nd, np, nve, ne];

fileID = open(filename, "w");
if mode==0 # binary
    write(fileID,Float64.(tmp));
    write(fileID,Float64.(p));
    write(fileID,Float64.(t));
else # text
    writedlm(fileID, tmp');
    writedlm(fileID, p');
    writedlm(fileID, t');
end
close(fileID);

end
