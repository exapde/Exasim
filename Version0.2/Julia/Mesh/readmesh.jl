function readmesh(filename, mode)

#fileID = fopen(filename,'r');
if mode==0 # binary
    tmp = reinterpret(Float64,read(filename));
    nd = Int(tmp[1]);
    np = Int(tmp[2]);
    nve = Int(tmp[3]);
    ne = Int(tmp[4]);
    p = reshape(tmp[5:(4+nd*np)], nd, np);
    t = reshape(tmp[(5+nd*np):(4+nd*np+nve*ne)], nve, ne);
    t = Int.(t);
else # text
    tmp = readdlm(filename);
    nd = Int(tmp[1,1]);
    np = Int(tmp[1,2]);
    nve = Int(tmp[1,3]);
    ne = Int(tmp[1,4]);
    p = Float64.(reshape(tmp[2:np+1,1:nd], np, nd));
    p = p';
    t = Int.(reshape(tmp[np+2:end,1:nve], ne, nve));
    t = t';
end

return p,t

end
