using Gencode, Preprocessing

function gmshcall(pde, filename, nd, elemtype)

opts = "-format msh3";

# find gmsh executable
gmsh = Preprocessing.findexec(pde.gmsh, pde.version);

print("Gmsh mesh generator...\n");
mystr = gmsh * " " * filename * ".geo -" * string(nd) * " " * opts;
run(Gencode.string2cmd(mystr));

tm = readlines(filename * ".msh");

nlines = length(tm);
j = 0;
for i = 1:nlines
    if tm[i] == "\$Nodes"
        j = i;
        break;
    end
end
np = parse(Int64,tm[j+1]);
p = zeros(nd,np);
i = j;
for ii = 1:np
    j = i+1+ii;
    pii = parse.(Float64,split(tm[j]));
    p[:,ii] = pii[2:(1+nd)];
end
i = j+1;
for ii = i:nlines
    if tm[ii] == "\$Elements"
        j = ii;
        break;
    end
end
ne = parse(Int64,tm[j+1]);
if nd==1 # line
    nve = 2;
    wcase = 1;
elseif nd==2 && elemtype==0 # tri
    nve = 3;
    wcase = 2;
elseif nd==2 && elemtype==1 # quad
    nve = 4;
    wcase = 3;
elseif nd==3 && elemtype==0 # tet
    nve = 4;
    wcase = 4;
elseif nd==3 && elemtype==1 # hex
    nve = 8;
    wcase = 5;
end
t = zeros(Int,nve,ne);
m = 0;
i = j;
for ii = 1:ne
    j = i+1+ii;
    tii = parse.(Int,split(tm[j]));
    if (tii[2] == wcase)
        m = m + 1;
        t[:,m] = tii[5:end];
    end
end
t = t[:,1:m];

return p, t

end
