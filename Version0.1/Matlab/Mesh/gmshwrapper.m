function [p,t] = gmshwrapper(gmsh,filename,opts,nd,elemtype)

disp('Gmsh mesh generator...');

[gmshstatus0,~] = system("which " + gmsh);
[gmshstatus1,~] = system("which gmsh");
[gmshstatus2,~] = system("which /usr/bin/gmsh");
[gmshstatus3,~] = system("which /usr/local/bin/gmsh");
[gmshstatus4,~] = system("which /opt/local/bin/gmsh");        

if gmshstatus0==0
elseif gmshstatus1==0        
    gmsh = "gmsh";        
elseif gmshstatus2==0        
    gmsh = "/usr/bin/gmsh";    
elseif gmshstatus3==0        
    gmsh = "/usr/local/bin/gmsh";    
elseif gmshstatus4==0        
    gmsh = "/opt/local/bin/gmsh";    
else            
    error("Exasim search in /usr/bin, /usr/local/bin, and /opt/local/bin and could not find Gmsh. Please see the documentation to install it. After installation, please set its path to app.gmsh"); 
end

% call gmsh
str = ['!' char(gmsh) ' ' char(filename) '.geo -' char(num2str(nd)) ' ' char(opts)];
eval(str);

if (elemtype == 0) && (nd==2)
    wcase = 2;
elseif (elemtype == 1) && (nd==2)
    wcase = 3;
elseif (elemtype == 0) && (nd==3)
    wcase = 4;
elseif (elemtype == 1) && (nd==3)    
    wcase = 5;
end

[p,t] = gmsh2pt([char(filename) '.msh'], wcase);
p = p';
t = t';

end

function [p,t] = gmsh2pt(fname, wcase)

fid = fopen(fname,'r');
if fid == -1
    error('Cannot open file');
end

% 1	2-node line	Edge Lagrange P1
% 2	3-node triangle	Triangle Lagrange P1
% 3	4-node quadrangle	Quadrilateral Lagrange P1
% 4	4-node tetrahedron	Tetrahedron Lagrange P1
% 5	8-node hexahedron	Hexahedron Lagrange P1
% 6	6-node prism	Wedge Lagrange P1
if (wcase==2) 
    nd = 2; nn = 3;
elseif (wcase==3)
    nd = 2; nn = 4;
elseif (wcase==4)
    nd = 3; nn = 4;    
elseif (wcase==5)
    nd = 3; nn = 8;    
end
    
readuntil(fid, '$Nodes');
np = fscanf(fid, '%d', 1);
p = zeros(np,nd);
for ii = 1:np
    foo = fscanf(fid, '%d', 1);
    p(ii,1:nd) = fscanf(fid, '%f', nd)';
    fgetl(fid);
end

readuntil(fid, '$Elements');
nt = fscanf(fid, '%d', 1);
t = zeros(nt,nn);
for ii = 1:nt
    foo = fscanf(fid, '%d', 1);
    eltype = fscanf(fid, '%d', 1);
    if eltype == wcase
        foo = fscanf(fid, '%d', 2);
        t(ii,:) = fscanf(fid, '%d', nn)';
    end
    fgetl(fid);
end

t = t(any(t~=0,2),:);

fclose(fid);

end

function readuntil(fid,str)
    
    while ~feof(fid)
        fline = fgetl(fid);
        if ~isempty(fline) && ~isempty(strmatch(str, fline))
            break;
        end
    end
end
