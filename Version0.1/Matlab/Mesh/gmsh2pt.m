function [p,t] = gmsh2pt(fname,elemtype)

fid = fopen(fname,'r');
if fid == -1
    error('Can''t open file');
end

% 1	2-node line	Edge Lagrange P1
% 2	3-node triangle	Triangle Lagrange P1
% 3	4-node quadrangle	Quadrilateral Lagrange P1
% 4	4-node tetrahedron	Tetrahedron Lagrange P1
% 5	8-node hexahedron	Hexahedron Lagrange P1
% 6	6-node prism	Wedge Lagrange P1
if (elemtype==2) 
    nd = 2; nn = 3;
elseif (elemtype==3)
    nd = 2; nn = 4;
elseif (elemtype==4)
    nd = 3; nn = 4;    
elseif (elemtype==5)
    nd = 3; nn = 8;    
end
    
readuntil(fid, '$Nodes');
np = fscanf(fid, '%d', 1);
p = zeros(np,nd);
for ii = 1:np
    foo = fscanf(fid, '%d', 1);
    p(ii,:) = fscanf(fid, '%f', nd)';
    fgetl(fid);
end

readuntil(fid, '$Elements');
nt = fscanf(fid, '%d', 1);
t = zeros(nt,nn);
for ii = 1:nt
    foo = fscanf(fid, '%d', 1);
    eltype = fscanf(fid, '%d', 1);
    if eltype == elemtype
        foo = fscanf(fid, '%d', 3);
        t(ii,:) = fscanf(fid, '%d', nn)';
    end
    t(ii,:)
    fgetl(fid);
    pause
end

t = t(any(t~=0,2),:);

fclose(fid);

end

function readuntil(fid,str)
    
    while ~feof(fid)
        fline = fgetl(fid);
        if ~isempty(fline) & ~isempty(strmatch(str, fline))
            break;
        end
    end
end
