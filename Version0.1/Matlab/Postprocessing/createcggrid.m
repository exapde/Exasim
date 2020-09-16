function [cgnodes, cgelcon, cgcells, cell_t] = createcggrid(dgnodes,tlocal)
    
% Get some dimensions of the mesh
[~, nd, ne] = size(dgnodes);
[~, nve]     = size(tlocal);


[cgnodes,cgelcon,~,~,~] = mkcgent2dgent(dgnodes,1e-7);
% np = size(cgnodes,1);
% if nd==2
%     cgnodes = cat(2,cgnodes,zeros(np,1));
% end
cgcells = createcgcells(cgelcon,tlocal,ne);

elemtype = 0;
if (nd==2) && (nve==4)
    elemtype=1;    
end
if (nd==3) && (nve==8)
    elemtype=1;    
end

cell_t = getcelltype(nd,elemtype);
