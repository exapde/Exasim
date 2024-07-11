function mesh = hdgmesh(mesh, porder)

nd = size(mesh.p,1);
nve = size(mesh.t,1);
ne = size(mesh.t, 2);

elemtype = 0;
if (nd==2) && (nve==4)
    elemtype=1;    
end
if (nd==3) && (nve==8)
    elemtype = 1;    
end
mesh.elemtype = elemtype;

[mesh.f, mesh.tprd, mesh.t2t] = facenumbering(mesh.p,mesh.t,elemtype,mesh.boundaryexpr,mesh.periodicexpr);

bcm = mesh.boundarycondition;
mesh.bf = mesh.f;
for j=1:length(bcm)
  ind = mesh.f==j;  % fix bug here
  mesh.bf(ind) = bcm(j);
end
f2t = mkf2e(mesh.t,elemtype,nd);

% reorder so that boundary faces are last
ina = find(f2t(3,:)>0); % interior faces
inb = find(f2t(3,:)==0); % boundary faces
inc = sub2ind(size(mesh.f), f2t(2,inb), f2t(1,inb));
fb = mesh.f(inc); % list of boundary indices

fa = unique(fb); % boundary indices    
bcn = unique(bcm(fa)); % a list of boundary conditions
nbc = length(bcn);

ind = zeros(1,length(fb));
m = 1;
for j=1:nbc % for each boundary condition bcn(j)
    bj = find(bcm==bcn(j)); % find all boundaries that have condition bcn(j)
    n = 0;
    for k = 1:length(bj) % for each boundary that has condition bcn(j)
        ii = find(fb == bj(k)); % indices of the boundary bj(k)
        l = length(ii);
        n = n + l;
        ind(m:(m+l-1)) = ii;
        m = m + l;
    end
end

% [interior faces, boundary faces]
f2t = f2t(:,[ina inb(ind)]);
mesh.f2t = f2t;
mesh.t2f = mke2f(f2t);
[mesh.facecon,mesh.elcon] = faceconnectivity2(mesh.t,f2t,nd,elemtype,porder);    
mesh.elcon = reshape(mesh.elcon, [], ne);

if isfield(mesh, "dgnodes")==0
  mesh.dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,porder);    
end

mesh.ne = ne;
mesh.nsiz = max(mesh.elcon(:));
mesh.bf = -mesh.bf;
