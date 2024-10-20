function fhdg = hdgface(mesh,a)

if nargin<2 
    a = 0.1;
end
if isempty(a)
    a = 0.1;
end

he = meshsize(mesh);
he = round(he/1e-10)*1e-10;
hemin = min(he,[],1);
hestd = sort(hemin);
nhdg = round(a*mesh.ne);
h = hestd(nhdg);
ihdg = hemin<=h;
fhdg = mesh.t2f(ihdg,:);
fhdg = unique(fhdg(:));








