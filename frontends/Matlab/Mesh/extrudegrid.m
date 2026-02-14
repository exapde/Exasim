function [p,t] = extrudegrid(p, t, zz)

nz = length(zz)-1;
nxy = size(p',1);

pz = repmat(zz,[nxy 1]);
p = [repmat(p',[nz+1 1]) pz(:)]';

t2d = t';
[ne2d, nv2d] = size(t2d);
t = zeros(ne2d*nz, nv2d*2);
for i = 1:nz
    t(ne2d*(i-1)+1:ne2d*i,:) = [t2d+(i-1)*nxy t2d+i*nxy];    
end
t = t';

