function [p,t]=block2pt(X)

XT=block2tets(X);
dim = size(XT,1);
p0=reshape(XT,dim,[])';
t0=reshape(1:prod(size(XT))/dim,dim+1,size(XT,3))';
[p,t]=fixmesh(p0,t0);
