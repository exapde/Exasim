function div = divergence(udg, href)

r    = udg(:,1,:);
ru   = udg(:,2,:);
rv   = udg(:,3,:);
rx   = udg(:,5,:);
rux  = udg(:,6,:);
ry   = udg(:,9,:);    
rvy  = udg(:,11,:);    

r1   = 1./r;
u    = ru.*r1;
v    = rv.*r1;

ux  = (rux - rx.*u).*r1;    
vy  = (rvy - ry.*v).*r1;            
div = href*(ux+vy);    

