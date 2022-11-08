function [p,t] = transitionlayer2to1(xv, yv, zv)

nx = length(xv)-1;
ny = length(yv)-1;
nz = length(zv)-1;

if ny<3
    y1 = yv(1);
    y2 = yv(2);
else
    error("yv must have at least two entries.");
end

x1 = 0; x2 = 1; z1 = 0; z2 = 1;
[p1,t1] = templateTL2ref1(x1, y1, z1, x2, y2, z2);
[p2,t2] = templateTL2ref2(x1, y1, z1, x2, y2, z2);
[p3,t3] = templateTL2ref3(x1, y1, z1, x2, y2, z2);
[p4,t4] = templateTL2ref4(x1, y1, z1, x2, y2, z2);

nt = 0;
for i = 1:nx
    x1 = xv(i); x2 = xv(i+1);
    if rem(i,2)==1
        pi = p3;        
        ti = t3;
    else
        pi = p4;
        ti = t4;
    end
    pi(:,1) = x1 + pi(:,1)*(x2-x1);
    if i==1
        px = pi;
        tx = ti;
    else        
        px = [px; pi];        
        tx = [tx; ti+nt];
    end    
    nt = nt + size(pi,1); 
end
[px4,tx4]=fixmesh2(px,tx);

nt = 0;
for i = 1:nx
    x1 = xv(i); x2 = xv(i+1);
    if rem(i,2)==1
        pi = p1;        
        ti = t1;
    else
        pi = p2;
        ti = t2;
    end
    pi(:,1) = x1 + pi(:,1)*(x2-x1);
    if i==1
        px = pi;
        tx = ti;
    else        
        px = [px; pi];        
        tx = [tx; ti+nt];
    end    
    nt = nt + size(pi,1); 
end
[px3,tx3]=fixmesh2(px,tx);

nt = 0;
for i = 1:nz
    z1 = zv(i); z2 = zv(i+1);
    if rem(i,2)==1
        pi = px3;        
        ti = tx3;
    else
        pi = px4;
        ti = tx4;
    end
    pi(:,3) = z1 + pi(:,3)*(z2-z1);
    if i==1
        p = pi;
        t = ti;
    else        
        p = [p; pi];        
        t = [t; ti+nt];
    end    
    nt = nt + size(pi,1); 
end
[p,t]=fixmesh2(p,t);

% ptplot(p,t);
% xlabel('x');
% ylabel('y');
% zlabel('z');
% axis on;
