function po = outletpressure(udg, udg0, udg1, udg2, nl, dd, pout, alpha, gam, stage) 

gam1 = gam-1;

nd = length(nl);
if nd==2
    r0 = udg0(1);
    ru0 = udg0(2);
    rv0 = udg0(3);
    rE0 = udg0(4);
    u0 = ru0/r0;
    v0 = rv0/r0;
    q0 = 0.5*(u0*u0+v0*v0);
    p0 = gam1*(rE0-r0*q0);
    c0 = sqrt(gam*p0/r0);
    un0 = u0*nl(1)+v0*nl(2);
    
    r = udg(1);
    ru = udg(2);
    rv = udg(3);
    %rE = udg(4);    
    u = ru/r;
    v = rv/r;
%     q = 0.5*(u*u+v*v);
%     p = gam1*(rE-r*q);
%     c = sqrt(gam*p/r);
    un = u*nl(1)+v*nl(2);    
else
    r0 = udg0(1);
    ru0 = udg0(2);
    rv0 = udg0(3);
    rw0 = udg0(4);
    rE0 = udg0(5);    
    u0 = ru0/r0;
    v0 = rv0/r0;
    w0 = rw0/r0;
    q0 = 0.5*(u0*u0+v0*v0+w0*w0);
    p0 = gam1*(rE0-r0*q0);
    c0 = sqrt(gam*p0/r0);
    un0 = u0*nl(1)+v0*nl(2)+w0*nl(3);    
    
    r = udg(1);
    ru = udg(2);
    rv = udg(3);
    rw = udg(4);
%    rE = udg(5);    
    u = ru/r;
    v = rv/r;
    w = rw/r;
%     q = 0.5*(u*u+v*v+w*w);
%     p = gam1*(rE-r*q);
%     c = sqrt(gam*p/r);
    un = u*nl(1)+v*nl(2)+w*nl(3);        
end

if (stage==2)
if nd==2
    r1 = udg1(1);
    ru1 = udg1(2);
    rv1 = udg1(3);
    rE1 = udg1(4);
    u1 = ru1/r1;
    v1 = rv1/r1;
    q1 = 0.5*(u1*u1+v1*v1);
    p1 = gam1*(rE1-r1*q1);
    %c1 = sqrt(gam*p1/r1);
    un1 = u1*nl(1)+v1*nl(2);    
else
    r1 = udg1(1);
    ru1 = udg1(2);
    rv1 = udg1(3);
    rw1 = udg1(4);
    rE1 = udg1(5);    
    u1 = ru1/r1;
    v1 = rv1/r1;
    w1 = rw1/r1;
    q1 = 0.5*(u1*u1+v1*v1+w1*w1);
    p1 = gam1*(rE1-r1*q1);
    %c1 = sqrt(gam*p1/r1);
    un1 = u1*nl(1)+v1*nl(2)+w1*nl(3);        
end
end


if (stage==3)
if nd==2
    r2 = udg2(1);
    ru2 = udg2(2);
    rv2 = udg2(3);
    rE2 = udg2(4);
    u2 = ru2/r2;
    v2 = rv2/r2;
    q2 = 0.5*(u2*u2+v2*v2);
    p2 = gam1*(rE2-r2*q2);
    %c2 = sqrt(gam*p2/r2);
    un2 = u2*nl(1)+v2*nl(2);    
else
    r2 = udg2(1);
    ru2 = udg2(2);
    rv2 = udg2(3);
    rw2 = udg2(4);
    rE2 = udg2(5);    
    u2 = ru2/r2;
    v2 = rv2/r2;
    w2 = rw2/r2;
    q2 = 0.5*(u2*u2+v2*v2+w2*w2);
    p2 = gam1*(rE2-r2*q2);
    %c2 = sqrt(gam*p2/r2);
    un2 = u2*nl(1)+v2*nl(2)+w2*nl(3);        
end
end

if (stage==1)
    po = (alpha*pout + dd(1,1)*p0 + dd(1,1)*r0*c0*(un-un0))/(dd(1,1)+alpha);
elseif (stage==2)
    po = (alpha*pout + dd(2,2)*p0 - d(2,1)*(p1-p0) + dd(2,2)*r0*c0*(un-un0) + dd(2,1)*r0*c0*(un1-un0))/(dd(2,2)+alpha);
elseif (stage==3)
    po = (alpha*pout + dd(3,3)*p0 - d(3,2)*(p2-p0) - d(3,1)*(p1-p0) + dd(3,3)*r0*c0*(un-un0) + dd(3,2)*r0*c0*(un2-un0) + dd(3,1)*r0*c0*(un1-un0))/(dd(2,2)+alpha);
end

end

function po = outletpressure2d1(udg, udg0, nl, dd, pout, alpha, gam) 

r0 = udg0(1);
ru0 = udg0(2);
rv0 = udg0(3);
rE0 = udg0(4);
u0 = ru0/r0;
v0 = rv0/r0;
q0 = 0.5*(u0*u0+v0*v0);
p0 = gam1*(rE0-r0*q0);
c0 = sqrt(gam*p0/r0);
un0 = u0*nl(1)+v0*nl(2);

r = udg(1);
ru = udg(2);
rv = udg(3);
u = ru/r;
v = rv/r;
un = u*nl(1)+v*nl(2);    
    
po = (alpha*pout + dd(1,1)*p0 + dd(1,1)*r0*c0*(un-un0))/(dd(1,1)+alpha);

end


function po = outletpressure2d2(udg, udg0, udg1, nl, dd, pout, alpha, gam) 

r0 = udg0(1);
ru0 = udg0(2);
rv0 = udg0(3);
rE0 = udg0(4);
u0 = ru0/r0;
v0 = rv0/r0;
q0 = 0.5*(u0*u0+v0*v0);
p0 = gam1*(rE0-r0*q0);
c0 = sqrt(gam*p0/r0);
un0 = u0*nl(1)+v0*nl(2);

r1 = udg1(1);
ru1 = udg1(2);
rv1 = udg1(3);
rE1 = udg1(4);
u1 = ru1/r1;
v1 = rv1/r1;
q1 = 0.5*(u1*u1+v1*v1);
p1 = gam1*(rE1-r1*q1);
un1 = u1*nl(1)+v1*nl(2);    

r = udg(1);
ru = udg(2);
rv = udg(3);
u = ru/r;
v = rv/r;
un = u*nl(1)+v*nl(2);    
    
po = (alpha*pout + dd(2,2)*p0 - d(2,1)*(p1-p0) + dd(2,2)*r0*c0*(un-un0) + dd(2,1)*r0*c0*(un1-un0))/(dd(2,2)+alpha);

end

function po = outletpressure2d3(udg, udg0, udg1, udg2, nl, dd, pout, alpha, gam) 

r0 = udg0(1);
ru0 = udg0(2);
rv0 = udg0(3);
rE0 = udg0(4);
u0 = ru0/r0;
v0 = rv0/r0;
q0 = 0.5*(u0*u0+v0*v0);
p0 = gam1*(rE0-r0*q0);
c0 = sqrt(gam*p0/r0);
un0 = u0*nl(1)+v0*nl(2);

r1 = udg1(1);
ru1 = udg1(2);
rv1 = udg1(3);
rE1 = udg1(4);
u1 = ru1/r1;
v1 = rv1/r1;
q1 = 0.5*(u1*u1+v1*v1);
p1 = gam1*(rE1-r1*q1);
un1 = u1*nl(1)+v1*nl(2);    

r2 = udg2(1);
ru2 = udg2(2);
rv2 = udg2(3);
rE2 = udg2(4);
u2 = ru2/r2;
v2 = rv2/r2;
q2 = 0.5*(u2*u2+v2*v2);
p2 = gam1*(rE2-r2*q2);
%c2 = sqrt(gam*p2/r2);
un2 = u2*nl(1)+v2*nl(2);    

r = udg(1);
ru = udg(2);
rv = udg(3);
u = ru/r;
v = rv/r;
un = u*nl(1)+v*nl(2);    
    
po = (alpha*pout + dd(3,3)*p0 - d(3,2)*(p2-p0) - d(3,1)*(p1-p0) + dd(3,3)*r0*c0*(un-un0) + dd(3,2)*r0*c0*(un2-un0) + dd(3,1)*r0*c0*(un1-un0))/(dd(2,2)+alpha);

end

function po = outletpressure3d1(udg, udg0, nl, dd, pout, alpha, gam) 

r0 = udg0(1);
ru0 = udg0(2);
rv0 = udg0(3);
rw0 = udg0(4);
rE0 = udg0(5);    
u0 = ru0/r0;
v0 = rv0/r0;
w0 = rw0/r0;
q0 = 0.5*(u0*u0+v0*v0+w0*w0);
p0 = gam1*(rE0-r0*q0);
c0 = sqrt(gam*p0/r0);
un0 = u0*nl(1)+v0*nl(2)+w0*nl(3);    

r = udg(1);
ru = udg(2);
rv = udg(3);
rw = udg(4);
u = ru/r;
v = rv/r;
w = rw/r;
un = u*nl(1)+v*nl(2)+w*nl(3);        
    
po = (alpha*pout + dd(1,1)*p0 + dd(1,1)*r0*c0*(un-un0))/(dd(1,1)+alpha);

end

function po = outletpressure3d2(udg, udg0, udg1, nl, dd, pout, alpha, gam) 

r0 = udg0(1);
ru0 = udg0(2);
rv0 = udg0(3);
rw0 = udg0(4);
rE0 = udg0(5);    
u0 = ru0/r0;
v0 = rv0/r0;
w0 = rw0/r0;
q0 = 0.5*(u0*u0+v0*v0+w0*w0);
p0 = gam1*(rE0-r0*q0);
c0 = sqrt(gam*p0/r0);
un0 = u0*nl(1)+v0*nl(2)+w0*nl(3);    

r1 = udg1(1);
ru1 = udg1(2);
rv1 = udg1(3);
rw1 = udg1(4);
rE1 = udg1(5);    
u1 = ru1/r1;
v1 = rv1/r1;
w1 = rw1/r1;
q1 = 0.5*(u1*u1+v1*v1+w1*w1);
p1 = gam1*(rE1-r1*q1);
%c1 = sqrt(gam*p1/r1);
un1 = u1*nl(1)+v1*nl(2)+w1*nl(3);        

r = udg(1);
ru = udg(2);
rv = udg(3);
rw = udg(4);
u = ru/r;
v = rv/r;
w = rw/r;
un = u*nl(1)+v*nl(2)+w*nl(3);        

po = (alpha*pout + dd(2,2)*p0 - d(2,1)*(p1-p0) + dd(2,2)*r0*c0*(un-un0) + dd(2,1)*r0*c0*(un1-un0))/(dd(2,2)+alpha);

end

function po = outletpressure3d3(udg, udg0, udg1, udg2, nl, dd, pout, alpha, gam) 

r0 = udg0(1);
ru0 = udg0(2);
rv0 = udg0(3);
rw0 = udg0(4);
rE0 = udg0(5);    
u0 = ru0/r0;
v0 = rv0/r0;
w0 = rw0/r0;
q0 = 0.5*(u0*u0+v0*v0+w0*w0);
p0 = gam1*(rE0-r0*q0);
c0 = sqrt(gam*p0/r0);
un0 = u0*nl(1)+v0*nl(2)+w0*nl(3);    

r1 = udg1(1);
ru1 = udg1(2);
rv1 = udg1(3);
rw1 = udg1(4);
rE1 = udg1(5);    
u1 = ru1/r1;
v1 = rv1/r1;
w1 = rw1/r1;
q1 = 0.5*(u1*u1+v1*v1+w1*w1);
p1 = gam1*(rE1-r1*q1);
%c1 = sqrt(gam*p1/r1);
un1 = u1*nl(1)+v1*nl(2)+w1*nl(3);        

r2 = udg2(1);
ru2 = udg2(2);
rv2 = udg2(3);
rw2 = udg2(4);
rE2 = udg2(5);    
u2 = ru2/r2;
v2 = rv2/r2;
w2 = rw2/r2;
q2 = 0.5*(u2*u2+v2*v2+w2*w2);
p2 = gam1*(rE2-r2*q2);
%c2 = sqrt(gam*p2/r2);
un2 = u2*nl(1)+v2*nl(2)+w2*nl(3);        

r = udg(1);
ru = udg(2);
rv = udg(3);
rw = udg(4);
u = ru/r;
v = rv/r;
w = rw/r;
un = u*nl(1)+v*nl(2)+w*nl(3);        

po = (alpha*pout + dd(3,3)*p0 - d(3,2)*(p2-p0) - d(3,1)*(p1-p0) + dd(3,3)*r0*c0*(un-un0) + dd(3,2)*r0*c0*(un2-un0) + dd(3,1)*r0*c0*(un1-un0))/(dd(2,2)+alpha);

end
