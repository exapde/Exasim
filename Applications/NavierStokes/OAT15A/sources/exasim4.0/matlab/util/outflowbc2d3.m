function [ubou,uC,pbou3] = outflowbc2d3(udg, udg0, pbou0, pbou1, pbou2, nl, dirka, pout, alpha, gam) 

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
q = 0.5*(u*u+v*v);
p = gam1*(rE-r*q);
c = sqrt(gam*p/r);
un = u*nl(1)+v*nl(2);    
M = abs(un)/c;

% interior primitive variables
uL = udg;
uL(1) = r;
uL(2) = u;
uL(3) = v;
uL(4) = p;

% exterior primitive variables
% (pbou3-p0) - r0*c0*(un-un0) + a31*dt*alpha*(pbou1-pout) + a32*dt*alpha*(pbou2-pout) + a33*dt*alpha*(pbou3-pout) = 0
% (1+a33*dt*alpha)*pbou3 = p0 + r0*c0*(un-un0) - a31*dt*alpha*(pbou1-pout) - a32*dt*alpha*(pbou2-pout) + a33*dt*alpha*pout
pbou3 = (pbou0 + r0*c0*(un-un0) - dirka(3,1)*dt*alpha*(pbou1-pout) - dirka(3,2)*dt*alpha*(pbou2-pout) + dirka(3,3)*dt*alpha*pout)/(1+dirka(3,3)*dt*alpha);
uR = udg;
uR(1) = (pbou3/p)*r;
uR(2) = u;
uR(3) = v;
uR(4) = pbou3;

% blended primitive variables
a = 0.5*(1-tanh(50*(M-1.0)));
uC = a*uR + (1-a)*uL;

% conservative variables
ubou = udg;
ubou(1) = uC(1);
ubou(2) = uC(1)*uC(2);
ubou(3) = uC(1)*uC(3);
ubou(4) = uC(4)/gam1 + 0.5*uC(1)*(uC(2).^2+uC(3).^2);

end

