function [U,Vprime] = stgcomp(xdg, t, gam, Minf, Uinf, lemax, kcute, waveno, randno, Umean, Amean)

% velocity field
Vprime = fluctvel(xdg, t,  Uinf, lemax, kcute, waveno, randno, Amean);
u = Umean(:,2) + Vprime(:,1);
v = Umean(:,3) + Vprime(:,2);
w = Umean(:,4) + Vprime(:,3);

% density and temperature
gam1 = gam-1;
Tprime = - (gam1*Minf*Minf)*Umean(:,2).*Vprime(:,1);
Rprime = - Tprime .* (Umean(:,1)./Umean(:,5));
r = Umean(:,1) + Rprime;    
T = Umean(:,5) + Tprime;   

% conservative variables
U(:,1) = r;
U(:,2) = r .* u;
U(:,3) = r .* v;
U(:,4) = r .* w;
U(:,5) = r .* T/((gam-1)*gam*Minf*Minf) + 0.5 * r .* (u.*u + v.*v + w.*w);            






