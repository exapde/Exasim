
function pde = pdemodel
pde.mass = @mass;
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.fbouhdg = @fbouhdg;
pde.ubou = @ubou;
pde.initu = @initu;
end

function m = mass(u, q, w, v, x, t, mu, eta)
m = sym([1.0; 1.0]); 
end

function f = flux(u, q, w, v, x, t, mu, eta)

T  = u(1);
Tx = q(1);
Cx = q(2);

L  = mu(1);
T0 = mu(3);
kB = mu(5);
D0 = mu(10);
H = mu(11);
tc = mu(12);
DT = mu(14);

Tphys = T*T0;
D = D0*exp(-H/(kB*Tphys));
DC = tc*D/(L*L);

f = [DT*Tx; DC*Cx];

end

function s = source(u, q, w, v, x, t, mu, eta)

T = u(1);
C = u(2);
Cl = lmax(C,1e3); % regularized concentration

% mu = [L Twall T0 R kB k H0 mu0 C0 D0 H tc Dc DTHf CvHf];

T0 = mu(3);
R = mu(4);
k = mu(6);  
H0 = mu(7);
mu0 = mu(8);
C0 = mu(9);
tc = mu(12);
CvHf = mu(15);

k1 = k*H0/CvHf; 
Tphysical = T * T0;
mup = mu0 + R * Tphysical * log(Cl); % J / mol
k2 = C0*k*mup/CvHf;   % K / s
KT = tc*(k1 + k2)/T0; % dimensionless rate
KC = k*tc;       % dimensionless reaction rate

s = [KT * Cl; -KC * Cl];
end


function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
fb = sym([0.0; 0.0]);
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
ub = sym([0.0; 0.0]);
end

function u0 = initu(x, mu, eta)

Twall = mu(2); % wall temperature
T0 = mu(3);
u0 = [(Twall/T0); exp(-1e4*x*x)];
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)

T  = u(1); 
C  = u(2);
f = flux(u, q, w, v, x, t, mu, eta);

Twall = mu(2);
T0 = mu(3);
heatflux = mu(16);

fb1 = [f(1)*n(1)+tau*(T-uhat(1))+heatflux tau*(1.0 - uhat(2))];
fb2 = [tau*(Twall/T0 - uhat(1))  f(2)*n(1)+tau*(C - uhat(2))];

fb = [fb1(:) fb2(:)];
end



