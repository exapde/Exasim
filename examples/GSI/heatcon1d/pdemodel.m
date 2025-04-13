
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

%C  = u(1);
T  = u(2);
Cx = q(1);
Tx = q(2);
% Cy = q(3);
% Ty = q(4);

% mu = [L Twall T0 R kB k H0 mu0 C0 D0 H tc Dc DTHf CvHf];

% L = mu(1);     % physical length
% Twall = mu(2); % wall temperature
T0 = mu(3);
kB = mu(5);
H = mu(11);
Dc = mu(13);
DTHf = mu(14);

Tphysical = T * T0;
DC = Dc*exp(-H/(kB*Tphysical));% dimensionless diffusion coefficient 
DT = DTHf;

f = [DC*Cx DT*Tx];
f = reshape(f, [2 1]);

%f = [DC*Cx DT*Tx DC*Cy DT*Ty];
%f = reshape(f, [2 2]);

end

function s = source(u, q, w, v, x, t, mu, eta)

C  = u(1);
T = u(2);

% mu = [L Twall T0 R kB k H0 mu0 C0 D0 H tc Dc DTHf CvHf];

% L = mu(1);   % physical length
% Twall = mu(2); % wall temperature
T0 = mu(3);
R = mu(4);
k = mu(6);  
H0 = mu(7);
mu0 = mu(8);
C0 = mu(9);
tc = mu(12);  
CvHf = mu(15);

KC = k*tc;   % dimensionless reaction rate

k1 = k*H0/CvHf; 
Tphysical = T * T0;
mup = mu0 + R * Tphysical * log(C); % J / mol
k2 = C0*k*mup/CvHf; % K / s

KT = tc*(k1 + k2)/T0; % dimensionless rate

s = [-KC*C; KT*C];
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
fb = [sym(0.0); sym(0.0)];
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
ub = [sym(0.0); sym(0.0)];
end

function u0 = initu(x, mu, eta)

Twall = mu(2); % wall temperature
T0 = mu(3);
u0 = [exp(-1.0e4*x*x); (Twall/T0)*sym(1.0)];
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)

C  = u(1);
T  = u(2); 
f = flux(u, q, w, v, x, t, mu, eta);

Twall = mu(2);
T0 = mu(3);
Cb = 1; 

fb1 = [tau*(Cb - uhat(1));  f(2)*n(1)+tau*(T-uhat(2))];
fb2 = [f(1)+tau*(C - uhat(1)); tau*(Twall/T0 - uhat(2))];

fb = [fb1(:) fb2(:)];
end



