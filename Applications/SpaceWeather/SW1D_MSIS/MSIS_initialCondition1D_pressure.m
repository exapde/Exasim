function u = MSIS_initialCondition1D_pressure(xdg,mu,indices,mass)

%number of dimensions/components
nc = 3;
nd = 1;

%parameters
r0 = mu(1);
lat0 = mu(2);
long = wrapTo180(mu(3));
year = mu(4);
doy = mu(5);
sec = mu(6);
F10p7 = mu(7);
F10p7a = mu(8);

hbot = mu(9);
H = mu(10);
T0 = mu(11);
rho0 = mu(12);

Fr = mu(13);
m = mu(14);

%computation
h = (xdg-r0)*H + hbot;
npoints = length(h);

lat = lat0*ones(npoints,1);
long = long*ones(npoints,1);
year = year*ones(npoints,1);
doy = doy*ones(npoints,1);
sec = sec*ones(npoints,1);
F10p7 = F10p7*ones(npoints,1);
F10p7a = F10p7a*ones(npoints,1);
LST = long/15 + sec/86400;

aph = zeros(npoints,7);
aph(:,1) = 4;
flags = ones(23,1);
flags(9) = -1;

dr = 500;
[TAll,rhoAll] = atmosnrlmsise00(h,lat,long,year,doy,sec,LST,F10p7a,F10p7,aph,flags);
[TAm,rhoAm] = atmosnrlmsise00(h-dr,lat,long,year,doy,sec,LST,F10p7a,F10p7,aph,flags);
[TAp,rhoAp] = atmosnrlmsise00(h+dr,lat,long,year,doy,sec,LST,F10p7a,F10p7,aph,flags);

rho = rhoAll(:,indices)*mass*m/rho0;
rhom = rhoAm(:,indices)*mass*m/rho0;
rhop = rhoAp(:,indices)*mass*m/rho0;


mass0 = (rhoAll(:,indices).*mass')./(rhoAll(:,indices)*mass)*mass;
massm = (rhoAm(:,indices).*mass')./(rhoAm(:,indices)*mass)*mass;
massp = (rhoAp(:,indices).*mass')./(rhoAp(:,indices)*mass)*mass;

% rho = rho/rho(1);
% rhop = rhop/rho(1);
% rhom = rhom/rho(1);

T = TAll(:,2)/T0;
Tm = TAm(:,2)/T0;
Tp = TAp(:,2)/T0;

rT = rho.*T./mass0;
rTm = rhom.*Tm./massm;
rTp = rhop.*Tp./massp;

drT = rTp - rTm;
drTdr = H*drT/(2*dr);

% acc = (Fr^2*xdg - (r0./xdg).^2);
acc = (Fr^2*xdg*cos(lat0)^2 - (r0./xdg).^2);

rho = drTdr./acc;
% rho = rho/rho(1);

T = mass0.*rT./rho;
T = T - T(1) + 1;
rho = rho/rho(1);
r = log(rho);
srT = sqrt(rho).*T;

drdx = gradient(r)./gradient(xdg);
dsrTdx = gradient(srT)./gradient(xdg);

u = zeros(npoints,nc*(nd+1));
iu = [1,nc,nc+1,nc*(nd+1)];
u(:,iu) = [r,srT,drdx,dsrTdx];

