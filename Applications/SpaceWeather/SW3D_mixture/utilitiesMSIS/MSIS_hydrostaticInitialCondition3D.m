function u = MSIS_hydrostaticInitialCondition3D(xdg,mu,indices,mass)

%number of dimensions/components
nc = 5;
nd = 3;

%parameters
r0 = mu(1);
r1 = mu(2);
year0 = mu(3);
doy0 = mu(4);
sec0 = mu(5);
F10p70 = mu(6);
F10p7a0 = mu(7);

hbot = mu(8);
H = mu(9);
T0 = mu(10);
rho0 = mu(11);

Fr = mu(12);
m = mu(13);

Ap = mu(14:20);

%Change of coordinates
x1 = xdg(:,1);
x2 = xdg(:,2);
x3 = xdg(:,3);

radius = sqrt(x1.^2 + x2.^2 + x3.^2);      %radial position  

sgnx1 = tanh(100*x1);
sgnx2 = tanh(100*x2);
sgnx2pow2 = tanh(100*x2.*x2);
sgnx1x2 = tanh(100*x1.*x2);
ax1x2 = abs(x1) + abs(x2);
sA = tanh(100*(ax1x2 - 1e-16));
ax1x2p = ax1x2.*sA + 1e-16*(1-sA);
long = ((pi - pi/2*(1+sgnx1).*(1-sgnx2pow2) - pi/4*(2+sgnx1).*sgnx2 - sgnx1x2.*atan((abs(x1) - abs(x2))./(ax1x2p))));
lat = asin(x3./radius);

%computation
h = (radius-r0)*H + hbot;
lat = lat*180/pi;
long = wrapToPi(long)*180/pi;
npoints = length(h);

year = year0*ones(npoints,1);
doy = doy0*ones(npoints,1);
sec = sec0*ones(npoints,1);
F10p7 = F10p70*ones(npoints,1);
F10p7a = F10p7a0*ones(npoints,1);
LST = long/15 + sec/86400;

aph = ones(npoints,1)*Ap;

h0 = hbot*ones(npoints,1);
dr = 500;
[TAll,rhoAll] = atmosnrlmsise00(h,lat,long,year,doy,sec,LST,F10p7a,F10p7,aph);
[TAm,rhoAm] = atmosnrlmsise00(h-dr,lat,long,year,doy,sec,LST,F10p7a,F10p7,aph);
[TAp,rhoAp] = atmosnrlmsise00(h+dr,lat,long,year,doy,sec,LST,F10p7a,F10p7,aph);
[T0all,rho0all] = atmosnrlmsise00(h0,lat,long,year,doy,sec,LST,F10p7a,F10p7,aph);

rho = rhoAll(:,indices)*mass*m/rho0;
rhom = rhoAm(:,indices)*mass*m/rho0;
rhop = rhoAp(:,indices)*mass*m/rho0;
rhoh0 = rho0all(:,indices)*mass*m/rho0;

mass0 = (rhoAll(:,indices)*mass)./sum(rhoAll(:,indices),2);
massm = (rhoAm(:,indices)*mass)./sum(rhoAm(:,indices),2);
massp = (rhoAp(:,indices)*mass)./sum(rhoAp(:,indices),2);

T = TAll(:,2)/T0;
Tm = TAm(:,2)/T0;
Tp = TAp(:,2)/T0;
Th0 = T0all(:,2)/T0;

p = rho.*T./mass0;
pm = rhom.*Tm./massm;
pp = rhop.*Tp./massp;

dp = pp - pm;
dpdr = H*dp/(2*dr);

acc = (Fr^2*radius - (r0./radius).^2);
rho = dpdr./acc;

T = mass0.*p./rho;
T = T - Th0 + 1;
rho = rho./rhoh0;
r = log(rho);
srT = sqrt(rho).*T;

drdx = gradient(r)./gradient(x1);
drdy = gradient(r)./gradient(x2);
drdz = gradient(r)./gradient(x3);

dsrTdx = gradient(srT)./gradient(x1);
dsrTdy = gradient(srT)./gradient(x2);
dsrTdz = gradient(srT)./gradient(x3);

u = zeros(npoints,nc*(nd+1));
iu = sort([(nc*(0:nd)+1) nc*(1:(nd+1))]);
u(:,iu) = [r,srT,drdx,dsrTdx,drdy,dsrTdy,drdz,dsrTdz];

