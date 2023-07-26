function u = MSIS_initialCondition1D(xdg,mu,indices,mass)

%number of dimensions/components
nc = 3;
nd = 1;

%parameters
r0 = mu(1);
lat = mu(2);
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
m = mu(14);

%computation
h = (xdg-r0)*H + hbot;
npoints = length(h);

lat = lat*ones(npoints,1);
long = long*ones(npoints,1);
year = year*ones(npoints,1);
doy = doy*ones(npoints,1);
sec = sec*ones(npoints,1);
F10p7 = F10p7*ones(npoints,1);
F10p7a = F10p7a*ones(npoints,1);
LST = long/15;

aph = zeros(npoints,7);
aph(:,1) = 4;
flags = ones(23,1);
flags(9) = -1;

[TAll,rhoAll] = atmosnrlmsise00(h,lat,long,year,doy,sec,LST,F10p7a,F10p7,aph,flags);
rho = rhoAll(:,indices)*mass*m/rho0;
r = log(rho) - log(rho(1));
srT = sqrt(exp(r)).*(TAll(:,2)-TAll(1,2)+T0)/T0;

drdx = gradient(r)./gradient(xdg);
dsrTdx = gradient(srT)./gradient(xdg);

u = zeros(npoints,nc*(nd+1));
iu = [1,nc,nc+1,nc*(nd+1)];
u(:,iu) = [r,srT,drdx,dsrTdx];

