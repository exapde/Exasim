function u = MSIS_initialCondition3D(xdg,mu,indices,mass)

%number of dimensions/components
nc = 5;
nd = 3;

%parameters
r0 = mu(1);
year = mu(3);
doy = mu(4);
sec = mu(5);
F10p7 = mu(6);
F10p7a = mu(7);

hbot = mu(8);
H = mu(9);
T0 = mu(10);
rho0 = mu(11);

m = mu(13);

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

h0 = hbot*ones(npoints,1);

[TAll,rhoAll] = atmosnrlmsise00(h,lat,long,year,doy,sec,LST,F10p7a,F10p7,aph,flags);
[T0all,rho0all] = atmosnrlmsise00(h0,lat,long,year,doy,sec,LST,F10p7a,F10p7,aph,flags);

rhoh = rhoAll(:,indices)*mass*m;
rhoh0 = rho0all(:,indices)*mass*m;

r = log(rhoh/rho0) - log(rhoh0/rho0);
rho = exp(r);
T = (TAll(:,2)-T0all(:,2)+T0)/T0;

srT = sqrt(rho).*T;

drdx = gradient(r)./gradient(x1);
drdy = gradient(r)./gradient(x2);
drdz = gradient(r)./gradient(x3);

dsrTdx = gradient(srT)./gradient(x1);
dsrTdy = gradient(srT)./gradient(x2);
dsrTdz = gradient(srT)./gradient(x3);

u = zeros(npoints,nc*(nd+1));
iu = sort([(nc*(0:nd)+1) nc*(1:(nd+1))]);
u(:,iu) = [r,srT,0*drdx,0*dsrTdx,0*drdy,0*dsrTdy,0*drdz,0*dsrTdz];

