function u = MSIS_initialCondition3D_integral(xdg,mu,indices,mass)

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
htop = (r1-r0)*H + hbot;
lat = lat*180/pi;
long = wrapToPi(long)*180/pi;
npoints = length(h);

year = year0*ones(npoints,1);
doy = doy0*ones(npoints,1);
sec = sec0*ones(npoints,1);
F10p7 = F10p70*ones(npoints,1);
F10p7a = F10p7a0*ones(npoints,1);
LST = long/15;

aph = zeros(npoints,7);
aph(:,1) = 4;
flags = ones(23,1);
flags(9) = -1;

h0 = hbot*ones(npoints,1);
ndiv = 10;
[T0all,rho0all] = atmosnrlmsise00(h0,lat,long,year,doy,sec,LST,F10p7a,F10p7,aph,flags);
Tbase = T0all(:,2)/T0;
rhobase = rho0all(:,indices)*mass*m/rho0;

Tm = Tbase - Tbase + 1;
massm = (rho0all(:,indices).*mass')./(rho0all(:,indices)*mass)*mass;
aim = Fr^2*r0 - 1;
rim = r0;
logpm = aim.*massm./Tm;

logp_p0 = zeros(npoints,1);

tic
for iDiv = 1:ndiv
    hall = hbot + iDiv*(h-hbot)/ndiv;
    [Tiall,rhoiall] = atmosnrlmsise00(hall,lat,long,year,doy,sec,LST,F10p7a,F10p7,aph,flags);
    Tp = Tiall(:,2)/T0;
    massp = (rhoiall(:,indices).*mass')./(rhoiall(:,indices)*mass)*mass;
    rip = (hall-hbot)/H+r0;
    aip = Fr^2*rip - (r0./rip).^2;
    
    logpp = aip.*massp./Tp;
    
    dr = rip-rim;
    logp_p0 = logp_p0 + 0.5*dr.*(logpp+logpm);
    
    rim = rip;
    logpm = logpp;
end

r = logp_p0 - log(Tp) + log(massp);

T = Tp - Tbase + 1;
rho = exp(r)./rhobase;
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

