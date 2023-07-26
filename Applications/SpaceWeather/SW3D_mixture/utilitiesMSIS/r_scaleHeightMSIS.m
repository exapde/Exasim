function xH = r_scaleHeightMSIS(mu,indices,mass)

%parameters
r0 = mu(1);
r1 = mu(2);
year0 = mu(3);
doy0 = mu(4);
F10p70 = mu(6);
F10p7a0 = mu(7);

hbot = mu(8);
H = mu(9);
T0 = mu(10);

Fr = mu(12);

%computation
h0 = logspace(-2,log10(r1-r0),10000);
h = h0*H + hbot;
lat = 0;
long = 0;
npoints = length(h);

year = year0*ones(npoints,1);
doy = doy0*ones(npoints,1);
sec = 43200*ones(npoints,1);
F10p7 = F10p70*ones(npoints,1);
F10p7a = F10p7a0*ones(npoints,1);
LST = long/15 + sec/86400;

aph = zeros(npoints,7);
aph(:,1) = 4;
flags = ones(23,1);
flags(9) = -1;

[TAll,rhoAll] = atmosnrlmsise00(h,lat,long,year,doy,sec,LST,F10p7a,F10p7,aph,flags);
mass0 = (rhoAll(:,indices).*mass')./(rhoAll(:,indices)*mass)*mass;
T = TAll(:,2)/T0;

acc = (Fr^2*(h0+r0) - (r0./(h0+r0)).^2);
dp_p = acc'.*mass0./T;
dr = h0(2:end)-h0(1:end-1);
dI = 0.5*(dp_p(2:end)+dp_p(1:end-1)).*dr';

logp_p0 = [0;cumsum(dI)];

drop = 2;
div = 1;
nr = ceil(-logp_p0(end)/log(drop^(1/div)))+1;
r = exp(-logp_p0(end)/(nr-1));
xH = zeros(nr,1);
rate = log(r.^(-(1:1:nr-1)));
pp0 = logp_p0;

for ix=1:nr-1
    ri = rate(ix);
    [~,ind] = min(abs(pp0-ri));
    f1 = pp0(ind)-ri;
    if f1>1e-4
        ind2 = ind+1;  
    else
        ind2 = ind-1;
    end
    f2 = pp0(ind2)-ri;
    h1 = h0(ind);
    h2 = h0(ind2);

    hi = h1 + f1*(h2-h1)/(f1-f2);
    xH(ix+1) = hi/(r1-r0);
end
