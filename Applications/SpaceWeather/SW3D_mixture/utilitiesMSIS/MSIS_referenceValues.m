function [rho0,T0,chi,cchi] = MSIS_referenceValues(indices,mass,H0,H1,year,doy,sec,F10p7,F10p7a,Ap)

nspecies = length(indices);

amu = 1.66e-27;
angleres = 5;
h = (H0:5000:H1);
lat = -90+angleres:angleres:90-angleres;
long = -180:angleres:180-angleres;

nH = length(h);
nlat = length(lat);
nlong = length(long);
nsamplesH = nlat*nlong;
nsamplesTot = nsamplesH*nH;

Lat = repmat(reshape(ones(nlong,1)*lat,nsamplesH,1),[nH,1]);
Long = repmat(repmat(long',[nlat,1]),[nH,1]);
LST = Long/15;
H = reshape(ones(nsamplesH,1)*h,nsamplesTot,1);

year = year*ones(nsamplesTot,1);
doy = doy*ones(nsamplesTot,1);
sec = sec*ones(nsamplesTot,1);

F10p7 = F10p7*ones(nsamplesTot,1);
F10p7a = F10p7a*ones(nsamplesTot,1);
aph = ones(nsamplesTot,1)*Ap;

[TAll,rhoAll] = atmosnrlmsise00(H,Lat,Long,year,doy,sec,LST,F10p7a,F10p7,aph);

Th = zeros(nH,1);
nh = zeros(nH,nspecies);

for iH=1:nH
    Ti = TAll(1+nsamplesH*(iH-1):nsamplesH*iH,:);
    rhoi = rhoAll(1+nsamplesH*(iH-1):nsamplesH*iH,:);
    
    Th(iH) = mean(Ti(:,2));
    nh(iH,:) = mean(rhoi(:,indices),1);
end

rhoh = nh*mass;
chi = (nh.*mass')./rhoh;

rho0 = rhoh(1);
T0 = Th(1);

cchi = zeros(nspecies-1,4);
options = fitoptions('exp2');

mindensity = 1;
densityO = 1;
for ispe = 1:nspecies-1
    cchi(ispe,:) = coeffvalues(fit(h'-H0,chi(:,ispe+1),'exp2',options));
    density_i = cchi(ispe,1)*exp(cchi(ispe,2)*(h'-H0)) + cchi(ispe,3)*exp(cchi(ispe,4)*(h'-H0));
    densityO = densityO - density_i;
    mindensity = min(mindensity,min(density_i));
end
mindensity = min(mindensity,min(densityO));

if min(mindensity) < 0
    mindensity = 1;
    densityO = 1;
    options.Robust = 'LAR';
    for ispe = 1:nspecies-1
        cchi(ispe,:) = coeffvalues(fit(h'-H0,chi(:,ispe+1),'exp2',options));
        density_i = cchi(ispe,1)*exp(cchi(ispe,2)*(h'-H0)) + cchi(ispe,3)*exp(cchi(ispe,4)*(h'-H0));
        densityO = densityO - density_i;
        mindensity = min(mindensity,min(density_i));
    end
    mindensity = min(mindensity,min(densityO));
end

if min(mindensity) < 0
    error('Negative density detected in exponential fit'); 
end


