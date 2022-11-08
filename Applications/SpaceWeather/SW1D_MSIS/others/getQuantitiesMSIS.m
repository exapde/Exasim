clear all
close all

tic
angleres = 5;
h = (100:5:600)*1e3;
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

year = 2015*ones(nsamplesTot,1);
doy = 1*ones(nsamplesTot,1);
sec = 86400/2*ones(nsamplesTot,1);

F10p7 = 150*ones(nsamplesTot,1);
F10p7a = 150*ones(nsamplesTot,1);
aph = zeros(nsamplesTot,7);
aph(:,1) = 4;
flags = ones(23,1);
flags(9) = -1;

[TAll,rhoAll] = atmosnrlmsise00(H,Lat,Long,year,doy,sec,LST,F10p7a,F10p7,aph,flags);

Th = zeros(nH,1);
rhoh = zeros(nH,1);
nOh = zeros(nH,1);
nO2h = zeros(nH,1);
nN2h = zeros(nH,1);
nHeh = zeros(nH,1);
nArh = zeros(nH,1);
nHh = zeros(nH,1);
nNh = zeros(nH,1);

for iH=1:nH
    Ti = TAll(1+nsamplesH*(iH-1):nsamplesH*iH,:);
    rhoi = rhoAll(1+nsamplesH*(iH-1):nsamplesH*iH,:);
    
    Th(iH) = mean(Ti(:,2));
    rhoh(iH) = mean(rhoi(:,6));
    nOh(iH) = mean(rhoi(:,2));
    nO2h(iH) = mean(rhoi(:,4));
    nN2h(iH) = mean(rhoi(:,3));
    
    nHeh(iH) = mean(rhoi(:,1));
    nArh(iH) = mean(rhoi(:,5));
    nHh(iH) = mean(rhoi(:,7));
    nNh(iH) = mean(rhoi(:,8));
end

nh = nOh + nO2h + nN2h + nHeh;
chiOh = nOh./nh;
chiO2h = nO2h./nh;
chiN2h = nN2h./nh;
chiHeh = nHeh./nh;

nh2 = nh + nArh + nHh + nNh;
chiOh = nOh./nh2;
chiO2h = nO2h./nh2;
chiN2h = nN2h./nh2;
chiHeh = nHeh./nh2;
chiArh = nArh./nh2;
chiHh = nHh./nh2;
chiNh = nNh./nh2;


rhoh2 = 16*nOh + 32*nO2h + 28*nN2h + 4*nHeh;
rhoiOh = 16*nOh./rhoh2;
rhoiO2h = 32*nO2h./rhoh2;
rhoiN2h = 28*nN2h./rhoh2;
rhoiHeh = 4*nHeh./rhoh2;

fO2 = fit(h'-100e3,rhoiO2h,'exp2');
fN2 = fit(h'-100e3,rhoiN2h,'exp2');
fHe = fit(h'-100e3,rhoiHeh,'exp2');
toc
fT = fit(h'-100e3,Th,'exp2');

figure
plot(rhoiN2h,h/1000,'k',rhoiO2h,h/1000,'b',rhoiOh,h/1000,'r',rhoiHeh,h/1000,'m','LineWidth',2)
hold on
plot(fN2(h-100e3),h/1000,'-.k',fO2(h-100e3),h/1000,'-.b',1-fN2(h-100e3)-fO2(h-100e3)-fHe(h-100e3),h/1000,'-.r',fHe(h-100e3),h/1000,'-.m','LineWidth',2)

set(gca,'FontSize',16,'TickLabelInterpreter','latex')
leg = legend('N2','O2','O','He','location','southeast');
xlabel('$\chi_i \, (\rho_i/\rho)$','Interpreter','latex','FontName','mwa_cmr10','FontSize',16+4)
ylabel('$h$ (km)','Interpreter','latex','FontName','mwa_cmr10','FontSize',16+4)
set(leg,'Interpreter','latex','Fontsize',16);
grid on
grid minor
set(gcf,'color','w');

mO = 16*rhoiOh;
mO2 = 32*rhoiO2h;
mN2 = 28*rhoiN2h;
mHe = 4*rhoiHeh;
m = (mO + mO2 + mN2 + mHe)/16;

figure
plot(m,h/1000,'b','LineWidth',2)
set(gca,'FontSize',16,'TickLabelInterpreter','latex')
xlabel('$m/m_{O}$','Interpreter','latex','FontName','mwa_cmr10','FontSize',16+4)
ylabel('$h$ (km)','Interpreter','latex','FontName','mwa_cmr10','FontSize',16+4)
grid on
grid minor
set(gcf,'color','w');

% figure
% plot(h,chiN2h,'k',h,chiO2h,'b',h,chiOh,'r',h,chiHeh,'m',h,chiArh,'g',h,chiHh,'c',h,chiNh,'y','LineWidth',2)
% legend('N2','O2','O','He','Ar','H','N')

% figure
% plot(h,chiN2h,'k',h,chiO2h,'b',h,chiOh,'r',h,chiHeh,'m',h,chiArh,'g',h,chiHh,'c',h,chiNh,'y','LineWidth',2)
% legend('N2','O2','O','He','Ar','H','N')