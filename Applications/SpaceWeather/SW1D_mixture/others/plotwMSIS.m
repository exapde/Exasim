close all
clear all

colours = ['b','g','k','r'];
nothing = NaN;
fontSize = 16;
hcomp = [150 200 250 400]'*1e3; %hm

% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% generate input files and store them in datain folder
load('../inputFiles');
indices = [2;3;4;1];

r0 = pde.physicsparam(16);
r1 = pde.physicsparam(17);
H = pde.physicsparam(18);
T0 = pde.physicsparam(19);
rho0 = pde.physicsparam(20);
t0 = pde.physicsparam(21);

lat = pde.physicsparam(23);
long = wrapTo180(pde.physicsparam(24));
F10p7 = pde.physicsparam(9);
F10p7a = pde.physicsparam(10);

Nday = pde.physicsparam(11);
year0 = pde.physicsparam(26);
doy0 = floor(Nday);
sec0 = (Nday - doy0)*86400;

ndaysyear = 365+(mod(year0,4)<1);


pde.soltime = pde.soltime + pde.timestepOffset;
% pde.soltime = 20:20:3500;
% ttotal = cumsum(pde.dt);
time = pde.dt(1)*pde.soltime*t0;
sec = mod(sec0+time,86400);
days1 = floor((sec0+time)/86400)+doy0;
days = mod(days1,ndaysyear);
year = floor(days1/ndaysyear)+year0;

% get solution from output files in dataout folder
sol = fetchsolution(pde,master,dmd,'dataout');
res = fetchresidual(pde);

npoints = length(hcomp);
ntime = length(pde.soltime);

hnond = r0 + (hcomp-1e5)/H;
dgnodes = mesh.dgnodes;

[ielem,xlocal] = findx1D(hnond,dgnodes);
plocal = masternodes(pde.porder,pde.nd,pde.elemtype);

V = vander(plocal);
yi = V'\(xlocal.^(pde.porder:-1:0))';
ui = permute(sol(:,:,ielem,:),[1 3 2 4]);

vecSize = ones(ntime,1);
aph = zeros(ntime,7);
aph(:,1) = 4;
flags = ones(23,1);
flags(9) = -1;

rho = zeros(ntime,npoints);
T = zeros(ntime,npoints);
rhoMSIS = zeros(ntime,npoints);
TMSIS = zeros(ntime,npoints);
for ipoint = 1:npoints
    
    hi = hcomp(ipoint)*vecSize;
    [TMSISi,rhoMSISi] = atmosnrlmsise00(hi,lat*vecSize,long*vecSize,year,days,sec,sec/3600 + long/15,F10p7a*vecSize,F10p7*vecSize,aph,flags);
    TMSIS(:,ipoint) = TMSISi(:,2);
    rhoMSIS(:,ipoint) = rhoMSISi(:,6);
    
    ux = pagemtimes(yi(:,ipoint)',ui(:,ipoint,:,:));
    
    r = exp(ux(:,:,1,:));
    rho(:,ipoint) = reshape(r,[ntime,1]);
    T(:,ipoint) = T0*reshape(ux(:,:,3,:)./sqrt(r),[ntime,1]);
end

time = time - time(1);
timeh = time/86400;
ntime = length(timeh);

figure
plot(nothing,nothing,':k',nothing,nothing,'-k','LineWidth', 1.5);
hold on
plot(nothing,nothing,sprintf('-%s',colours(1)),nothing,nothing,sprintf('-%s',colours(2)),nothing,nothing,sprintf('-%s',colours(3)),nothing,nothing,sprintf('-%s',colours(4)),'LineWidth', 2);
plot(timeh,T(:,1),sprintf('-%s',colours(1)),timeh,T(:,2),sprintf('-%s',colours(2)),timeh,T(:,3),sprintf('-%s',colours(3)),timeh,T(:,4),sprintf('-%s',colours(4)),'LineWidth', 2);
plot(timeh,TMSIS(:,1),sprintf(':%s',colours(1)),timeh,TMSIS(:,2),sprintf(':%s',colours(2)),timeh,TMSIS(:,3),sprintf(':%s',colours(3)),timeh,TMSIS(:,4),sprintf(':%s',colours(4)),'LineWidth', 2);
set(gca,'FontSize',16,'TickLabelInterpreter','latex')
leg = legend('MSIS','Exasim','h=150km','h=200km','h=250km','h=400km','location','north','NumColumns',3);
xlabel('Time (days)','Interpreter','latex','FontName','mwa_cmr10','FontSize',16+4)
ylabel('$T$ (K)','Interpreter','latex','FontName','mwa_cmr10','FontSize',16+4)
set(leg,'Interpreter','latex','Fontsize',16);
axis([0 10 625 1425])
xticks(0:1:10)
grid on
grid minor
set(gcf,'color','w');


figure
semilogy(nothing,nothing,':k',nothing,nothing,'-k','LineWidth', 1.5);
hold on
semilogy(nothing,nothing,sprintf('-%s',colours(1)),nothing,nothing,sprintf('-%s',colours(2)),nothing,nothing,sprintf('-%s',colours(3)),nothing,nothing,sprintf('-%s',colours(4)),'LineWidth', 2);
semilogy(timeh,rho0*rho(:,1),sprintf('-%s',colours(1)),timeh,rho0*rho(:,2),sprintf('-%s',colours(2)),timeh,rho0*rho(:,3),sprintf('-%s',colours(3)),timeh,rho0*rho(:,4),sprintf('-%s',colours(4)),'LineWidth', 2);
semilogy(timeh,rhoMSIS(:,1),sprintf(':%s',colours(1)),timeh,rhoMSIS(:,2),sprintf(':%s',colours(2)),timeh,rhoMSIS(:,3),sprintf(':%s',colours(3)),timeh,rhoMSIS(:,4),sprintf(':%s',colours(4)),'LineWidth', 2);
set(gca,'FontSize',16,'TickLabelInterpreter','latex')
leg = legend('MSIS','Exasim','h=150km','h=200km','h=250km','h=400km','location','north','NumColumns',3);
xlabel('Time (days)','Interpreter','latex','FontName','mwa_cmr10','FontSize',16+4)
ylabel('$\rho$ (kg/m$^3$)','Interpreter','latex','FontName','mwa_cmr10','FontSize',16+4)
set(leg,'Interpreter','latex','Fontsize',16);
xticks(0:1:10)
grid on
grid minor
set(gcf,'color','w');



% figure
% % plot(time,T,'-','LineWidth',2)
% semilogy(timeh,rho0*rho,'-',timeh,rhoMSIS,':','LineWidth',2)

altitude = (dgnodes(:)-r0)*H + 1e5;
vecSize = ones(size(altitude));
aph = zeros(length(altitude),7);
aph(:,1) = 4;

% 
% for i=1:ntime
% % for i=1:1
%     rho = sol(:,1,:,i); 
%     vr = sol(:,2,:,i)./sqrt(exp(rho));
%     T = sol(:,3,:,i)./sqrt(exp(rho));
%     
%     [TMSISt,rhoMSISt] = atmosnrlmsise00(altitude,lat*vecSize,long*vecSize,year(i)*vecSize,days(i)*vecSize,sec(i)*vecSize,sec(i)/3600 + long*vecSize/15,F10p7a*vecSize,F10p7*vecSize,aph,flags);
%     TMSISt0 = TMSISt(:,2)/T0;
%     rhoMSISt0 = log(rhoMSISt(:,6)/rho0);
%     figure(3)
%     hold on
%     plot(dgnodes(:),T(:),'-b',dgnodes(:),rho(:),'-r','linewidth',2);
%     hold on
%     plot(dgnodes(:),TMSISt0,'-.b',dgnodes(:),rhoMSISt0,'-.r','linewidth',2);
%     axis([min(dgnodes(:)) max(dgnodes(:)) -18 10])
%     hold on
%     grid on
%     grid minor
%     
% %     plot(dgnodes(:),vr(:),'linewidth',2);
% %     plot(dgnodes(:),rho(:),'linewidth',2);
%     
%     text(min(dgnodes(:)) + 10,-11,sprintf('t = %d',i))
% 
%     pause(.1)
%     hold off
% 
% end
