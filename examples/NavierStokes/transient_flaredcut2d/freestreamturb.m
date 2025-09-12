

% ind = find((mesh.f(:,end)==-2) | (mesh.f(:,end)==-4));
% elem = mesh.f(ind,end-1);
% n = length(elem);
% npf = size(mesh.plocfc,1);
% xdg = zeros(npf,mesh.nd,n);
% for i = 1:n
%     ei = mesh.t2f(elem(i),:);
%     for j = 1:4        
%         if (mesh.f(ei(j),end)==-2) || (mesh.f(ei(j),end)==-4)
%             xdg(:,:,i) = mesh.dgnodes(mesh.perm(:,j),:,elem(i));
%         end
%     end
% end
% xdg = reshape(permute(xdg, [1 3 2]),[npf*n mesh.nd]);
% xdg = unique(xdg,'rows');

%running with line case% 
X = linspace(0,1,1000);
Y = linspace(0,1,1000);
gam = 1.4;
Minf = 7;
xdg = [0*Y(:) Y(:) 0*Y(:)];

sdg = [0; sqrt(sum((xdg(2:end,:)-xdg(1:end-1,:)).^2,2))];
sdg = cumsum(sdg);

stgNmode = 200;
gridLength = 1e-03;
turbLengthFactor = 10;
visc = 3.3551679839273703e-06;%dynamic (I think pascals*s)
turbIntensity = 1/100;
Ustg = 1010.38;%m/s
stgdata = stghomogeneousturbulence(gridLength, turbLengthFactor, visc, turbIntensity, Ustg, stgNmode+1);

xdg(:,3) = 0;
rmean = 1;
umean = 1;
vmean = 0;
wmean = 0;
pmean = 1/(gam*Minf^2);
U = 2.5; 
dt = 5e-5;
for i = 1:100    
    [un, vn, wn] = homogenousisotropicconvectivefluctuations(xdg, stgdata, U, i*dt);
    
    rn = -2*un;
    pn = 0;
    
    r = rmean + rn;
    ru = (umean + un).*r;
    rv = (vmean + vn).*r;
    rw = (wmean + wn).*r;    
    rE = (pmean+pn)/(gam-1) + 0.5*(ru.*ru + rv.*rv + rw.*rw)./r;
    
    figure(6); clf;
    plot(sdg, r, sdg, ru, sdg, rmean*ones(size(sdg)));
    
%     figure(6); clf;
%     plot(sdg, rv, sdg, rw, sdg, r.*un, sdg, 0*un);
%     
%     figure(7); clf;
%     plot(sdg, rE, sdg, (pmean/(gam-1) + 0.5*(umean^2)*rmean)*ones(size(sdg)));
%     
    if i<=3
        pause
    end
    pause(0.1)    
end
return
n = 1e4;
dt = 4e-5;
a=zeros(n,3);
for i = 1:n
    [un, vn, wn] = homogenousisotropicconvectivefluctuations(xdg(800,:), stgdata, 1, i*dt);
    a(i,1) = (umean + un);
    [un, vn, wn] = homogenousisotropicconvectivefluctuations(xdg(800,:), stgdata, 2, i*dt);
    a(i,2) = (umean + un);
    [un, vn, wn] = homogenousisotropicconvectivefluctuations(xdg(800,:), stgdata, 2.5, i*dt);
    a(i,3) = (umean + un);
end

n=200;
figure(1); clf;
plot((1:n)*dt,a(1:n,1),(1:n)*dt,a(1:n,2),(1:n)*dt,a(1:n,3));

n = 2000;
b = 1e3*dt/Ustg;
figure(1); clf; plot((1:n)*b,a(1:n,1),(1:n)*b,a(1:n,2),(1:n)*b,a(1:n,3));


% gam = 1.4;
% Re = 10.8224e6;
% Pr = 0.71;
% Minf = 6.0;
% Tref  = 51.2195;
% Twall = 300;
% pinf = 1/(gam*Minf^2);
% Uinf = [ 1.0, 1.0, 0.0, 0.0, 0.5+pinf/(gam-1)];
% rmean = 1;
% umean = 1;
% vmean = 0;
% wmean = 0;
% Tmean = 1;
% U = 1;

% minLength = 0.5/400;
% %visc = 1e-6;
% visc = 7.955401172478808e-05*2;
% turbIntensity = 1/100;
% N = 50;
% dt = 1e-4/2;
% 
% %[fluctamp,pars,Ek,Em,awaveno,dwaveno,waveno] = homogeneousturbulence(minLength, visc, turbIntensity, N);
% %[fluctamp,pars,Ek,Em,awaveno,dwaveno,waveno] = homogeneousturbulence(minLength, visc, turbIntensity*861, N);
% [fluctamp,pars,Ek,Em,awaveno,dwaveno,waveno] = homogeneousturbulence(minLength, 20*minLength, visc, turbIntensity*861, N+1);
% fluctamp = fluctamp/861;
% 
% randno = randomgen(length(dwaveno));
% randno(:,8) = 0;
% for n=1:length(dwaveno)
%     c = turbIntensity*awaveno(n);
%     randno(n,8) = normrnd(c,c);
% end
% xdg(:,3) = 0;
% %randno(:,8)=0;
% 
% for i = 1:1000    
%     [un, vn, wn] = homogenousisotropicconvectivefluctuations(xdg, fluctamp, awaveno, randno, U, i*dt);
%     
%     %up = sqrt((un.*un + vn.*vn + wn.*wn)/3);
%     %Tp = -((gam-1)*Minf*Minf*umean)*un;
%     %rp = -(rmean/Tmean)*Tp;
%     rn = 0*un;
%     pn = 2*umean*un/(gam*Minf*Minf);
%     
%     r = rmean + rn;
%     ru = (umean + un).*r;
%     rv = (vmean + vn).*r;
%     rw = (wmean + wn).*r;
%     %T = Tmean + Tp;
%     %rE = r.*T./((gamma-1)*gam*Minf*Minf) + 0.5*(ru.*ru + rv.*rv + rw.*rw)./r;
%     rE = (pinf+pn)/(gam-1) + 0.5*(ru.*ru + rv.*rv + rw.*rw)./r;
%     
%     figure(2); clf;
%     plot(sdg, r, sdg, ru);
%     
%     figure(3); clf;
%     plot(sdg, rv, sdg, rw, sdg, r.*un, sdg, 0*un);
%     
%     figure(4); clf;
%     plot(sdg, rE, sdg, Uinf(end)*ones(size(un)));
%     
%     pause(3)
%     
% %     figure(3); clf;
% %     plot(sdg, un, '-o', sdg, vn, sdg, wn);
% %     i
% %     pause(0.1)
% end
% 
% 
