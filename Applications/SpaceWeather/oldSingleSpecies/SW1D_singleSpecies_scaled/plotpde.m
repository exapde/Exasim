% close all

% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% generate input files and store them in datain folder
load('inputFiles');
% tfixed = 172800;
% freq = 5;
% tinit = 0;
% 
pde.soltime = pde.soltime + pde.timestepOffset;

% get solution from output files in dataout folder
sol = fetchsolution(pde,master,dmd);
res = fetchresidual(pde);

dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);    
figure(1);clf;
for i=1:size(sol,4)
    rho = sol(:,1,:,i); 
    vr = sol(:,2,:,i)./sqrt(exp(rho));
    T = sol(:,3,:,i)./sqrt(exp(rho));
    
    plot(dgnodes(:),T(:),'linewidth',2);
    axis([min(dgnodes(:)) max(dgnodes(:)) -18 10])
    hold on
    grid on
    grid minor
    
    plot(dgnodes(:),vr(:),'linewidth',2);
    plot(dgnodes(:),rho(:),'linewidth',2);
    
    text(min(dgnodes(:)) + 10,-11,sprintf('t = %d',i))

    pause(.1)
    hold off

end
% 
% figure(2)
% mu = pde.physicsparam;
% x1 = dgnodes(:);
% gam = mu(1);
% Minf = mu(4);
% M2 = Minf^2;
% Fr2 = mu(5);
% omega = mu(6);
% 
% Tbot = mu(8);
% Ttop = mu(9);
% R0 = mu(10);
% Ldim = mu(14);
% h0 = 35000/Ldim;
% 
% a0 = gam*M2*(-Fr2 + omega^2*R0);
%     
% T1 = Ttop - (Ttop-Tbot)*exp(-(x1-R0)/h0);
% logp_p0 = a0*h0/Ttop*log(1+Ttop/Tbot*(exp((x1-R0)/h0)-1));
% rho1 = logp_p0 - log(T1);
% % rho1 = exp(rtilde);
% 
% vr1 = zeros(size(x1));
% rho1 = sol(:,1,:,1); 
% vr1 = sol(:,2,:,1)./sqrt(exp(rho1));
% T1 = sol(:,3,:,1)./sqrt(exp(rho1));
% 
% rhoend = sol(:,1,:,end); 
% vrend = sol(:,2,:,end)./sqrt(exp(rhoend));
% Tend = sol(:,3,:,end)./sqrt(exp(rhoend));
% 
% plot(dgnodes(:),T1(:),dgnodes(:),Tend(:),'linewidth',2);
% axis([min(dgnodes(:)) max(dgnodes(:)) -18 10])
% hold on
% grid on
% grid minor
% plot(dgnodes(:),vr1(:),dgnodes(:),vrend(:),'linewidth',2);
% plot(dgnodes(:),rho1(:),dgnodes(:),rhoend(:),'linewidth',2);
% 
% legend('$T_1$','$T_\infty$','$v_1$','$v_\infty$','$\rho_1$','$\rho_\infty$','location','southwest','Interpreter','latex','FontSize',18)
% 
% hold off
% 
% disp("Done!");
