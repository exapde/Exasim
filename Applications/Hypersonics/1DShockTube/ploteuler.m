%% plot solution
% sol = fetchsolution(pde,master,dmd, 'dataout');
for ti = pde.soltime
% ti = nt;
sol = getsolution(['dataout/out_t' num2str(ti)],dmd,master.npe);
% sol = getsolution(['dataout/out'],dmd,master.npe);
nspecies=1;
dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);    
for i = 1:nspecies 
    rho = sol(:,i,:);
%     subplot(1,3,1)
%     hold on
%     plot(dgnodes(:),u(:),'LineWidth',1.3)
end
rhou = sol(:,nspecies+1,:);
% subplot(1,3,2)
% plot(dgnodes(:),u(:),'LineWidth',1.3);

rhoE = sol(:,nspecies+2,:);
  
% u = rhou./rho;
% subplot(1,3,3)
p = eulereval(sol,'p',1.4,0.55);
% p = rho;   

u = rhou(:)./rho(:);
figure(1);
subplot(1,3,1); plot(dgnodes(:),rho(:),'LineWidth',1.3); title("rho")
subplot(1,3,2); plot(dgnodes(:),rhou(:),'LineWidth',1.3) ; title("rho u")
subplot(1,3,3); plot(dgnodes(:),p(:),'LineWidth',1.3) ; title("p")

gamma=1.4;
a = sqrt(gamma * p ./ rho);
figure(2);
plot(dgnodes(:),u(:) + a(:),dgnodes(:),u(:), dgnodes(:),u(:)-a(:),"LineWidth",1.3)
legend("u+a","u","u-a")

drawnow; 
% waitforbuttonpress; 
end
% plot


