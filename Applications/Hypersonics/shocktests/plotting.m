%% plot solution
% sol = fetchsolution(pde,master,dmd, 'dataout');
% for ti = pde.saveSolFreq:pde.saveSolFreq:(length(pde.dt)*6)
for ti = pde.saveSolFreq:pde.saveSolFreq:length(pde.dt)
    sol = getsolution(['dataout/out_t' num2str(ti)],dmd,master.npe);
    sol_output = getsolution(['dataout/out_outputCG_t' num2str(ti)],dmd,master.npe);
%     solw = getsolution(['dataout/out_wdg_t' num2str(ti)],dmd,master.npe);
    dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);    
            
%     rho = sum(sol(:,1:5,:),2);
%     AV = sol_output(:,1,:);
    u = sol(:,1,:);
    figure(1)
    clf
    hold on
    plot(dgnodes(:),u(:),'LineWidth',1.3); 
%         figure(2)
%     clf
%     hold on
%     plot(dgnodes(:),AV(:),'LineWidth',1.3); 
% set(gca, 'Xscale','log'); 
%         waitforbuttonpress
%     end
%     rhou = sol(:,6,:);
%     plot(dgnodes(:),rhou(:)./rho(:),'LineWidth',1.3); set(gca, 'Xscale','log'); 

%     rhoE = sol(:,7,:);
%     p = solw(:,6,:);
%     for i = 1:5
%         u = solw(:,i,:);
%         hold on
%         plot(dgnodes(:),u(:),'LineWidth',1.3); set(gca, 'Xscale','log'); 
%         title("\omega_i")
%     end
%     subplot(1,3,2); 
%     plot(dgnodes(:), rhou(:)./rho(:),'LineWidth',1.3);    set(gca, 'Xscale','log')
%     subplot(1,3,3); 
%     plot(dgnodes(:), rhoE(:),'LineWidth',1.3); 
%     set(gca, 'Xscale','log')

%     u = sol(:,6,:);
%     subplot(1,3,2)
%     plot(dgnodes(:),u(:),'LineWidth',1.3);
%     
%     u = sol(:,7,:);
%     subplot(1,3,3)
%     plot(dgnodes(:), u(:),'LineWidth',1.3); 
%     
%     solw = getsolution(['dataout/out_wdg_t' num2str(ti)],dmd,master.npe);
%     w = solw(:,6,:) * rho_inf * u_inf^2;
%     disp(w(1));
%     figure(2);
%     semilogx(dgnodes(:),w(:));
    drawnow
%     waitforbuttonpress
    disp(ti)
end
% final
%%
% ttst = 27*120;
% ttst = 324
ttst = 500;
sol = getsolution(['dataout/out_t' num2str(ttst)],dmd,master.npe);
sol_output = getsolution(['dataout/out_outputCG_t' num2str(ttst)],dmd,master.npe);

htr_out
figure(10)
rho = sum(sol(:,1:5,:),2);
rhou = sol(:,6,:);
figure(10); hold on
plot(dgnodes(:),rho(:),'LineWidth',1.3); 
scatter(rhoTrue(:,1), rhoTrue(:,2))

figure(11); hold on
plot(dgnodes(:),rhou(:)./rho(:),'LineWidth',1.3);
scatter(uTrue(:,1), uTrue(:,2))


figure(12); hold on
p = sol_output(:,1,:);
plot(dgnodes(:),p(:),'LineWidth',1.3);
scatter(pTrue(:,1), pTrue(:,2)*1e4)

figure(13); hold on
rho1 = sol(:,1,:);
rho2 = sol(:,2,:);
rho3 = sol(:,3,:);
rho4 = sol(:,4,:);
rho5 = sol(:,5,:);
plot(dgnodes(:),rho1(:)./rho(:),'LineWidth',1.3);
plot(dgnodes(:),rho2(:)./rho(:),'LineWidth',1.3);
plot(dgnodes(:),rho3(:)./rho(:),'LineWidth',1.3);
plot(dgnodes(:),rho4(:)./rho(:),'LineWidth',1.3);
plot(dgnodes(:),rho5(:)./rho(:),'LineWidth',1.3);

        
