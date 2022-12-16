%% plot solution
% sol = fetchsolution(pde,master,dmd, 'dataout');
% for ti = pde.saveSolFreq:pde.saveSolFreq:(length(pde.dt)*6)
for ti = pde.saveSolFreq:pde.saveSolFreq:(50000)
%     sol = getsolution(['dataout/out_t' num2str(ti)],dmd,master.npe);
sol = getsolution(['dataout/out_outputCG_t' num2str(ti)],dmd,master.npe);
%     solw = getsolution(['dataout/out_wdg_t' num2str(ti)],dmd,master.npe);
    dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);    
            
%     rho = sum(sol(:,1:5,:),2);
    for i = 1
        u = sol(:,i,:);
%         subplot(1,3,1); 
        figure(i)
        clf

        
        hold on
        
%         plot(dgnodes(:),u(:)./rho(:),'LineWidth',1.3); 
% 
plot(dgnodes(:),rho(:),'LineWidth',1.3); 
set(gca, 'Xscale','log'); 
%         waitforbuttonpress
    end
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
    waitforbuttonpress
%     disp(ti)
end


        
