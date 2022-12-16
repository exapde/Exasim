%% plot solution
% sol = fetchsolution(pde,master,dmd, 'dataout');
% for ti = pde.saveSolFreq:pde.saveSolFreq:(length(pde.dt)*6)
hoverp = mesh.vdg(1,2,1) / pde.porder;
hmat = (mesh.dgnodes(2:end) - mesh.dgnodes(1:end-1));
for ti = pde.saveSolFreq:pde.saveSolFreq:length(pde.dt)
% for ti = 100:100:10000
% for ti = 10000:10:12000
    sol = getsolution(['dataout/out_t' num2str(ti)],dmd,master.npe);
    sol_output = getsolution(['dataout/out_outputCG_t' num2str(ti)],dmd,master.npe);
%     solw = getsolution(['dataout/out_wdg_t' num2str(ti)],dmd,master.npe);
    dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);    
%     rho = sum(sol(:,1:5,:),2);
%     drho = sum(sol(:,8:12,:),2);
% %     AV = sol_output(:,1,:);
%     u = sol(:,1,:);
%     du1 = sol(:, 8,:);
%     figure(11); clf; hold on
%     plot(dgnodes(:),drho(:),'LineWidth',1.3); 
%     plot(dgnodes(:),abs(sT(:))/max(abs(sT(:))),'LineWidth',1.3); 
%     ylim([0 1.2])
%     w1 = sol_output(:,3,:);     
%     w2 = sol_output(:,4,:);
%     w3 = sol_output(:,5,:);
%     w4 = sol_output(:,6,:);
%     w5 = sol_output(:,end,:);

% Why don't we check the nondimensional version. 
%     sb = sol_output(:,1,:);
    sT = sol_output(:,2,:);
    T = sol_output(:,8,:);
    dTdx = sT / hoverp;
    dTdxFD = (T(2:end) - T(1:end-1)) ./ (hmat);% * mesh.vdg(1,2,1);
    figure(1); clf
    plot(dgnodes(:),dTdx(:),'LineWidth',1.3); 
    
    figure(2); clf
    plot(dgnodes(:),T(:),'LineWidth',1.3); 

    figure(3); clf; hold on
    plot(dgnodes(1:end-1),dTdx(1:end-1),'LineWidth',1.3); 
    plot(dgnodes(1:end-1),dTdxFD(:),'LineWidth',1.3); 

    Pr = sol_output(:,10,:);
    figure(4); clf; hold on
    plot(dgnodes(:), Pr(:), 'LineWidth',1.4)
%     drawnow
%     waitforbuttonpress

%     sY1 = sol_output(:,3,:);
%     sY2 = sol_output(:,4,:);
%     sY3 = sol_output(:,5,:);
%     sY4 = sol_output(:,6,:);
%     sY5 = sol_output(:,7,:);
%     figure(1); clf; hold on
%     plot(dgnodes(:),abs(sb(:))/max(abs(sb(:))),'LineWidth',1.3); 
%     plot(dgnodes(:),abs(sT(:))/max(abs(sT(:))),'-*r','LineWidth',1.3); 
%     ylim([0 1.2])
% 
% %     figure(3); clf; hold on;
%     plot(dgnodes(:),abs(sY1(:))/max(abs(sY1(:))),'-.k','LineWidth',1.3);
%     plot(dgnodes(:),abs(sY2(:))/max(abs(sY2(:))),'-.k','LineWidth',1.3); 
%     plot(dgnodes(:),abs(sY3(:))/max(abs(sY3(:))),'-.k','LineWidth',1.3); 
%     plot(dgnodes(:),abs(sY4(:))/max(abs(sY4(:))),'-.k','LineWidth',1.3); 
%     plot(dgnodes(:),abs(sY5(:))/max(abs(sY5(:))),'-.k','LineWidth',1.3); 
%     ylim([0, 10000])
% xlim([0.45, 0.65])

    drawnow
    waitforbuttonpress
    disp(ti)
end
% final



        
