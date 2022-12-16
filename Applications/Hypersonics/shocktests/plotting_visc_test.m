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

    D1 = sol_output(:,1,:);
    D2 = sol_output(:,2,:);
    D3 = sol_output(:,3,:);
    D4 = sol_output(:,4,:);
    D5 = sol_output(:,5,:);
    h1 = sol_output(:,6,:);
    h2 = sol_output(:,7,:);
    h3 = sol_output(:,8,:);
    h4 = sol_output(:,9,:);
    h5 = sol_output(:,10,:);
    mu = sol_output(:,11,:);
    kappa = sol_output(:,12,:);


%     T = sol_output(:,8,:);
    dTdx = sT / hoverp;
    dTdxFD = (T(2:end) - T(1:end-1)) ./ (hmat);% * mesh.vdg(1,2,1);
    figure(1); clf
    subplot(5, 1, 1)
    plot(dgnodes(:),D1(:),'LineWidth',1.3); 
    subplot(5, 1, 2)
    plot(dgnodes(:),D2(:),'LineWidth',1.3); 
    subplot(5, 1, 3)
    plot(dgnodes(:),D3(:),'LineWidth',1.3); 
    subplot(5, 1, 4)
    plot(dgnodes(:),D4(:),'LineWidth',1.3); 
    subplot(5, 1, 5)
    plot(dgnodes(:),D5(:),'LineWidth',1.3); 

    figure(2); clf
    subplot(5, 1, 1)
    plot(dgnodes(:),h1(:),'LineWidth',1.3); 
    subplot(5, 1, 2)
    plot(dgnodes(:),h2(:),'LineWidth',1.3); 
    subplot(5, 1, 3)
    plot(dgnodes(:),h3(:),'LineWidth',1.3); 
    subplot(5, 1, 4)
    plot(dgnodes(:),h4(:),'LineWidth',1.3); 
    subplot(5, 1, 5)
    plot(dgnodes(:),h5(:),'LineWidth',1.3); 

    figure(3); clf
    plot(dgnodes(:),mu(:),'LineWidth',1.3);
    figure(4); clf
    plot(dgnodes(:), kappa(:),'LineWidth',1.3);

%     
%     figure(2); clf
%     plot(dgnodes(:),T(:),'LineWidth',1.3); 
% 
%     figure(3); clf; hold on
%     plot(dgnodes(1:end-1),dTdx(1:end-1),'LineWidth',1.3); 
%     plot(dgnodes(1:end-1),dTdxFD(:),'LineWidth',1.3); 
% 
%     Pr = sol_output(:,10,:);
%     figure(4); clf; hold on
%     plot(dgnodes(:), Pr(:), 'LineWidth',1.4)

    drawnow
%     waitforbuttonpress
    disp(ti)
end
% final



        
