%% plot solution
% sol = fetchsolution(pde,master,dmd, 'dataout');
% for ti = pde.saveSolFreq:pde.saveSolFreq:(length(pde.dt)*6)
for ti = pde.saveSolFreq:pde.saveSolFreq:length(pde.dt)
% for ti = 100:100:10000
% for ti = 10000:10:12000
    sol = getsolution(['dataout/out_t' num2str(ti)],dmd,master.npe);
    sol_output = getsolution(['dataout/out_outputCG_t' num2str(ti)],dmd,master.npe);
%     solw = getsolution(['dataout/out_wdg_t' num2str(ti)],dmd,master.npe);
    dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);    
    rho = sum(sol(:,1:5,:),2);
    rho1 = sol(:,1,:);     
    rho2 = sol(:,2,:);
    rho3 = sol(:,3,:);
    rho4 = sol(:,4,:);
    rho5 = sol(:,5,:);
    rhou = sol(:,6,:);
    rhoE = sol(:,7,:);
    drho = sum(sol(:,8:12,:),2);
%     AV = sol_output(:,1,:);
    uplt = sol(:,6,:);
    du1 = sol(:, 8,:);
    figure(11); clf; hold on
    plot(dgnodes(:),uplt(:),'LineWidth',1.3); 
%     plot(dgnodes(:),abs(sT(:))/max(abs(sT(:))),'LineWidth',1.3); 
%     ylim([0 1.2])
%     w1 = sol_output(:,3,:);     
%     w2 = sol_output(:,4,:);
%     w3 = sol_output(:,5,:);
%     w4 = sol_output(:,6,:);
%     w5 = sol_output(:,end,:);


    sb = sol_output(:,1,:);
    sT = sol_output(:,2,:);
    sY1 = sol_output(:,3,:);
    sY2 = sol_output(:,4,:);
    sY3 = sol_output(:,5,:);
    sY4 = sol_output(:,6,:);
    sY5 = sol_output(:,7,:);
    figure(1); clf; hold on
    plot(dgnodes(:),abs(sb(:)),'LineWidth',1.3); 
    plot(dgnodes(:),abs(sT(:))/5,'-r','LineWidth',1.3); 
    ylim([0 2])

%     figure(3); clf; hold on;
    plot(dgnodes(:),abs(sY1(:)),'-k','LineWidth',1.3);
    plot(dgnodes(:),abs(sY2(:)),'-k','LineWidth',1.3); 
    plot(dgnodes(:),abs(sY3(:)),'-k','LineWidth',1.3); 
    plot(dgnodes(:),abs(sY4(:)),'-k','LineWidth',1.3); 
    plot(dgnodes(:),abs(sY5(:)),'-k','LineWidth',1.3); 
%     ylim([0, 10000])
xlim([0.45, 0.65])
legend("div v", "dTdx", "dYdx")

%     scatter(rhoTrue(:,1), rhoTrue(:,2))
    figure(50); 
    clf;
    hold on;
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
    ylim([-0.1 1])


% figure(12); clf; hold on
% p = sol_output(:,1,:);
% plot(dgnodes(:),p(:),'LineWidth',1.3);
% scatter(pTrue(:,1), pTrue(:,2)*1e4)

    drawnow
%     waitforbuttonpress
    disp(ti)
end
% final
%%
% ttst = 27*120;
% ttst = 324
ttst = ti;
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

figure(14); 
subplot(1,3,1); hold on
plot(dgnodes(:),rho(:)*pde.externalparam(1),'LineWidth',1.3); 
scatter(rhoTrue(:,1), rhoTrue(:,2)*pde.externalparam(1))
title("\rho")

subplot(1,3,2); hold on
plot(dgnodes(:),rhou(:)./rho(:),'LineWidth',1.3);
scatter(uTrue(:,1), uTrue(:,2))
title("u")

subplot(1,3,3); hold on
plot(dgnodes(:),p(:),'LineWidth',1.3);
scatter(pTrue(:,1), pTrue(:,2)*1e4)
title("p")


        
