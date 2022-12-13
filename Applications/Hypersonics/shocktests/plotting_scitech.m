%% laplacian
ttst = 0.96 * length(pde.dt);
solLa = getsolution(['laplacian_snapshots/out_t' num2str(ttst)],dmd,master.npe);
solLa_output = getsolution(['laplacian_snapshots/out_outputCG_t' num2str(ttst)],dmd,master.npe);

solPb1_5 = getsolution(['PB_snapshots_av15/out_t' num2str(ttst-20)],dmd,master.npe);
solPb1_5_output = getsolution(['PB_snapshots_av15/out_outputCG_t' num2str(ttst-20)],dmd,master.npe);

solPb10 = getsolution(['PB_snapshots/out_t' num2str(ttst-20)],dmd,master.npe);
solPb10_output = getsolution(['PB_snapshots/out_outputCG_t' num2str(ttst-20)],dmd,master.npe);
% 
solPbDi = getsolution(['PB_Di_snapshots/out_t' num2str(ttst-20)],dmd,master.npe);
solPbDi_output = getsolution(['PB_Di_snapshots/out_outputCG_t' num2str(ttst-20)],dmd,master.npe);
htr_out

%% Plot Lap
figure(10)
rhoLa = sum(solLa(:,1:5,:),2);
rhouLa = solLa(:,6,:);

pLa = solLa_output(:,1,:);

figure(13); hold on
rho1La = solLa(:,1,:);
rho2La = solLa(:,2,:);
rho3La = solLa(:,3,:);
rho4La = solLa(:,4,:);
rho5La = solLa(:,5,:);
plot(dgnodes(:),rho1La(:)./rhoLa(:),'LineWidth',1.3);
plot(dgnodes(:),rho2La(:)./rhoLa(:),'LineWidth',1.3);
plot(dgnodes(:),rho3La(:)./rhoLa(:),'LineWidth',1.3);
plot(dgnodes(:),rho4La(:)./rhoLa(:),'LineWidth',1.3);
plot(dgnodes(:),rho5La(:)./rhoLa(:),'LineWidth',1.3);
title("Species Mass Fractions");
grid on
legend(["Y_N" "Y_O" "Y_{NO}" "Y_{N2}" "Y_{O2}"])

figure(14); 
subplot(1,3,1); hold on
plot(dgnodes(:),rhoLa(:)*pde.externalparam(1),'LineWidth',1.3); 
scatter(rhoTrue(:,1), rhoTrue(:,2)*pde.externalparam(1))
xlim([0 1])
ylabel("\rho (kg / m^3)")
% title("Density")
legend(["p=3 sol", "reference sol"])
grid on

subplot(1,3,2); hold on
plot(dgnodes(:),rhouLa(:)./rhoLa(:),'LineWidth',1.3);
scatter(uTrue(:,1), uTrue(:,2))
xlim([0 1])
ylabel("u (m/s)")
% title("Velocity")
grid on

subplot(1,3,3); hold on
plot(dgnodes(:),pLa(:),'LineWidth',1.3);
scatter(pTrue(:,1), pTrue(:,2)*1e4)
xlim([0 1])
grid on
% title("Pressure")
ylabel("p (Pa)")

%% Plot PB 15
figure(100)
rhoPb1_5 = sum(solPb1_5(:,1:5,:),2);
rhouPb1_5 = solPb1_5(:,6,:);

pPb1_5 = solPb1_5_output(:,1,:);

figure(130); hold on
rho1Pb1_5 = solPb1_5(:,1,:);
rho2Pb1_5 = solPb1_5(:,2,:);
rho3Pb1_5 = solPb1_5(:,3,:);
rho4Pb1_5 = solPb1_5(:,4,:);
rho5Pb1_5 = solPb1_5(:,5,:);
plot(dgnodes(:),rho1Pb1_5(:)./rhoPb1_5(:),'LineWidth',1.3);
plot(dgnodes(:),rho2Pb1_5(:)./rhoPb1_5(:),'LineWidth',1.3);
plot(dgnodes(:),rho3Pb1_5(:)./rhoPb1_5(:),'LineWidth',1.3);
plot(dgnodes(:),rho4Pb1_5(:)./rhoPb1_5(:),'LineWidth',1.3);
plot(dgnodes(:),rho5Pb1_5(:)./rhoPb1_5(:),'LineWidth',1.3);
title("Species Mass Fractions");
grid on
legend(["Y_N" "Y_O" "Y_{NO}" "Y_{N2}" "Y_{O2}"])

figure(140); 
subplot(1,3,1); hold on
plot(dgnodes(:),rhoPb1_5(:)*pde.externalparam(1),'LineWidth',1.3); 
scatter(rhoTrue(:,1), rhoTrue(:,2)*pde.externalparam(1))
xlim([0 1])
ylabel("\rho (kg / m^3)")
% title("Density")
legend(["p=3 sol", "reference sol"])
grid on

subplot(1,3,2); hold on
plot(dgnodes(:),rhouPb1_5(:)./rhoPb1_5(:),'LineWidth',1.3);
scatter(uTrue(:,1), uTrue(:,2))
xlim([0 1])
ylabel("u (m/s)")
% title("Velocity")
grid on

subplot(1,3,3); hold on
plot(dgnodes(:),pPb1_5(:),'LineWidth',1.3);
scatter(pTrue(:,1), pTrue(:,2)*1e4)
xlim([0 1])
grid on
% title("Pressure")
ylabel("p (Pa)")

%% Plot PB 10
figure(1000)
rhoPb10 = sum(solPb10(:,1:5,:),2);
rhouPb10 = solPb10(:,6,:);

pPb10 = solPb10_output(:,1,:);

figure(1300); hold on
rho1Pb10 = solPb10(:,1,:);
rho2Pb10 = solPb10(:,2,:);
rho3Pb10 = solPb10(:,3,:);
rho4Pb10 = solPb10(:,4,:);
rho5Pb10 = solPb10(:,5,:);
plot(dgnodes(:),rho1Pb10(:)./rhoPb10(:),'LineWidth',1.3);
plot(dgnodes(:),rho2Pb10(:)./rhoPb10(:),'LineWidth',1.3);
plot(dgnodes(:),rho3Pb10(:)./rhoPb10(:),'LineWidth',1.3);
plot(dgnodes(:),rho4Pb10(:)./rhoPb10(:),'LineWidth',1.3);
plot(dgnodes(:),rho5Pb10(:)./rhoPb10(:),'LineWidth',1.3);
title("Species Mass Fractions");
grid on
legend(["Y_N" "Y_O" "Y_{NO}" "Y_{N2}" "Y_{O2}"])

figure(1400); 
subplot(1,3,1); hold on
plot(dgnodes(:),rhoPb10(:)*pde.externalparam(1),'LineWidth',1.3); 
scatter(rhoTrue(:,1), rhoTrue(:,2)*pde.externalparam(1))
xlim([0.2 0.7])
ylabel("\rho (kg / m^3)")
% title("Density")
legend(["p=3 sol", "reference sol"])
grid on

subplot(1,3,2); hold on
plot(dgnodes(:),rhouPb10(:)./rhoPb10(:),'LineWidth',1.3);
scatter(uTrue(:,1), uTrue(:,2))
xlim([0.2 0.7])
ylabel("u (m/s)")
% title("Velocity")
grid on

subplot(1,3,3); hold on
plot(dgnodes(:),pPb10(:),'LineWidth',1.3);
scatter(pTrue(:,1), pTrue(:,2)*1e4)
xlim([0.2 0.7])
grid on
% title("Pressure")
ylabel("p (Pa)")

%% Plot PB 10
figure(10000)
rhoPbDi = sum(solPbDi(:,1:5,:),2);
rhouPbDi = solPbDi(:,6,:);

pPbDi = solPbDi_output(:,1,:);

figure(13000); hold on
rho1PbDi = solPbDi(:,1,:);
rho2PbDi = solPbDi(:,2,:);
rho3PbDi = solPbDi(:,3,:);
rho4PbDi = solPbDi(:,4,:);
rho5PbDi = solPbDi(:,5,:);
plot(dgnodes(:),rho1PbDi(:)./rhoPbDi(:),'LineWidth',1.3);
plot(dgnodes(:),rho2PbDi(:)./rhoPbDi(:),'LineWidth',1.3);
plot(dgnodes(:),rho3PbDi(:)./rhoPbDi(:),'LineWidth',1.3);
plot(dgnodes(:),rho4PbDi(:)./rhoPbDi(:),'LineWidth',1.3);
plot(dgnodes(:),rho5PbDi(:)./rhoPbDi(:),'LineWidth',1.3);
title("Species Mass Fractions");
grid on
legend(["Y_N" "Y_O" "Y_{NO}" "Y_{N2}" "Y_{O2}"])

figure(1400); 
subplot(1,3,1); hold on
plot(dgnodes(:),rhoPbDi(:)*pde.externalparam(1),'--r','LineWidth',1.3); 
scatter(rhoTrue(:,1), rhoTrue(:,2)*pde.externalparam(1))
% xlim([0 1])
ylabel("\rho (kg / m^3)")
% title("Density")
legend(["p=3 sol", "reference sol"])
grid on

subplot(1,3,2); hold on
plot(dgnodes(:),rhouPbDi(:)./rhoPbDi(:),'--r','LineWidth',1.3);
scatter(uTrue(:,1), uTrue(:,2))
% xlim([0 1])
ylabel("u (m/s)")
% title("Velocity")
grid on

subplot(1,3,3); hold on
plot(dgnodes(:),pPbDi(:),'--r','LineWidth',1.3);
scatter(pTrue(:,1), pTrue(:,2)*1e4)
% xlim([0 1])
grid on
% title("Pressure")
ylabel("p (Pa)")


%% misc

figure(15)
subplot(1,2,1);
plot(dgnodes(:),rhouPb1_5(:)./rhoPb1_5(:),'LineWidth',1.3);
xlim([0 1])
grid on
ylabel("u (m/s)")
title("k_{\beta} = 1.5")
subplot(1,2,2)
plot(dgnodes(:),rhouPb10(:)./rhoPb10(:),'LineWidth',1.3);
xlim([0 1])
grid on
title("k_{\beta} = 10")
ylabel("u (m/s)")


figure(16);
subplot(1,3,1);
plot(dgnodes(:),rho2La(:)./rhoLa(:),'LineWidth',1.3);
% xlim([0.25 0.75])
grid on
ylabel("Y_{O}")
ylim([0 0.33])
subplot(1,3,2);
plot(dgnodes(:),rho2Pb10(:)./rhoPb10(:),'LineWidth',1.3);
ylim([0 0.33])
ylabel("Y_{O}")

grid on
% xlim([0.25 0.75])
subplot(1,3,3);
plot(dgnodes(:),rho2PbDi(:)./rhoPbDi(:),'LineWidth',1.3);
ylim([0 0.33])
ylabel("Y_{O}")

grid on
% xlim([0.25 0.75])

figure(17);
subplot(1,3,1);
plot(dgnodes(:),rho3La(:)./rhoLa(:),'LineWidth',1.3);
% xlim([0.25 0.75])
grid on
ylim([0 0.05])
ylabel("Y_{NO}")
subplot(1,3,2);
plot(dgnodes(:),rho3Pb10(:)./rhoPb10(:),'LineWidth',1.3);
ylim([0 0.05])
ylabel("Y_{NO}")
grid on
% xlim([0.25 0.75])
subplot(1,3,3);
plot(dgnodes(:),rho3PbDi(:)./rhoPbDi(:),'LineWidth',1.3);
ylim([0 0.05])
ylabel("Y_{NO}")
grid on

%%
figure(19); 
Atmp = get(gca,'colororder');
subplot(1,3,1); hold on
plot(dgnodes(:),rhoPb10(:)*pde.externalparam(1),'LineWidth',1.3); 
plot(dgnodes(:),rhoPbDi(:)*pde.externalparam(1),'--r','LineWidth',1.3); 
scatter(rhoTrue(:,1), rhoTrue(:,2)*pde.externalparam(1),'MarkerEdgeColor', Atmp(2,:))
xlim([0.2 0.7])
ylabel("\rho (kg / m^3)")
% title("Density")
legend(["AV \beta^*", "AV \beta^*, D_i^*","reference sol"])
grid on

subplot(1,3,2); hold on
plot(dgnodes(:),rhouPb10(:)./rhoPb10(:),'LineWidth',1.3);
plot(dgnodes(:),rhouPbDi(:)./rhoPbDi(:),'--r','LineWidth',1.3);
scatter(uTrue(:,1), uTrue(:,2),'MarkerEdgeColor', Atmp(2,:))
xlim([0.2 0.7])
ylabel("u (m/s)")
% title("Velocity")
grid on

subplot(1,3,3); hold on
plot(dgnodes(:),pPb10(:),'LineWidth',1.3);
plot(dgnodes(:),pPbDi(:),'--r','LineWidth',1.3);
scatter(pTrue(:,1), pTrue(:,2)*1e4,'MarkerEdgeColor', Atmp(2,:))
xlim([0.2 0.7])
grid on
% title("Pressure")
ylabel("p (Pa)")
% 
% subplot(1,3,1); hold on
% plot(dgnodes(:),rhoPbDi(:)*pde.externalparam(1),'--r','LineWidth',1.3); 
% scatter(rhoTrue(:,1), rhoTrue(:,2)*pde.externalparam(1))
% % xlim([0 1])
% ylabel("\rho (kg / m^3)")
% % title("Density")
% legend(["p=3 sol", "reference sol"])
% grid on

% subplot(1,3,2); hold on
% plot(dgnodes(:),rhouPbDi(:)./rhoPbDi(:),'--r','LineWidth',1.3);
% scatter(uTrue(:,1), uTrue(:,2))
% % xlim([0 1])
% ylabel("u (m/s)")
% % title("Velocity")
% grid on
% 
% subplot(1,3,3); hold on
% plot(dgnodes(:),pPbDi(:),'--r','LineWidth',1.3);
% scatter(pTrue(:,1), pTrue(:,2)*1e4)
% % xlim([0 1])
% grid on
% % title("Pressure")
% ylabel("p (Pa)")
