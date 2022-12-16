figure(10)
ti = 129000;

dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);    
sol = getsolution(['dataout/out_t' num2str(ti)],dmd,master.npe);
rho = sum(sol(:,1:5,:),2);
Y1 = sol(:,1,:)./rho;
Y2 = sol(:,2,:)./rho;
Y3 = sol(:,3,:)./rho;
Y4 = sol(:,4,:)./rho;
Y5 = sol(:,5,:)./rho;
cgsol = getsolution(['dataout/out_outputCG_t' num2str(ti)],dmd,master.npe);
Pout = cgsol(:,1,:);
Tout = cgsol(:,2,:);


subplot(1,3,1)
hold on
plot(dgnodes(:),Y1(:),'b','LineWidth',1.3); 
plot(dgnodes(:),Y2(:),'r','LineWidth',1.3); 
plot(dgnodes(:),Y3(:),'g','LineWidth',1.3); 
plot(dgnodes(:),Y4(:),'k','LineWidth',1.3); 
plot(dgnodes(:),Y5(:),'m','LineWidth',1.3); 

scatter(xY1true,Y1true2,'bo','LineWidth',1.3)
scatter(xY2true,Y2true2,'ro','LineWidth',1.3)
scatter(xY3true,Y3true2,'go','LineWidth',1.3)
scatter(xY4true,Y4true2,'ko','LineWidth',1.3)
scatter(xY5true,Y5true2,'mo','LineWidth',1.3)
set(gca,'Xscale','log');
xlim([3e-5, 1.0])
ylabel("Y_i")
xlabel("x (m)")
title("Mass Fractions")
legend("Y_N","Y_O","Y_{NO}","Y_{N2}","Y_{O2}")
grid on

subplot(1,3,2)
hold on
plot(dgnodes(:),Tout(:),'b','LineWidth',1.3)
scatter(xTtrue,Ttrue2,'bo','LineWidth',1.3)
set(gca,'Xscale','log');
xlim([3e-5, 1.0])
ylabel("T (K)")
xlabel("x (m)")
title("Temperature")
legend("Temperature", "Reference Temperature")
grid on

subplot(1,3,3)
hold on
plot(dgnodes(:),Pout(:),'b','LineWidth',1.3)
scatter(xPtrue,Ptrue2,'bo','LineWidth',1.3)
set(gca,'Xscale','log');
xlim([3e-5, 1.0])
title("Pressure")
xlabel("x (m)")
legend("Pressure", "Reference Pressure")
ylabel("P (Pa)")
grid on

