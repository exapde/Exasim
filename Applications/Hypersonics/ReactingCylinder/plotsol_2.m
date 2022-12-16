mesh.nd = master.dim;
mesh.plocal = master.xpe;
mesh.tlocal = master.telem;
mesh.porder = pde.porder;
% soldir = 'nanSnapshots/dataout';
% soldir = 'check10kOf100k';
% soldir = 'checkSkChange';
% soldir = 'checkRefineBoundaryLayer';
% soldir = 'checkUniRef';
% soldir = 'checkHhChange';
% note: first 100k dt of 2.5e-5
%       next 10k dt of 2.5e-4
% 
% soldir = 'overnightN80Run';
% soldir = 'checkChangeLscale';
% soldir = 'checkDomainChange';
% soldir = 'checkDomainChangeAndLrefChange';
% soldir = 'checkNoSource'
% soldir = 'checkp1'
% soldir = 'checkp1NoSource'; % note: i might have also changed hh here
% soldir = 'checkp1NoSourceN100';
% soldir = 'checkp1NoSourceEinitChange';
% soldir = 'checkp1NoSourceEinitChangeKbChange';
% soldir = 'p1N100RegularGrid';
% soldir = 'p1N200RegularGrid';
% soldir = 'checkAppScp_AccidentalSource';
% soldir = 'checkAppScp';
% soldir = 'sourceTermTest'
% ti = 10000;
% soldir = 'sourceTermTest_halfdt'
% soldir = 'sourceTermTest_halfdt_sb0change'
% soldir = 'sourceTermTest_halfdt_n100'
% soldir = 'sourceTermTest_adaptdt'
% soldir = 'sourceTermTest_gmrestol'
% soldir = 'sourceTermTest_Mach59';
% soldir = 'sourceTermTest_Mach59_testEps';
% soldir = 'sourceTermTest_Mach59_testEps1e9';
% soldir = 'sourceTermTest_Mach59_overnight'
% soldir = 'p2n80Ma8';
% soldir = 'p2n80Ma8_smalldt';
% soldir = 'p1n80Ma8';
soldir = 'dataout'

ti = 5;

% 50k of 2.5e-5
% 10k of 1.5e-4

% Uout = getsolution(['dataout/out'],dmd,master.npe);
Uout = getsolution([strcat(soldir,'/out_t') num2str(ti)],dmd,master.npe);
UCGout = getsolution([strcat(soldir,'/out_outputCG_t') num2str(ti)],dmd,master.npe);
% out = eulereval(Uout,'r',gam,Minf);
rho = sum(Uout(:,1:5,:),2);
rhou = Uout(:,6,:);
rhov = Uout(:,7,:);
rhoE = Uout(:,8,:);
drho_dx = sum(Uout(:,9:13,:),2);
drho_dy = sum(Uout(:,17:21,:),2);
for i = 17:21
    out = Uout(:,i,:);
    disp(max(out(:))/min(out(:)));
end
disp(max(drho_dy(:))/min(drho_dy(:)) )

out = Uout(:,24,:)./rho;
% out = drho_dy;
% out = rho;
% out = rho;
% [om,ind] = max(abs(out(:)));
figure(100); clf; scaplot(mesh,out,[],1);
ylim([-1.1, 1.1])
xlim([-1 0.2])
% x = mesh.dgnodes(:,1,:);
% y = mesh.dgnodes(:,2,:);
% hold on;
% plot(x(ind),y(ind),'o');
colormap('jet');
% %%
%%
figure(50);
species = ["N","O","NO","N2","O2"];
for i = 1:5
%     if i ~= 9
        figure(i); clf
% subplot(1,5,i)
        outcurr = Uout(:,i,:);
        outcurr = outcurr(:)./rho(:);
        scaplot(mesh,outcurr,[],1); colormap('jet')
%         hold on; plot(linspace(-3, 1), 0*linspace(-3,1))
%         xlim([-3 0])
%         ylim([0 3])
%         set(gca, "visible", "off")
axis("off")
       title(strcat("Y_{", string(species(i)), string("}")))
       drawnow

%         waitforbuttonpress
%     end
end

%%
titles = ["w1","w2","w3","w4","w5","p","t","gam","rhoe","sb0"];
for i = 1:5
%     if i ~= 9
        figure(i); clf
        wcurr = UCGout(:,i,:);
        scaplot(mesh,wcurr*0.045,[],1);
        title(titles(i))
%         waitforbuttonpress
        drawnow
%     end
end
%%
P = UCGout(:,6,:);
T = UCGout(:,7,:);
figure(200); scaplot(mesh,T,[],1,0);
% hold on; plot(linspace(-3, 1), 0*linspace(-3,1))
colormap('jet')
title("Temperature")
gam = UCGout(:,8,:);
sb0 = UCGout(:,10,:);
rhoe = UCGout(:,9,:);
figure(300); scaplot(mesh, gam, [], 1, 0);
colormap('jet')
figure(301); scaplot(mesh, sb0, [], 1, 0);
colormap('jet')
figure(400); scaplot(mesh, rhoe, [], 1, 0); title('rhoe')
colormap('jet')

%%
figure(5000); scaplot(mesh, rho, [], 1, 0); title("rho")
colormap('jet')
figure(6000); scaplot(mesh, rhou./rho, [], 1, 0); title("u");
colormap('jet')
figure(7000); scaplot(mesh, rhov./rho, [], 1, 0); title("v");
colormap('jet')
figure(8000); scaplot(mesh, Uout(:,8,:), [], 1, 0); title("rhoE");

disp("rho min/max")
disp(min(rho(:)))
disp(max(rho(:)))

disp("u min/max")
disp(min(rhou(:)./rho(:)))
disp(max(rhou(:)./rho(:)))

disp("v min/max")
disp(min(rhov(:)./rho(:)))
disp(max(rhov(:)./rho(:)))

disp("rhoE min/max")
disp(min(rhoE(:)))
disp(max(rhoE(:)))

disp("rhoe min/max")
disp(min(rhoe(:)))
disp(max(rhoe(:)))

% figure(1)
% out = Uout(:,8,:);
% % outCG = 
% [om,ind] = max(abs(out(:)));
% figure(1); clf; scaplot(mesh,out,[],1);
% x = mesh.dgnodes(:,1,:);
% y = mesh.dgnodes(:,2,:);
% hold on;
% plot(x(ind),y(ind),'o');
% colormap('jet');

% 
% Uout = getsolution(['dataout/out_odg' num2str(ti)],dmd,master.npe);
% figure(2); clf; scaplot(mesh,Uout(:,1,:),[],1); colormap('jet');
% out=Uout(:,1,:); [om,ind] = max(abs(out(:)));
% x = mesh.dgnodes(:,1,:);
% y = mesh.dgnodes(:,2,:);
% hold on;
% plot(x(ind),y(ind),'o');
% e=Uout(:,2,:)-mesh.vdg(:,2,:); max(abs(e(:)))
% 
% Uout = getsolution(['dataout/out_smoothedodg' num2str(ti)],dmd,master.npe);
% figure(4); clf; scaplot(mesh,Uout(:,1,:),[],1); colormap('jet');
% out=Uout(:,1,:); [om,ind] = max(abs(out(:)));
% x = mesh.dgnodes(:,1,:);
% y = mesh.dgnodes(:,2,:);
% hold on;
% plot(x(ind),y(ind),'o');
% e=Uout(:,2,:)-mesh.vdg(:,2,:); max(abs(e(:)))
% 
% udg = getsolution(['dataout/out_udg' num2str(ti)],dmd,master.npe);
% vdg = getsolution(['dataout/out_odg' num2str(ti)],dmd,master.npe);
% sdg = getsolution(['dataout/out_smoothedodg' num2str(ti)],dmd,master.npe);
% avField = getavfield2dmatlab(udg(:,1:4,:),udg(:,5:12,:),vdg,pde.physicsparam);
% avc = getavfield2dc(udg,vdg,pde.physicsparam);
% 
% figure(1); clf; scaplot(mesh,avField(:,1,:),[],1); colormap('jet');
% figure(2); clf; scaplot(mesh,avc(:,1,:),[],1); colormap('jet');
% figure(3); clf; scaplot(mesh,vdg(:,1,:),[],1); colormap('jet');
% figure(4); clf; scaplot(mesh,sdg(:,1,:),[],1); colormap('jet');
% 
% 
% % figure(1); clf; hold on;
% % for ib=1:3
% %     f = mkf(mesh.t,mesh.f,2);
% %     if ib==1
% %         pars={'facecolor',[.8,1,.8],'edgecolor','b','Linew',1,'FaceAlpha',1,'EdgeAlpha',1};
% %     elseif ib==2
% %         pars={'facecolor',[.8,1,.8],'edgecolor','r','Linew',1,'FaceAlpha',1,'EdgeAlpha',1};
% %     else
% %         pars={'facecolor',[.8,1,.8],'edgecolor','k','Linew',1,'FaceAlpha',1,'EdgeAlpha',1};
% %     end
% %     ind = find(f(end,:)==-ib);
% %     patch('faces',f(1:end-4,ind)','vertices',mesh.p',pars{:});                            
% % end
% % axis equal; axis tight;
% % 
% % 
% 
% 
% %%
% xtmp = []
% tmp = []
% dgtst = reshape(mesh.dgnodes, [6400*4 2]);
% Ttst = reshape(T, [6400*4 1]);
% for i = 1:size(dgtst,1)
%     if abs(dgtst(i,2)) < 1e-4
%         disp(dgtst(i,1))
%         xtmp = [xtmp, dgtst(i,1)];
%         tmp = [tmp, Ttst(i)];
%     end
% end

%%%
% no viscous: 
% rho min/max
%    0.999999946808838
% 
%    1.000590779570255
% 
% u min/max
%      1.208341179345986e-04
% 
%    1.000000019931180
% 
% v min/max
%     -4.684172697272210e-05
% 
%      4.684173134558622e-05
% 
% rhoE min/max
%    0.111810873353252
% 
%    0.611594698012404
% 
% rhoe min/max
%      6.117602744265359e+03
% 
%      6.231298010549617e+03

% rho min/max
%    0.999999913817226
% 
%    1.000420041679332
% 
% u min/max
%    0.041600470690700
% 
%    1.000085250338504
% 
% v min/max
%   -0.023404378404112
% 
%    0.023404375771517
% 
% rhoE min/max
%    0.267275124566171
% 
%    0.611679971000627
% 
% rhoe min/max
%      6.122020659204081e+03
% 
%      1.771234360440572e+04


