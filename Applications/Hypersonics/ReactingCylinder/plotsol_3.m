mesh.nd = master.dim;
mesh.plocal = master.xpe;
mesh.tlocal = master.telem;
mesh.porder = pde.porder;

soldir = 'dataout'; 
% ti = 52000; 
ti = 4900+382;
% 50k of 2.5e-5
% 10k of 1.5e-4

% Uout = getsolution(['dataout/out'],dmd,master.npe);
Uout = getsolution([strcat(soldir,'/out_t') num2str(ti)],dmd,master.npe);
UCGout = getsolution([strcat(soldir,'/out_outputCG_t') num2str(ti)],dmd,master.npe);
% out = eulereval(Uout,'r',gam,Minf);
rho = sum(Uout(:,1:5,:),2);
drho_dx = sum(Uout(:,9:13,:),2);
drho_dy = sum(Uout(:,17:21,:),2);
for i = 17:21
    out = Uout(:,i,:);
    disp(max(out(:))/min(out(:)));
end
disp(max(drho_dy(:))/min(drho_dy(:)) )

% out = Uout(:,8,:);
% out = drho_dy;
out = rho;
% out = rho;
% [om,ind] = max(abs(out(:)));
figure(100); clf; scaplot(mesh,out,[]);
% ylim([-1.1, 1.1])
% xlim([-1 0.2])
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
        scaplot(mesh,outcurr,[]); colormap('jet')
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
for i = 1:10
%     if i ~= 9
        figure(i); clf
        wcurr = UCGout(:,i,:);
        scaplot(mesh,wcurr,[],1);
        colormap('jet')
        title(titles(i))
        drawnow
%         waitforbuttonpress
%     end
end
%%
P = UCGout(:,6,:);
T = UCGout(:,7,:);
figure(200); scaplot(mesh,T,[],1,0); title("T")
% hold on; plot(linspace(-3, 1), 0*linspace(-3,1))
colormap('jet')
title("Temperature")
gam = UCGout(:,8,:);
sb0 = UCGout(:,10,:);
rhoe = UCGout(:,9,:);
figure(300); scaplot(mesh, gam, [], 1, 0); title("gamma"); drawnow
colormap('jet')
figure(301); scaplot(mesh, sb0, [], 1, 0); title("sb0"); drawnow
colormap('jet')
figure(400); scaplot(mesh, rhoe, [], 1, 0); title("re"); drawnow
colormap('jet')


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
%%
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
ncu = 5 + 2 + 1;
rhoinv = 1./rho;
uv = Uout(:,6,:);
vv = Uout(:,7,:);
drhou_dx = Uout(:,ncu+5+1,:);
drhov_dx = Uout(:,ncu+5+2,:);
drhou_dy = Uout(:,ncu+ncu+5+1,:);
drhov_dy = Uout(:,ncu+ncu+5+2,:);
du_dx = (drhou_dx - drho_dx.*uv).*rhoinv;
dv_dx = (drhov_dx - drho_dx.*vv).*rhoinv;
du_dy = (drhou_dy - drho_dy.*uv).*rhoinv;
dv_dy = (drhov_dy - drho_dy.*vv).*rhoinv;
div_v = du_dx + dv_dy;
vort = (dv_dx - du_dy);
vort = sqrt(vort.*vort);
 DucrosRatio = div_v.*div_v ./ (div_v.*div_v + vort.*vort + 1.0e-16);
figure(1); clf; scaplot(mesh, div_v, [], 1, 0);
figure(2); clf; scaplot(mesh, vort, [], 1, 0);
figure(3); clf; scaplot(mesh, DucrosRatio, [], 1, 0);