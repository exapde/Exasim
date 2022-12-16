mesh.nd = master.dim;
mesh.plocal = master.xpe;
mesh.tlocal = master.telem;
mesh.porder = pde.porder;
ti = 10;
% with dt=2.5e-5: ti = 10->-0.95
%                 ti = 100->-0.7558
% with dt=1e-5:   ti = 20->-0.9626
%                 ti = 200->-0.7855
% with dt=1e-6:   ti = 200-> -0.9629
%                 ti = 2000->-0.7873
% 
% soldir = 'nanSnapshots/dataout';
% soldir = 'check10kOf100k';
% soldir = 'checkSkChange';
% soldir = 'checkRefineBoundaryLayer';
soldir = 'dataout';
% note: first 100k dt of 2.5e-5
%       next 10k dt of 2.5e-4

%    -0.9635
% 
%    -0.8788
% 
%    -0.8804
% 
%    -0.9596
% 
%    -0.9616
% 
%    -0.9616



% Uout = getsolution(['dataout/out'],dmd,master.npe);
Uout = getsolution([strcat(soldir,'/out_t') num2str(ti)],dmd,master.npe);
UCGout = getsolution([strcat(soldir,'/out_outputCG_t') num2str(ti)],dmd,master.npe);
% out = eulereval(Uout,'r',gam,Minf);
rho = sum(Uout(:,1:5,:),2);
drho_dx = sum(Uout(:,9:13,:),2);

drho_dy = sum(Uout(:,17:21,:),2);
for i = 17:21
    out = Uout(:,i,:);
    disp(max(out(:))/min(out(:)) )
end
% out = drho_dy;
% out = rho;
% out = rho;
% outCG = 
[om,ind] = max(abs(out(:)));
disp(max(out(:))/min(out(:)) )
%%
% disp(max(drho_dx(:))/min(drho_dx(:)) )
% disp(max(drho_dy(:))/min(drho_dy(:)) )
figure(100); clf; scaplot(mesh,out,[],1);
% ylim([-1.1, 1.1])
% xlim([-1 0.2])
x = mesh.dgnodes(:,1,:);
y = mesh.dgnodes(:,2,:);
hold on;
plot(x(ind),y(ind),'o');
colormap('jet');
%%
titles = ["w1","w2","w3","w4","w5","p","t","gam","rhoe","sb0"];
% for i = 1:10
% %     if i ~= 9
%         figure(i); clf
%         wcurr = UCGout(:,i,:);
%         scaplot(mesh,wcurr,[],1);
%         title(titles(i))
% %         waitforbuttonpress
% %     end
% end
% figure(7); clf
P = UCGout(:,6,:);
T = UCGout(:,7,:);
figure(2); scaplot(mesh,T / 901,[],1,0);
colormap('jet')
title("T")
gam = UCGout(:,8,:);
sb0 = UCGout(:,10,:);
figure(3); scaplot(mesh, sb0, [], 1, 0);
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
