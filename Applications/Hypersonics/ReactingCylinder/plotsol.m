mesh.nd = master.dim;
mesh.plocal = master.xpe;
mesh.tlocal = master.telem;
mesh.porder = pde.porder;
ti = 2000;
% Uout = getsolution(['dataout/out'],dmd,master.npe);
Uout = getsolution(['dataout/out_t' num2str(ti)],dmd,master.npe);
UCGout = getsolution(['dataout/out_outputCG_t' num2str(ti)],dmd,master.npe);
% out = eulereval(Uout,'r',gam,Minf);
rho = sum(Uout(:,1:5,:),2);
out = Uout(:,1,:);
% out = rho;
% out = rho;
% outCG = 
[om,ind] = max(abs(out(:)));
figure(100); clf; scaplot(mesh,out,[],1);
x = mesh.dgnodes(:,1,:);
y = mesh.dgnodes(:,2,:);
hold on;
plot(x(ind),y(ind),'o');
colormap('jet');
%%
titles = ["w1","w2","w3","w4","w5","p","t","gam","rhoe","sb0"];
for i = 1:10
%     if i ~= 9
        figure(i); clf
        wcurr = UCGout(:,i,:);
        scaplot(mesh,wcurr,[],1);
        title(titles(i))
%         waitforbuttonpress
%     end
end
% figure(7); clf
P = UCGout(:,6,:);
T = UCGout(:,7,:);
gam = UCGout(:,8,:);
sb0 = UCGout(:,10,:);
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
