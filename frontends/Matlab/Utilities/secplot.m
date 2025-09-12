function secplot(mesh,master,mastersubgrid,UH,udg)

nqf = size(mesh.perm,1);
nfe = size(mesh.perm,2);
ne  = size(mesh.dgnodes,3);
nd  = mesh.nd;

UH  = reshape(UH(mesh.elcon),[nqf nfe ne]);
xi     = linspace(0,1,6)';
shapmf = mkshape(master.porder,master.plocfc,xi,1);
shapmf = shapmf(:,:,1)';

s = length(mastersubgrid);
elemnodes = mastersubgrid{s}.geomnodes;
ng = size(elemnodes,3);
for k=1:ng
    tmp = mkshape(mesh.porder,mesh.plocal,elemnodes(:,:,k),mesh.elemtype);
    shp(:,:,k) = tmp(:,:,1)';
    tmp = mkshape(mastersubgrid{s}.porder,mastersubgrid{s}.plocvl,elemnodes(:,:,k),mesh.elemtype);
    shv(:,:,k) = tmp(:,:,1)';
end            

figure(1); clf;
hold on;
for k = 1:ne    
    for j=1:nfe
        dx = mesh.dgnodes(:,:,k);
        dg = dx(mesh.perm(:,j),:);        
        if max(abs(dg(:,2)))<1e-4 
            dg = shapmf*reshape(dg,[nqf nd]);        
            uh = shapmf*UH(:,j,k);
            plot(dg(:,1),uh,'-k','Linewidth',1);
            
            for q=1:ng 
                xg = shp(:,:,q)*dx;                   
                ug = shv(:,:,q)*udg{s}(:,1,(k-1)*ng+q);                   
                for r=1:nfe
                    dy = xg(mesh.perm(:,r),:);                    
                    if max(abs(dy(:,2)))<1e-4 && min(dx(:,2))<-1e-4                         
                        plot(dy(:,1),ug(mesh.perm(:,r)),'-r','Linewidth',1);                       
                    end
%                     if max(abs(dy(:,2)))<1e-4 && max(dx(:,2))>1e-4 
%                         plot(dy(:,1),ug(mesh.perm(:,r)),'-b','Linewidth',1);                       
%                     end
                end
            end            
        end
    end
end
hold off;
%axis equal;
axis tight;
xlabel('x','Fontsize',16);
ylabel('Pressure profile at y = 0','FontSize',16);
legend({'uhat','uplus','uminus'},'FontSize',14,'Location','NorthWest');
box on;

% for j = 1:ns % for each subgrid
%     % obtain element indices for subgrid j
%     ind = find(mesh.subgrids==j);    
%     if isempty(ind)==0
%         nj = length(ind);
%         ng = size(mastersubgrid{j}.elemnodes,3);
%         npm = size(mesh.dgnodes,1);        
%         npv = size(u{j},1);
%         dgnodes = mesh.dgnodes(:,:,ind);                                        
% %         udg = subgridprojection(mastersubgrid{j},mastersubgrid{1},mesht.dgnodes,reshape(u{j}(:,iu,:),[npv 1 ng nj]));                
% %         udg = reshape(udg,[npm 1 nj]);        
% %         scalarplot(mesht,udg,nref);        
%         
%         % obtain subgrid interpolation nodes
% %         if mastersubgrid{j}.porder==0
% %             [plocvl,tlocvl] = masternodes(1,nd,mastersubgrid{j}.elemtype,mastersubgrid{j}.nodetype);
% %             elemnodes = mknodes(mastersubgrid{j}.p,mastersubgrid{j}.t,plocvl);            
% %             mesht.plocal  = plocvl;
% %             mesht.tlocal  = tlocvl;                
% %         else
% %             elemnodes = mastersubgrid{j}.elemnodes;    
% %             mesht.plocal  = mastersubgrid{j}.plocvl;
% %             mesht.tlocal  = mastersubgrid{j}.tlocvl;                
% %         end
%         elemnodes = mastersubgrid{j}.geomnodes;
%         ng = size(elemnodes,3);
% 
%         % compute shape functions at the subgrid DG nodes
% %         tmp = mkshape(mesh.porder,mesh.plocal,elemnodes(:,:,1),mesh.elemtype);
% %         n1 = size(tmp,1); 
% %         n2 = size(tmp,2); 
% %         shp = zeros(npm,n1,ng);                
% %         shv = zeros(n2,n0,ng);
%         clear shp shv;
%         for k=1:ng
%             tmp = mkshape(mesh.porder,mesh.plocal,elemnodes(:,:,k),mesh.elemtype);
%             shp(:,:,k) = tmp(:,:,1)';
%             tmp = mkshape(mastersubgrid{j}.porder,mastersubgrid{j}.plocvl,elemnodes(:,:,k),mesh.elemtype);
%             shv(:,:,k) = tmp(:,:,1)';
%         end            
%                 
%         % compute subgrid DG nodes in the physical space
%         %nc  = size(u{j},2);        
%         npm = size(dgnodes,1);        
%         dgx = zeros(npm,nd,nj*ng);        
%         udg = zeros(npm,1,nj*ng);        
%         for i=1:nj % for each superelement                           
%             for k=1:ng % for each subgrid element                    
%                 dgx(:,:,(i-1)*ng+k) = shp(:,:,k)*dgnodes(:,:,i);                   
%                 udg(:,1,(i-1)*ng+k) = shv(:,:,k)*u{j}(:,iu,(i-1)*ng+k);                   
%             end
%         end        
%         
%         mesht.dgnodes = dgx; 
%         mesht.t       = mesh.t(ind,:);
%         %mesht.porder  = mastersubgrid{j}.porder;   
%         %npv = size(u,1);
%         %nc  = size(u,2);
%         nes = size(mastersubgrid{j}.t,1);        
%         %scalarplot(mesht,u{j}(:,:,(ind(1)-1)*nes+1:ind(end)*nes),nref)        
% %         size(dgx)
% %         size(u{j})
%         scalarplot(mesht,udg,nref);
%         hold on;
%     end
% end
