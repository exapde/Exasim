function mesh = createhighordermesh(mesh,app)

disp('run createnodes...');  
[mesh.dgnodes,mesh.elemtype,perm] = createnodes(mesh.p,mesh.t,app.porder);

disp('run facenumbering...');  
[mesh.f, mesh.tprd] = facenumbering(mesh.p,mesh.t,mesh.elemtype,mesh.bndexpr,mesh.periodicexpr);

% Project nodes on the curved boundaries
if ~isempty(mesh.curvedboundaryexpr) && app.porder>1 && max(mesh.curvedboundary)>0 
    disp('Project dgnodes onto the curved boundaries...');  
    fd = mesh.curvedboundaryexpr;
    nd=size(mesh.p,1);
    [nfe,ne] = size(mesh.f);  
    if nd == 2      
        for i = 1:ne
            for j = 1:nfe
                if mesh.f(j,i)~=0 % boundary element
                    k = abs(mesh.f(j,i)); % get boundary index
                    if mesh.curvedboundary(k)==1 % if this boundary is curved
                        p = mesh.dgnodes(perm(:,j),:,i);
                        deps = sqrt(eps)*max(max(p)-min(p));     
                        d = fd{k}(p);
                        dgradx = (fd{k}([p(:,1)+deps,p(:,2)])-d)/deps;
                        dgrady = (fd{k}([p(:,1),p(:,2)+deps])-d)/deps;
                        dgrad2 = dgradx.^2+dgrady.^2;
                        dgrad2(dgrad2==0) = 1;
                        p = p-[d.*dgradx./dgrad2,d.*dgrady./dgrad2];                  
                        mesh.dgnodes(perm(:,j),:,i) = p;                        
                    end
                end
            end
        end
    elseif nd==3
        for i = 1:ne
            for j = 1:nfe
                if mesh.f(j,i)~=0 % boundary element
                    k = abs(mesh.f(j,i)); % get boundary index
                    if mesh.curvedboundary(k)==1 % if this boundary is curved                
                        p = mesh.dgnodes(perm(:,j),:,i);
                        deps = sqrt(eps)*max(max(p)-min(p));     
                        d = fd{k}(p);
                        dgradx = (fd{k}([p(:,1)+deps,p(:,2),p(:,3)])-d)/deps;
                        dgrady = (fd{k}([p(:,1),p(:,2)+deps,p(:,3)])-d)/deps;
                        dgradz = (fd{k}([p(:,1),p(:,2),p(:,3)+deps])-d)/deps;
                        dgrad2 = dgradx.^2+dgrady.^2+dgradz.^2;
                        dgrad2(dgrad2==0) = 1;
                        p = p-[d.*dgradx./dgrad2,d.*dgrady./dgrad2,d.*dgradz./dgrad2];
                        mesh.dgnodes(perm(:,j),:,i) = p;         
                    end
                end
            end
        end        
    end        
end

end
