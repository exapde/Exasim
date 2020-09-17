x = mesh.dgnodes(:,1,:);
y = mesh.dgnodes(:,2,:);
UDG(:,1,:) = sin(pi*x).*sin(pi*y); % exact solution
UDG(:,2,:) = -(pi)*cos(pi*x).*sin(pi*y);
UDG(:,3,:) = -(pi)*sin(pi*x).*cos(pi*y);
filename = ['data/' app.appname 'out'];
tmp = getsolution(filename, mpiprocs, size(UDG,1), app.nc);
fprintf('Maximum error in UDG: %g\n',max(abs(tmp(:)-UDG(:))));

[cgnodes,cgelcon,rowent2elem,colent2elem,cgent2dgent] = mkcgent2dgent(mesh.dgnodes,1e-6);

filename = ['data/' app.appname 'out_outputCG'];
outdg = getsolution(filename, mpiprocs, size(UDG,1), app.nce);
e1 = outdg(:,1,:)-tmp(:,1,:);
e2 = outdg(:,2,:)-tmp(:,2,:)-tmp(:,3,:);
[max(abs(e1(:))) max(abs(e2(:)))]

% filename = ['data/' app.appname 'out_outputCG'];
% fileID = fopen([filename '_np0.bin'],'r');
% outcg = fread(fileID,'double');
% [ucg,vcg] = dg2cg(outdg, cgelcon, cgent2dgent, rowent2elem);
% max(abs(vcg(:)-outcg(:)))
% max(abs(ucg(:)-outdg(:)))

% filename = ['data/' app.appname 'out_output'];
% out = getsolution(filename, mpiprocs, size(UDG,1), app.nce);
% %e1 = out(:,1,:)-tmp(:,1,:);
% e1 = out(:,1,:)-tmp(:,2,:)-tmp(:,3,:);
% max(abs(e1(:)))
% 
% [cgnodes,cgelcon,rowent2elem,colent2elem,cgent2dgent] = mkcgent2dgent(mesh.dgnodes,1e-6);
% ucg = dg2cg(out, cgelcon, cgent2dgent, rowent2elem);
% %e1 = ucg(:,1,:)-tmp(:,1,:);
% e1 = ucg(:,1,:)-tmp(:,2,:)-tmp(:,3,:);
% max(abs(e1(:)))


% figure(1); clf; scaplot(mesh,UDG(:,1,:),[],2,1); axis off; colormap jet;
% figure(2); clf; scaplot(mesh,UDG(:,2,:),[],2,1); axis off; colormap jet;
% figure(3); clf; scaplot(mesh,UDG(:,3,:),[],2,1); axis off; colormap jet;
% master=mkmasterelement(mesh.nd,mesh.porder,mesh.porder,3*mesh.porder,3*mesh.porder,mesh.elemtype,mesh.nodetype);
% qdg = gradu(master.shapnt(:,:,2:end), mesh.dgnodes, UDG(:,1,:));
% figure(4); clf; scaplot(mesh,qdg(:,1,:),[],2,1); axis off; colormap jet;
% figure(5); clf; scaplot(mesh,qdg(:,2,:),[],2,1); axis off; colormap jet;
% [cgnodes,cgelcon,rowent2elem,colent2elem,cgent2dgent] = mkcgent2dgent(mesh.dgnodes,1e-6);
% ucg = dg2cg(tmp, cgelcon, cgent2dgent, rowent2elem);
%
% figure(2); clf; scaplot(mesh,ucg(:,1,:),[],2,1); axis off; colormap jet;
%
% fileID = fopen([filename 'ucg.bin'],'r');
% tm1 = fread(fileID,'double');
% fclose(fileID);
% tm2 = ucg(:,1,:);
% max(abs(tm1(:)-tm2(:)))
%
% fileID = fopen([filename 'udg.bin'],'r');
% tm1 = fread(fileID,'double');
% fclose(fileID);
% tm2 = tmp(:,1,:);
% max(abs(tm1(:)-tm2(:)))
