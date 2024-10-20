function mesh=mkcgmesh(mesh)

% remove duplicate nodes in mesh.p1
dim = size(mesh.p,2);
[ns,dim,nt]=size(mesh.dgnodes(:,1:dim,:));
A=reshape(permute(mesh.dgnodes(:,1:dim,:),[1,3,2]),[ns*nt,dim]);
snap=1e-8;
A=round(A/snap)*snap;
B=flipdim(A,1);
[~,I]=unique(B,'rows'); 
B=B(sort(I),:);
B=flipdim(B,1);
[~,b]=ismember(A,B,'rows');

% CG mesh
mesh.p2=B;
mesh.t2=reshape(b,[ns nt])';
