hybrid='hdg';
preconditioner = 0;
ngrid=20;
porder = 4;
nproc = 4;
elemtype=1;
nodetype=1;

% HDG
mesh = mkmesh_square(ngrid,ngrid,porder,1,1,1,elemtype,nodetype);
master = mkmaster(mesh,2*mesh.porder);
[master,mesh] = preprocess(master,mesh,hybrid);
t2t = mkt2t(mesh.t,mesh.elemtype);
[epart, npart, ent2ent, ent2entStart] = meshpart(mesh,nproc,mesh.p,mesh.t);

if strcmp(hybrid,'hdg')
    dmdmat = domaindecompositionmpi(mesh.t2f',t2t,mesh.t2f,nproc,preconditioner,mesh.elcon);    
else
    dmdmat = domaindecompositionmpi(mesh.elcon,t2t,mesh.t2f,nproc,preconditioner);
end

param = zeros(20,1);
param(1) = mesh.ne;
param(2) = mesh.nf;
if strcmp(hybrid, 'hdg')
    param(3) = mesh.nf;    
    param(13) = 1;
else
    param(3) = max(mesh.elcon(:));    
    param(13) = 0;
end
param(4) = size(mesh.perm,2);
param(5) = size(t2t,2);
param(6) = size(mesh.dgnodes,1);
param(7) = size(mesh.perm,1);
param(11) = preconditioner;
param(12) = preconditioner+2;
param(14) = nproc;

fileID = fopen('dmdin.bin','w');
fwrite(fileID,param,'double');
fwrite(fileID,mesh.elcon-1,'double');
fwrite(fileID,t2t-1,'double');
fwrite(fileID,mesh.t2f-1,'double');
fwrite(fileID,epart,'double');
fwrite(fileID,npart,'double');
fclose(fileID);

return;

for i = 1:nproc
    disp(['---------------------------']);
    disp(['CPU:  ' num2str(i)]);
    disp(['---------------------------']);
    fn = ['dmdout' num2str(i-1)];
    dmdcpp{i} = readDMDSTRUCTfromBinaryFile(fn);
    validate_dmd(dmdmat{i},dmdcpp{i});
end

