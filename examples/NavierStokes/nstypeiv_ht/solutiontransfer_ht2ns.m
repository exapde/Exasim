function [vdg, UHb, in, im] = solutiontransfer_ht2ns(pde, dmd, mesh, meshns, UHht)

% Get UDG and UH from the binary files
% UDGht = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
if nargin == 4
fileID = fopen(pde.buildpath + "/dataout/out_uhat_np0.bin",'r');
UHht = fread(fileID,'double');
fclose(fileID);
UHht = reshape(UHht,pde.ncu,[]);    
end

% Get UDG and UH on the interface
%UDGb = getsolutiononboundary(UDGht, mesh.f, mesh.perm, mesh.ibwall);
UHb = getuhonboundary(UHht, dmd{1}.elemcon, mesh.f, mesh.perm, mesh.ibwall);

% get the temprature on the interface
Ubht = 1.0*UHb(:,1,:);

% match the DG nodes on the interface
XDGbht = getsolutiononboundary(mesh.dgnodes, mesh.f, mesh.perm, mesh.ibwall);
XDGbns = getsolutiononboundary(meshns.dgnodes, meshns.f, meshns.perm, meshns.ibwall);
[in,im] = matchdgnodesonboundary(XDGbht, XDGbns);

% transfer the temprature from the heat domain to the fluid domain
Ubns = 0*Ubht;
for i = 1:length(in)
  i2 = in(i);
  Ubns(im(:,i),:,i2) = Ubht(:,:,i);
end

% store the temperature into vdg
vdg = zeros(size(meshns.dgnodes,1), size(Ubns,2) ,size(meshns.dgnodes,3));
vdg = putsolutiononboundary(vdg, Ubns, meshns.f, meshns.perm, meshns.ibwall);


