function [vdg, UHb, in, im] = solutiontransfer_ns2ht(pde, dmd, mesh, meshht)

% Get UDG and UH from the binary files
fileID = fopen(pde.buildpath + "/dataout/out_np0.bin",'r');
UDGns = fread(fileID,'double');
fclose(fileID);
npe = size(mesh.dgnodes,1);
ne = size(mesh.dgnodes,3);
nc = numel(UDGns)/(npe*ne);
UDGns = reshape(UDGns,npe,nc,ne);

fileID = fopen(pde.buildpath + "/dataout/out_uhat_np0.bin",'r');
UHns = fread(fileID,'double');
fclose(fileID);
UHns = reshape(UHns,pde.ncu,[]);    

% Get UDG and UH on the interface
UDGb = getsolutiononboundary(UDGns, mesh.f, mesh.perm, mesh.ibwall);
UHb = getuhonboundary(UHns, dmd{1}.elemcon, mesh.f, mesh.perm, mesh.ibwall);

% get the heat flux on the interface
[hfbx, hfby] = heatflux(UHb, UDGb(:,(pde.ncu+1):end,:), pde.physicsparam);
Ubns = 1.0*UHb;
Ubns(:,1,:) = UDGb(:,4,:);
Ubns(:,2,:) = hfbx;
Ubns(:,3,:) = hfby;
Ubns(:,4,:) = UHb(:,4,:);

% match the DG nodes on the interface
XDGbns = getsolutiononboundary(mesh.dgnodes, mesh.f, mesh.perm, mesh.ibwall);
XDGbht = getsolutiononboundary(meshht.dgnodes, meshht.f, meshht.perm, meshht.ibwall);
[in,im] = matchdgnodesonboundary(XDGbns, XDGbht, 1e-8);

% transfer the heat flux from the NS domain to the heat equation domain
Ubht = 0*Ubns;
for i = 1:length(in)
  i2 = in(i);
  Ubht(im(:,i),:,i2) = Ubns(:,:,i);
end

% store the heat flux into vdg
vdg = zeros(size(meshht.dgnodes,1), size(Ubht,2) ,size(meshht.dgnodes,3));
vdg = putsolutiononboundary(vdg, Ubht, meshht.f, meshht.perm, meshht.ibwall);


