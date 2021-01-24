function uhplot(mesh,master,UH)

nqf = size(mesh.perm,1);
nfe = size(mesh.perm,2);
ne  = size(mesh.dgnodes,3);
nd  = mesh.nd;

UH  = reshape(UH(mesh.elcon),[nqf nfe ne]);

xi     = linspace(0,1,40)';
shapmf = mkshape(master.porder,master.plocfc,xi,1);
shapmf = shapmf(:,:,1)';

figure(1); clf;
hold on;
for k = 1:ne
    uh = shapmf*UH(:,:,k);
    dg = mesh.dgnodes(mesh.perm,:,k);
    dg = shapmf*reshape(dg,[nqf nfe*nd]);
    dg = reshape(dg,[],nd);
    plot3(dg(:,1),dg(:,2),uh(:));
end
hold off;
axis equal;
axis tight;
view(3);

