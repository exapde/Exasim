
function measure = computeElemMeasure(mesh,master)

nd = mesh.nd;
ne = mesh.ne;
ngv = master.ngv;
dgnodes = mesh.dgnodes;

shapmv = master.shapmv;
gwvl = master.gwvl;

[~, ~, jacE] = volgeom(shapmv,permute(dgnodes(:,1:nd,:),[1,3,2]));
jacE = reshape(jacE, [ngv,ne]);

measure = zeros(ne,1);
for elem = 1:ne
    measure(elem) = sum(gwvl.*jacE(:,elem));
end

end
