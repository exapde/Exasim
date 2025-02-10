function [p_new] = dgnodes_to_p(mapping, mesh, master)
corner_dofs = master.perm(1 ,:);
mesh.p = mesh.p';
mesh.t = mesh.t';
% p_new = mesh.dgnodes(corner_dofs, :, :);

ne = size(mesh.dgnodes, 3);
nd = mesh.nd;
p_new = mesh.p;
for iel = 1:ne
    dg_el = mesh.dgnodes(corner_dofs,1:nd,iel);
    p_curr = mesh.p(mesh.t(iel,:),:);
    for idg = 1:size(dg_el,1)
        dg_curr = dg_el(idg,:);
        [tst, ii] = ismember(dg_curr, p_curr,'rows'); % find vertex that matches node
        if tst == 0
            disp("Something wrong");
        end
        p_match = mesh.t(iel, ii);
        p_new(p_match,:) = mapping(corner_dofs(idg), 1:nd, iel);
    end
end
end