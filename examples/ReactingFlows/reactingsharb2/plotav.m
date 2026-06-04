% ----- plot current AV field from the initial solution -----

ns  = 5;
nd  = 2;
ncu = ns + nd + 1;   % 8
nq  = ncu*nd;        % 16

% Check that mesh.udg contains both u and q
assert(size(mesh.udg,2) >= ncu+nq, ...
    'mesh.udg does not contain enough columns for [u,q].');

u0 = mesh.udg(:,1:ncu,:);
q0 = mesh.udg(:,ncu+1:ncu+nq,:);

hm_current = mesh.vdg(:,1,:);

avcoeff = 1e-4; %pde.physicsparam(13);

npv = size(mesh.dgnodes,1);
ne  = size(mesh.dgnodes,3);

av0_current = zeros(npv,1,ne);

for e = 1:ne
    for i = 1:npv
        y   = mesh.dgnodes(i,2,e);
        ui  = squeeze(u0(i,:,e)).';
        qi  = squeeze(q0(i,:,e)).';
        hmi = hm_current(i,1,e);

        av0_current(i,1,e) = getavfield2dchem(y, ui, qi, hmi, [], avcoeff, pde.porder);
    end
end

fprintf('Current AV min = %.3e, max = %.3e\n', min(av0_current(:)), max(av0_current(:)));

figure;
scaplot(mesh, av0_current, [-1e3 1e3], 1);
axis equal tight;
colorbar; colormap(jet);
title('Initial AV field using current transient input v(9)');