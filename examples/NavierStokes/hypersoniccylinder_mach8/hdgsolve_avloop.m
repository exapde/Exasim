function [UDG1, UH1, ACG1, mine1, minf1, ming1, mesh1] = hdgsolve_avloop(master, master_mat, mesh, mesh_mat, app, q, UDG, UH, lambda0, kappa0)

mesh1 = mesh;
% mesh1.ib = [];
% mesh1.in = 1:size(mesh1.p2,1);
mesh1.dgnodes(:,1:2,:) = q;
% mesh1.dgnodes(:,3,:) = 0;
mesh1.dist = tanh(mesh.dist*30);
% mesh1.dgnodes(:,3,:) = 0.25*mesh1.dist;

S0 = 0.2;
eta = 0.9; m = 9;
lambda = ones(m,1)*lambda0;
for i = 2:m
    lambda(i) = lambda(i-1)*eta;
end
kappa = ones(m,1)*kappa0;
for i = 2:m
    kappa(i) = 1 + (kappa(i-1)-1)*eta;
end

[UDG1, UH1, ACG1, mine1, minf1, ming1] = avloop(master, master_mat, mesh1, mesh_mat, app, UDG, UH, S0, lambda, kappa,0);
