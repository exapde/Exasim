porder = 1;
AV = 0;
RBdim = 0;
dt = 2e-3;
meshfile = "./euler_mesh_390.msh";
tau = 1e-3;
buildpath=FILL ME IN

pdeapp_small(porder, AV, RBdim, dt, meshfile, tau, buildpath);
