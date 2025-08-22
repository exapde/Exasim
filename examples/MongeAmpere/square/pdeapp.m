% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

[mesh, rho, drhodx, drhody] = mkmesh_square(50, 4, 1);

mesh1 = radaptivity(mesh, rho, drhodx, drhody, 1e-4);


