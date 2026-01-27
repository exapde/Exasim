% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

pde = initializeexasim();
pde.builtinmodelID = 0;
pde.modelfile = "pdemodel" + num2str(pde.builtinmodelID);  
pde.hybrid = 1;
pde.nd = 0;
pde.ncu = 0;
pde.ncq = 0;
pde.ncw = 0;
pde.ncv = 0;
pde.ntau = 0;
pde.nmu = 0;
pde.neta = 0;

kkgenmodel(pde);

