% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

ExasimPath = cdir(1:(ii+5));
srcdir = ExasimPath + "/"  + src + "/Matlab";
addpath(char(srcdir + "/Modeling/CNS5air/"));

pde = initializeexasim();
pde.builtinmodelID = 4;
pde.modelfile = "pdemodel" + num2str(pde.builtinmodelID);  
pde.hybrid = 1;
pde.nd = 2;
pde.ncu = 8;
pde.ncq = 16;
pde.ncw = 1;
pde.ncv = 1;
pde.ntau = 1;
pde.nmu = 12;
pde.neta = 18;

kkgenmodel(pde);

