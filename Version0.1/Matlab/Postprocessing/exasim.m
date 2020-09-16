function [sol,pde,mesh,master,dmd,compilerstr,runstr] = exasim(pde,mesh)

% search compilers and set options
pde = setcompilers(pde);       

% generate input files and store them in datain folder
[pde,mesh,master,dmd] = preprocessing(pde,mesh);

% generate source codes and store them in app folder
gencode(pde);

% compile source codes to build an executable file and store it in app folder
compilerstr = compilecode(pde);

% run executable file to compute solution and store it in dataout folder
runstr = runcode(pde);

% get solution from output files in dataout folder
sol = fetchsolution(pde,master,dmd);


