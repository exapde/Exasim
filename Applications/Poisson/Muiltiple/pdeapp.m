pdeapp1;
pdeapp2;

gencodeall(2);

% compile source codes to build an executable file and store it in app folder
compilerstr = compilecode(pde1);

% run executable file to compute solution and store it in dataout folder
runstr = runcode(pde1,2);
