% Add Exasim to Matlab search path
src = "frontends"; 
ExasimPath = cdir(1:(ii+5));
srcdir = ExasimPath + "/"  + src + "/Matlab";
addpath(char(ExasimPath + "/install"));
addpath(char(srcdir + "/Gencode"));
addpath(char(srcdir + "/master"));
addpath(char(srcdir + "/Mesh"));
addpath(char(srcdir + "/Mesh/mkmesh"));
addpath(char(srcdir + "/Mesh/cmesh"));
addpath(char(srcdir + "/Mesh/lesmesh"));
addpath(char(srcdir + "/Mesh/surfmesh"));
addpath(char(srcdir + "/HDG"));
addpath(char(srcdir + "/Preprocessing"));
addpath(char(srcdir + "/Postprocessing"));
addpath(char(srcdir + "/Utilities"));

% Set Matlab's PATH enviroment variable so that Exasim can call external packages    
setenv('PATH', '/usr/local/bin:/usr/bin:/opt/local/bin:/bin:/usr/sbin:/sbin:/home/linuxbrew/.linuxbrew/opt/llvm@11/bin/:/usr/local/cuda-11.7/bin/');
% Add more paths if neccesary 
setenv('PATH', [getenv('PATH') ':/Applications/ParaView-5.8.1.app/Contents/MacOS']);
setenv('LD_LIBRARY_PATH','/usr/bin:/usr/local/cuda-11.7/bin:/usr/local/cuda-11.7/lib64');

% create build folder if it does not exist
if exist(ExasimPath + "/build", "file") == 0
    mkdir(char(ExasimPath + "/build"));    
end

fprintf("==> Exasim " + " ...\n"');
