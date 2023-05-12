% Add Exasim to Matlab search path
src = "src"; 
srcdir = cdir(1:(ii+5)) + "/"  + src + "/Matlab";
addpath(char(cdir(1:(ii+5)) + "/Installation"));
addpath(char(srcdir + "/Gencode"));
addpath(char(srcdir + "/Mesh"));
addpath(char(srcdir + "/Mesh/lesmesh"));
addpath(char(srcdir + "/Preprocessing"));
addpath(char(srcdir + "/Postprocessing"));
addpath(char(srcdir + "/Utilities"));
addpath(char(srcdir + "/Utilities/rom"));

% Set Matlab's PATH enviroment variable so that Exasim can call external packages    
setenv('PATH', '/usr/local/bin:/usr/bin:/opt/local/bin:/bin:/usr/sbin:/sbin');
% Add more paths if neccesary 
setenv('PATH', [getenv('PATH') ':/Applications/ParaView-5.8.1.app/Contents/MacOS']);

fprintf("==> Exasim " + " ...\n"');
