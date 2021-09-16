% Add Exasim to Matlab search path
versiondir = cdir(1:(ii+5)) + "/"  + version + "/Matlab";
addpath(char(cdir(1:(ii+5)) + "/Installation"));
addpath(char(versiondir + "/Gencode"));
addpath(char(versiondir + "/Mesh"));
addpath(char(versiondir + "/Preprocessing"));
addpath(char(versiondir + "/Postprocessing"));
addpath(char(versiondir + "/Utilities"));

% Set Matlab's PATH enviroment variable so that Exasim can call external packages    
setenv('PATH', '/usr/local/bin:/usr/bin:/opt/local/bin:/bin:/usr/sbin:/sbin');
% Add more paths if neccesary 
setenv('PATH', [getenv('PATH') ':/Applications/ParaView-5.8.1.app/Contents/MacOS']);

fprintf("==> Exasim " + version + " ...\n"');
