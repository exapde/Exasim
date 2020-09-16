% Add Exasim to Matlab search path
versiondir = cdir(1:(ii+5)) + "/"  + version + "/Matlab";
addpath(versiondir + '/Gencode');
addpath(versiondir + '/Mesh');
addpath(versiondir + '/Preprocessing');
addpath(versiondir + '/Postprocessing');
addpath(versiondir + '/Utilities');

% Set Matlab's PATH enviroment variable so that Exasim can call external packages    
setenv('PATH', '/usr/local/bin:/usr/bin:/opt/local/bin:/bin:/usr/sbin:/sbin');
% Add more paths if neccesary 
setenv('PATH', [getenv('PATH') ':/Applications/ParaView-5.8.1.app/Contents/MacOS']);

fprintf("==> Exasim " + version + " ...\n"');
