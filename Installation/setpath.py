import sys, os

# Add Exasim to Python search path
versiondir = cdir[0:(ii+6)] + "/"  + version + "/Python";
sys.path.append(versiondir + '/Gencode');
sys.path.append(versiondir + '/Mesh');
sys.path.append(versiondir + '/Preprocessing');
sys.path.append(versiondir + '/Postprocessing');
sys.path.append(versiondir + '/Utilities');
sys.path.append(cdir);

# Set Python's PATH enviroment variable so that Exasim can call external programs
os.environ['PATH'] = "/usr/local/bin:/usr/bin:/opt/local/bin:/bin:/usr/sbin:/sbin";
# Add more paths if neccesary
os.environ['PATH'] = os.environ['PATH'] + ":/Applications/ParaView-5.8.1.app/Contents/MacOS"

print('==> Exasim ' + version + ' ...\n');
