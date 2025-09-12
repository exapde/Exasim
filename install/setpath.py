import sys, os

# Add Exasim to Python search path
src = "frontends"; 
srcdir = cdir[0:(ii+6)] + "/"  + src + "/Python";
sys.path.append(cdir[0:(ii+6)] + '/Installation');
sys.path.append(srcdir + '/Gencode');
sys.path.append(srcdir + '/Mesh');
sys.path.append(srcdir + '/Preprocessing');
sys.path.append(srcdir + '/Postprocessing');
sys.path.append(srcdir + '/Utilities');
sys.path.append(cdir);

# Set Python's PATH enviroment variable so that Exasim can call external programs
os.environ['PATH'] = "/usr/local/bin:/usr/bin:/opt/local/bin:/bin:/usr/sbin:/sbin:/opt/homebrew/bin";
# Add more paths if neccesary
os.environ['PATH'] = os.environ['PATH'] + ":/Applications/ParaView-6.0.0.app/Contents/MacOS"

print('==> Exasim ' + src + ' ...\n');

