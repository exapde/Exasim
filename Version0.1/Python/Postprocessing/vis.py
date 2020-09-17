import os, sys, shutil
from numpy import *
import Preprocessing
from createcggrid import createcggrid
from vtuwrite import vtuwrite
from pvdwrite import pvdwrite

def vis(visfields,app,mesh):

    nt = 0;
    if visfields.ndim==3:
        nt = 1;
    else:
        nt = visfields.shape[3];

    if app['viselem'] == []:
        ne = mesh['t'].shape[1];
        app['viselem'] = range(0,ne);

    dgnodes = Preprocessing.createdgnodes(mesh['p'],mesh['t'][:,app['viselem']],mesh['f'][:,app['viselem']],mesh['curvedboundary'],mesh['curvedboundaryexpr'],app['porder']);
    cgnodes, cgelcon, cgcells, celltype = createcggrid(dgnodes,mesh['telem'])[0:4];

    app['paraview'] = Preprocessing.findexec(app['paraview'],app['version']);
    # paraview = app['paraview'];
    # paraviewstatus0 = shutil.which(paraview);
    # paraviewstatus1 = shutil.which("paraview");
    # paraviewstatus2 = shutil.which("usr/bin/paraview");
    # paraviewstatus3 = shutil.which("/usr/local/bin/paraview");
    # paraviewstatus4 = shutil.which("/opt/local/bin/paraview");
    #
    # if paraviewstatus0 != None:
    #     paraview = paraview;
    # elif paraviewstatus1 != None:
    #     paraview = "paraview"
    # elif paraviewstatus2 != None:
    #     paraview = "/usr/bin/paraview";
    # elif paraviewstatus3 != None:
    #     paraview = "/usr/local/bin/paraview";
    # elif paraviewstatus4 != None:
    #     paraview = "/opt/local/bin/paraview";
    # else:
    #     sys.exit("Exasim search in /usr/bin, /usr/local/bin, and /opt/local/bin and could not find Paraview. Please see the documentation to install it. After installation, please set its path to app['paraview']");
    # app['paraview'] = paraview;

    if nt==1:
        vtuwrite(app['visfilename'], cgnodes, cgelcon, cgcells, celltype, app['visscalars'], app['visvectors'], visfields[:,:,app['viselem']]);
        if len(app['paraview'])>0:
            str = app['paraview'] + " --data=" + app['visfilename'] + ".vtu &";
            os.system(str);
    else:
        pvdwrite(app['visfilename'], cgnodes, cgelcon, cgcells, celltype, app['visscalars'], app['visvectors'], visfields[:,:,app['viselem'],:],app['visdt']);
        if len(app['paraview'])>0:
            str = app['paraview'] + " --data=" + app['visfilename'] + ".pvd &";
            os.system(str);

    return dgnodes;
