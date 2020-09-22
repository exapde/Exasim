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

    if app['porder']>1:
        visorder = min(2*app['porder'],8);
    else:
        visorder = app['porder'];

    mesh['xpe'] = Preprocessing.masternodes(app['porder'],app['nd'],app['elemtype'])[0];
    xpe,telem = Preprocessing.masternodes(visorder,app['nd'],app['elemtype'])[0:2];

    visshape = Preprocessing.mkshape(app['porder'],mesh['xpe'],xpe,app['elemtype']);
    visshape = visshape[:,:,0].T;

    dgnodes = Preprocessing.createdgnodes(mesh['p'],mesh['t'][:,app['viselem']],mesh['f'][:,app['viselem']],mesh['curvedboundary'],mesh['curvedboundaryexpr'],visorder);
    cgnodes, cgelcon, cgcells, celltype = createcggrid(dgnodes,telem)[0:4];
    dgnodes = Preprocessing.createdgnodes(mesh['p'],mesh['t'][:,app['viselem']],mesh['f'][:,app['viselem']],mesh['curvedboundary'],mesh['curvedboundaryexpr'],app['porder']);

    app['paraview'] = Preprocessing.findexec(app['paraview'],app['version']);

    if nt==1:
        tm = matmul(visshape,reshape(visfields[:,:,app['viselem']],(visshape.shape[1], visfields.shape[1]*len(app['viselem'])), 'F'));
        tm = reshape(tm,(visshape.shape[0], visfields.shape[1], len(app['viselem'])),'F');

        vtuwrite(app['visfilename'], cgnodes, cgelcon, cgcells, celltype, app['visscalars'], app['visvectors'], tm);
        if len(app['paraview'])>0:
            str = app['paraview'] + " --data=" + app['visfilename'] + ".vtu &";
            os.system(str);
    else:
        pvdwrite(app['visfilename'], cgnodes, cgelcon, cgcells, celltype, app['visscalars'], app['visvectors'], visfields[:,:,app['viselem'],:],app['visdt'],visshape);
        if len(app['paraview'])>0:
            str = app['paraview'] + " --data=" + app['visfilename'] + ".pvd &";
            os.system(str);

    return dgnodes;
