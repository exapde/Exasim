from initializeexasim import initializeexasim
from findexec import findexec
from masternodes import masternodes
from mkshape import mkshape
from mkmaster import mkmaster
from writeapp import writeapp
from readapp import readapp
from checkapp import checkapp
from writesol import writesol
from writebin import writebin
from writemaster import writemaster
from readmeshstruct import readmeshstruct
from checkmesh import checkmesh
from meshpartition2 import meshpartition2
from mkcgent2dgent import mkcgent2dgent
from mkelemblocks import mkelemblocks
from mkfaceblocks import mkfaceblocks
from mkdge2dgf import mkdge2dgf
from createdgnodes import createdgnodes
from facenumbering import facenumbering
import os, sys, sympy
from importlib import import_module
import importlib.util
from numpy import *

def preprocessing(app,mesh):

    if app['modelnumber']==0:
        strn = "";
    else:
        strn = str(app['modelnumber']);    
     
    if os.path.isdir("datain" + strn)==False:
        os.mkdir("datain" + strn)
    if os.path.isdir("dataout" + strn)==False:
        os.mkdir("dataout" + strn)

    filename = "datain" + strn + "/";
    fileapp = filename + "app.bin";
    filemaster = filename + "master.bin";

    if app['preprocessmode']==0:
        # update app structure
        app = writeapp(app,fileapp);
        return app;

    app['nd']  = mesh['p'].shape[0];
    app['ncx'] = app['nd'];

    nve,ne = mesh['t'].shape;

    app['elemtype'] = 0;
    if (app['nd']==2) and (nve==4):
        app['elemtype'] = 1;
    if (app['nd']==3) and (nve==8):
        app['elemtype'] = 1;
    app['pgauss'] = 2*app['porder'];

    master = mkmaster(app['nd'],app['porder'],app['pgauss'],app['elemtype'],app['nodetype']);
    writemaster(master,filemaster);

    pdemodel = import_module(app['modelfile']);
    # spec = importlib.util.spec_from_file_location('pdemodel', app['modelfile'])
    # pdemodel = importlib.util.module_from_spec(spec);
    # spec.loader.exec_module(pdemodel);

    app['boundaryconditions'] = mesh['boundarycondition'];
    app['uinf'] = app['externalparam'];

    nuinf = len(app['uinf']);
    nparam = len(app['physicsparam']);
    xdgsym = sympy.symbols("xdg1:" + str(app['ncx']+1));
    uinfsym = sympy.symbols("uinf1:" + str(nuinf+1));
    paramsym = sympy.symbols("param1:" + str(nparam+1));
    if hasattr(pdemodel, 'initu'):
        udgsym = pdemodel.initu(xdgsym,paramsym,uinfsym);
        app['ncu'] = len(udgsym);
    else:
        sys.exit('pdemodel.initu is not defined')
    if hasattr(pdemodel, 'initv'):
        vdgsym = pdemodel.initv(xdgsym,paramsym,uinfsym);
        app['nco'] = len(vdgsym);
    else:
        app['nco'] = 0;
    if hasattr(pdemodel, 'initw'):
        wdgsym = pdemodel.initw(xdgsym,paramsym,uinfsym);
        app['ncw'] = len(wdgsym);
    else:
        app['ncw'] = 0;

    if app['model']=="ModelC" or app['model']=="modelC":
        app['wave'] = 0;
        app['nc'] = app['ncu'];
    elif app['model']=="ModelD" or app['model']=="modelD":
        app['wave'] = 0;
        app['nc'] = round((app['ncu'])*(app['nd']+1));
    elif app['model']=="ModelW" or app['model']=="modelW":
        app['tdep'] = 1;
        app['wave'] = 1;
        app['nc'] = round((app['ncu'])*(app['nd']+1));

    app['ncq'] = app['nc'] - app['ncu'];
    app['nch'] = app['ncu'];

    if max(app['dt'])>0.0:
        app['tdep'] = 1;
    else:
        app['tdep'] = 0;

    print("run facenumbering...");
    mesh['f'], mesh['tprd'], t2t = facenumbering(mesh['p'],mesh['t'],app['elemtype'],mesh['boundaryexpr'],mesh['periodicexpr'])[0:3];

    mpiprocs = app['mpiprocs'];
    dmd = [dict() for x in range(mpiprocs)];
    #dmd = meshpartition(dmd,mesh['p'],mesh['t'],mesh['f'],mesh['tprd'],app['elemtype'],app['boundaryconditions'],mesh['boundaryexpr'],mesh['periodicexpr'],app['porder'],mpiprocs,app['metis']);
    dmd = meshpartition2(dmd,mesh['tprd'],mesh['f'],t2t,app['boundaryconditions'],app['nd'],app['elemtype'],app['porder'],mpiprocs,app['metis']);

    for i in range(0,mpiprocs):
        ii = i + 1;
        xdg = createdgnodes(mesh['p'],mesh['t'][:,dmd[i]['elempart']],mesh['f'][:,dmd[i]['elempart']],mesh['curvedboundary'],mesh['curvedboundaryexpr'],app['porder']);

        cgelcon,rowent2elem,colent2elem,cgent2dgent = mkcgent2dgent(xdg,1e-8)[0:4];

        if mpiprocs==1:
            writesol(filename + "/sol",0,xdg);
            ne = size(dmd[i]['elempart']);
            eblks,nbe = mkelemblocks(ne,app['neb'])[0:2];
            eblks = concatenate([eblks, zeros((1, eblks.shape[1]))]);
            mf = cumsum(concatenate([[0], dmd[i]['facepartpts'].flatten()]));
            fblks,nbf = mkfaceblocks(mf,dmd[i]['facepartbnd'],app['nfb']);
            neb = max(eblks[1,:]-eblks[0,:])+1;
            nfb = max(fblks[1,:]-fblks[0,:])+1;
        else:
            writesol(filename + "/sol",i+1,xdg);
            me = cumsum([0, dmd[i]['elempartpts'][0], dmd[i]['elempartpts'][1], dmd[i]['elempartpts'][2]]);
            eblks,nbe = mkfaceblocks(me,[0, 1, 2],app['neb'])[0:2];
            mf = cumsum(concatenate([[0], dmd[i]['facepartpts'].flatten('F')]));
            fblks,nbf = mkfaceblocks(mf,dmd[i]['facepartbnd'],app['nfb']);
            neb = max(eblks[1,:]-eblks[0,:])+1;
            nfb = max(fblks[1,:]-fblks[0,:])+1;

        npe = master['npe'];
        nfe = master['perm'].shape[1];
        facecon1 = reshape(dmd[i]['facecon'][:,0,:], (dmd[i]['facecon'].shape[0], dmd[i]['facecon'].shape[2]), 'F');
        facecon2 = reshape(dmd[i]['facecon'][:,1,:], (dmd[i]['facecon'].shape[0], dmd[i]['facecon'].shape[2]), 'F');
        ind = [];
        for ii in range (0,fblks.shape[1]):
            if fblks[2,ii]>0:
                ind = concatenate([ind, arange(fblks[0,ii]-1,fblks[1,ii])]);

        if len(ind)>0:
            facecon2 = facecon2[:, setdiff1d(arange(0,facecon1.shape[1]), ind).flatten('F')];
        rowe2f1,cole2f1,ent2ind1 = mkdge2dgf(facecon1,master['npe']*size(dmd[i]['elempart']))[0:3];
        rowe2f2,cole2f2,ent2ind2 = mkdge2dgf(facecon2,master['npe']*size(dmd[i]['elempart']))[0:3];

        dmd[i]['elempart'] = dmd[i]['elempart'] + 1;
        dmd[i]['facecon'] = dmd[i]['facecon'] + 1;
        if mpiprocs>1:
            dmd[i]['nbsd'] = dmd[i]['nbsd'] + 1;
            dmd[i]['elemrecv'] = dmd[i]['elemrecv'] + 1;
            dmd[i]['elemsend'] = dmd[i]['elemsend'] + 1;

        ndims = zeros(20);
        ndims[1-1] = mesh['p'].shape[0];
        ndims[2-1] = size(dmd[i]['elempart']);
        ndims[3-1] = sum(dmd[i]['facepartpts']);
        ndims[4-1] = max(mesh['t'].flatten())+1;
        ndims[5-1] = nfe;
        ndims[6-1] = nbe;
        ndims[7-1] = neb;
        ndims[8-1] = nbf;
        ndims[9-1] = nfb;

        nsize = zeros((21,1));
        nsize[1-1] = size(ndims);
        nsize[2-1] = size(dmd[i]['facecon']);
        nsize[3-1] = size(eblks);
        nsize[4-1] = size(fblks);
        nsize[5-1] = size(dmd[i]['nbsd']);
        nsize[6-1] = size(dmd[i]['elemsend']);
        nsize[7-1] = size(dmd[i]['elemrecv']);
        nsize[8-1] = size(dmd[i]['elemsendpts']);
        nsize[9-1] = size(dmd[i]['elemrecvpts']);
        nsize[10-1] = size(dmd[i]['elempart']);
        nsize[11-1] = size(dmd[i]['elempartpts']);
        nsize[12-1] = size(cgelcon);
        nsize[13-1] = size(rowent2elem);
        nsize[14-1] = size(cgent2dgent);
        nsize[15-1] = size(colent2elem);
        nsize[16-1] = size(rowe2f1);
        nsize[17-1] = size(cole2f1);
        nsize[18-1] = size(ent2ind1);
        nsize[19-1] = size(rowe2f2);
        nsize[20-1] = size(cole2f2);
        nsize[21-1] = size(ent2ind2);

        if (mpiprocs>1):
            print("Writing mesh into file " + str(i+1));
            fileID = open(filename + "/mesh" + str(i+1) + ".bin","w");
        else:
            print("Writing mesh into file ");
            fileID = open(filename + "/mesh" + ".bin","w");

        array(size(nsize), dtype=float64).tofile(fileID);
        nsize.astype('float64').tofile(fileID);
        ndims.astype('float64').tofile(fileID);
        dmd[i]['facecon'].transpose(1,0,2).flatten('F').astype(float64).tofile(fileID);
        eblks.flatten('F').astype(float64).tofile(fileID);
        fblks.flatten('F').astype(float64).tofile(fileID);
        if len(dmd[i]['nbsd'])>0:
            dmd[i]['nbsd'].flatten('F').astype(float64).tofile(fileID);
        if len(dmd[i]['elemsend'])>0:
            dmd[i]['elemsend'].flatten('F').astype(float64).tofile(fileID);
        if len(dmd[i]['elemrecv'])>0:
            dmd[i]['elemrecv'].flatten('F').astype(float64).tofile(fileID);
        if len(dmd[i]['elemsendpts'])>0:
            dmd[i]['elemsendpts'].flatten('F').astype(float64).tofile(fileID);
        if len(dmd[i]['elemrecvpts'])>0:
            dmd[i]['elemrecvpts'].flatten('F').astype(float64).tofile(fileID);
        if len(dmd[i]['elempart'])>0:
            dmd[i]['elempart'].flatten('F').astype(float64).tofile(fileID);
        if len(dmd[i]['elempartpts'])>0:
            dmd[i]['elempartpts'].flatten('F').astype(float64).tofile(fileID);
        cgelcon.flatten('F').astype(float64).tofile(fileID);
        rowent2elem.flatten('F').astype(float64).tofile(fileID);
        cgent2dgent.flatten('F').astype(float64).tofile(fileID);
        colent2elem.flatten('F').astype(float64).tofile(fileID);
        rowe2f1.flatten('F').astype(float64).tofile(fileID);
        cole2f1.flatten('F').astype(float64).tofile(fileID);
        ent2ind1.flatten('F').astype(float64).tofile(fileID);
        rowe2f2.flatten('F').astype(float64).tofile(fileID);
        cole2f2.flatten('F').astype(float64).tofile(fileID);
        ent2ind2.flatten('F').astype(float64).tofile(fileID);
        fileID.close();

    app = writeapp(app,fileapp);

    mesh['telem'] = master['telem'];
    mesh['tface'] = master['telem'];
    mesh['xpe'] = master['xpe'];
    mesh['xpf'] = master['xpf'];
    for i in range(0,mpiprocs):
        dmd[i]['nbsd'] = [];
        #dmd[i]['elem2cpu'] = [];
        dmd[i]['elemrecv'] = [];
        dmd[i]['elemsend'] = [];
        dmd[i]['elemrecvpts'] = [];
        dmd[i]['elemsendpts'] = [];
        dmd[i]['facepartpts'] = [];
        dmd[i]['facepartbnd'] = [];
        dmd[i]['facecon'] = [];

    return app, mesh, master, dmd
