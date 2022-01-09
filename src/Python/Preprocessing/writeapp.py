from numpy import *

def writeapp(app,filename):

    app['flag'] = array(app['flag']);
    app['problem'] = array(app['problem']);
    app['factor'] = array(app['factor']);
    app['solversparam'] = array(app['solversparam']);

    appname = 0;
    tmp = array([app['tdep'], app['wave'], app['linearproblem'], app['debugmode'], app['matvecorder'], app['GMRESortho'], app['preconditioner'], app['precMatrixType'], app['NLMatrixType'], app['runmode'], app['tdfunc'], app['source'], app['modelnumber']]);
    app['flag'] =  concatenate([tmp,app['flag']]);
    tmp = array([app['hybrid'], appname, app['temporalscheme'], app['torder'], app['nstage'], app['convStabMethod'], app['diffStabMethod'], app['rotatingFrame'], app['viscosityModel'], app['SGSmodel'], app['ALE'], app['AV'], app['linearsolver'], app['NLiter'], app['linearsolveriter'], app['GMRESrestart'], app['RBdim'], app['saveSolFreq'], app['saveSolOpt'], app['timestepOffset'], app['stgNmode'], app['saveSolBouFreq'], app['ibs'], app['dae_steps'], app['saveResNorm'], app['AVsmoothingIter'], app['frozenAVflag']]);
    app['problem'] = concatenate([tmp, app['problem']]);
    tmp = array([app['time'], app['dae_alpha'], app['dae_beta'], app['dae_gamma'], app['dae_epsilon']])    
    app['factor'] = concatenate([tmp, app['factor']]);
    tmp = array([app['NLtol'], app['linearsolvertol'], app['matvectol'], app['NLparam']]);
    app['solversparam'] = concatenate([tmp, app['solversparam']]);

    app['flag'] = array(app['flag']).flatten('F');
    app['problem'] = array(app['problem']);
    app['factor'] = array(app['factor']);
    app['solversparam'] = array(app['solversparam']);

    ndims = zeros((40,1));
    ndims[1-1] = app['mpiprocs'];  # number of processors
    ndims[2-1] = app['nd'];
    ndims[3-1] = 0;
    ndims[4-1] = 0;
    ndims[5-1] = 0;
    ndims[6-1] = app['nc'];
    ndims[7-1] = app['ncu'];
    ndims[8-1] = app['ncq'];
    ndims[9-1] = app['ncp'];
    ndims[10-1] = app['nco'];
    ndims[11-1] = app['nch'];
    ndims[12-1] = app['ncx'];
    ndims[13-1] = app['nce'];
    ndims[14-1] = app['ncw'];

    #if app['nco'] != app['vindx'].shape[0]:  #size(app.vindx,1):
    #    error("app.nco mus be equal to size(app.vindx,1)");

    nsize = zeros((20,1));
    nsize[1-1] = size(ndims);
    nsize[2-1] = size(app['flag']);  # size of flag
    nsize[3-1] = size(app['problem']); # size of physics
    nsize[4-1] = size(app['uinf']); # boundary data
    nsize[5-1] = size(app['dt']); # number of time steps
    nsize[6-1] = size(app['factor']); # size of factor
    nsize[7-1] = size(app['physicsparam']); # number of physical parameters
    nsize[8-1] = size(app['solversparam']); # number of solver parameters
    nsize[9-1] = size(app['tau']); # number of solver parameters
    nsize[10-1] = size(app['stgdata']);
    nsize[11-1] = size(app['stgparam']);
    nsize[12-1] = size(app['stgib']);
    nsize[13-1] = size(app['vindx']);
    nsize[14-1] = size(app['dae_dt']); # number of dual time steps

    print("Writing app into file...");
    fileID = open(filename, 'wb');
    array(size(nsize), dtype=float64).tofile(fileID)
    nsize.astype('float64').tofile(fileID)
    if nsize[1-1] > 0:
        ndims.astype('float64').tofile(fileID);
    if nsize[2-1] > 0:
        app['flag'] = array(app['flag']).flatten(order = 'F');
        app['flag'].astype('float64').tofile(fileID);
    if nsize[3-1] > 0:
        app['problem'] = array(app['problem']).flatten(order = 'F');
        app['problem'].astype('float64').tofile(fileID);
    if nsize[4-1] > 0:
        app['uinf'] = array(app['uinf']).flatten(order = 'F');
        app['uinf'].astype('float64').tofile(fileID);
    if nsize[5-1] > 0:
        app['dt'] = array(app['dt']).flatten(order = 'F');
        app['dt'].astype('float64').tofile(fileID);
    if nsize[6-1] > 0:
        app['factor'] = array(app['factor']).flatten(order = 'F');
        app['factor'].astype('float64').tofile(fileID);
    if nsize[7-1] > 0:
        app['physicsparam'] = array(app['physicsparam']).flatten(order = 'F');
        app['physicsparam'].astype('float64').tofile(fileID);
    if nsize[8-1] > 0:
        app['solversparam'] = array(app['solversparam']).flatten(order = 'F');
        app['solversparam'].astype('float64').tofile(fileID);
    if nsize[9-1] > 0:
        app['tau'] = array(app['tau']).flatten(order = 'F');
        app['tau'].astype('float64').tofile(fileID);
    if nsize[10-1] > 0:
        app['stgdata'] = array(app['stgdata']).flatten(order = 'F');
        app['stgdata'].astype('float64').tofile(fileID);
    if nsize[11-1] > 0:
        app['stgparam'] = array(app['stgparam']).flatten(order = 'F');
        app['stgparam'].astype('float64').tofile(fileID);
    if nsize[12-1] > 0:
        app['stgib'] = array(app['stgib']).flatten(order = 'F');
        app['stgib'].astype('float64').tofile(fileID);
    if nsize[13-1] > 0:
        app['vindx'] = array(app['vindx']).flatten(order = 'F')-1;
        app['vindx'].astype('float64').tofile(fileID);
    if nsize[14-1] > 0:
        app['dae_dt'] = array(app['dae_dt']).flatten(order = 'F');
        app['dae_dt'].astype('float64').tofile(fileID);

    fileID.close();

    return app;
