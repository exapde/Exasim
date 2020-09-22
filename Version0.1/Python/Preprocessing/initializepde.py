import numpy

def initializepde(version):

    pde = {'codename' : [], 'version' : []};

    pde['cpucompiler']  = "g++";
    pde['mpicompiler'] = "mpicxx";
    pde['gpucompiler'] = "nvcc";
    pde['mpirun'] = "mpirun";
    pde['gmsh'] = "gmsh";
    pde['metis'] = "mpmetis";
    pde['paraview'] = "paraview";

    pde['codename'] = "Exasim";
    pde['version'] = version;
    pde['appname'] = "pde";
    pde['model'] = "ModelD";
    pde['modelfile'] = "";
    pde['platform'] = "cpu";
    #pde['cpuflags'] = "-O2 -Wall -ldl -lm -lblas -llapack";
    pde['cpuflags'] = "-O2 -ldl -lm -lblas -llapack";
    pde['gpuflags'] = "-lcudart -lcublas";
    pde['preprocessmode'] = 1;
    pde['mpiprocs'] = 1;
    pde['nd'] = 1;
    pde['nc'] = 1;
    pde['ncu'] = 1;
    pde['ncq'] = 0;
    pde['ncp'] = 0;
    pde['nco'] = 0;
    pde['nch'] = 1;
    pde['ncx'] = 1;
    pde['ncw'] = 0;
    pde['nce'] = 0;
    pde['neb'] = 512*8;
    pde['nfb'] = 512*32;
    pde['elemtype'] = 1;
    pde['nodetype'] = 1;
    pde['hybrid'] = 0;
    pde['tdep'] = 0;
    pde['wave'] = 0;
    pde['linearproblem'] = 0;
    pde['debugmode'] = 0;
    pde['stgNmode'] = 0;
    pde['porder'] = 1;
    pde['pgauss'] = 2;
    pde['temporalscheme'] = 0;
    pde['torder'] = 1;
    pde['nstage'] = 1;
    pde['convStabMethod'] = 0;
    pde['diffStabMethod'] = 0;
    pde['rotatingFrame'] = 0;
    pde['viscosityModel'] = 0;
    pde['SGSmodel'] = 0;
    pde['ALE'] = 0;
    pde['AV'] = 0;
    pde['nonlinearsolver'] = 0;
    pde['linearsolver'] = 0;
    pde['NLiter'] = 20;
    pde['linearsolveriter'] = 200;
    pde['GMRESrestart'] = 25;
    pde['GMRESortho'] = 0;
    pde['preconditioner'] = 0;
    pde['precMatrixType'] = 0;
    pde['NLMatrixType'] = 0;
    pde['runmode'] = 0;
    pde['tdfunc'] = 1;
    pde['source'] = 1;
    pde['matvecorder'] = 1;
    pde['RBdim'] = 5;
    pde['saveSolFreq'] = 1;
    pde['saveSolOpt'] = 1;
    pde['timestepOffset'] = 0;
    pde['saveSolBouFreq'] = 0;
    pde['ibs'] = 0;
    pde['compudgavg'] = 0;
    pde['extFhat'] = 0;
    pde['extUhat'] = 0;

    pde['time'] = 0.0;
    pde['NLparam'] = 0.0;
    pde['NLtol'] = 1e-6;
    pde['linearsolvertol'] = 1e-3;
    pde['matvectol'] = 1e-3;

    pde['flag'] = [0, 0];
    pde['problem'] = [0, 0];
    pde['boundaryconditions'] = [0, 0];
    pde['stgib'] = [0, 0];

    pde['tau'] = 1.0; # stabilization parameters
    pde['uinf'] = numpy.array([0.0, 0.0]); # freestream values
    pde['externalparam'] = numpy.array([0.0, 0.0]); # external parameters
    pde['dt'] = numpy.array([0.0]);  # time steps
    pde['factor'] = numpy.array([0.0, 0.0]);  # factors
    pde['physicsparam'] = numpy.array([0.0, 0.0]); # physical parameters
    pde['solversparam'] = numpy.array([0.0, 0.0]); # solvers parameters
    pde['stgdata'] = numpy.array([0.0, 0.0]); # synthetic turbulence
    pde['stgparam'] = numpy.array([0.0, 0.0]); # synthetic turbulence

    pde['viselem'] = [];
    pde['visfilename'] = "dataout/output";  # filename for paraview output files

    return pde
