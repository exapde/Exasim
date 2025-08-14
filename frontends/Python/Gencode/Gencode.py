
import os, sys, numpy
from sympy import *
from importlib import import_module
import importlib.util
from checkcompilers import checkcompilers
from setcompilers import setcompilers
from compilecode import compilecode
from cmakecompile import cmakecompile
from compilepdemodel import compilepdemodel
from runcode import runcode
from gencodeall import gencodeall
from syminit import syminit
from gencodeelemface import gencodeelemface
from gencodeelem import gencodeelem
from nocodeelem import nocodeelem
from gencodeelem2 import gencodeelem2
from nocodeelem2 import nocodeelem2
from gencodeelem3 import gencodeelem3
from nocodeelem3 import nocodeelem3
from gencodeelem4 import gencodeelem4
from nocodeelem4 import nocodeelem4
from gencodeface import gencodeface
from nocodeface import nocodeface
from gencodeface2 import gencodeface2
from gencodeface3 import gencodeface3
from nocodeface2 import nocodeface2
from hdggencodeelem import hdggencodeelem
from hdgnocodeelem import hdgnocodeelem
from hdggencodeelem2 import hdggencodeelem2
from hdgnocodeelem2 import hdgnocodeelem2
from hdggencodeface import hdggencodeface
from hdgnocodeface import hdgnocodeface
from hdggencodeface2 import hdggencodeface2
from hdgnocodeface2 import hdgnocodeface2

def gencode(app):

    if app['codegenerator'] == "text2code":
        runstr = os.path.join(app['exasimpath'], "text2code", "text2code", "text2code")
        inputfile = app['modelfile'] + ".txt"
        full_cmd = [runstr, inputfile]
        subprocess.run(full_cmd, check=True)
        return

    # foldername = os.path.join(app['exasimpath'], "build", "model")
    foldername = os.path.join(app['backendpath'], "Model")

    # Ensure the target directory exists
    os.makedirs(foldername, exist_ok=True)

    # Read the content of the file
    with open(os.path.join(app['backendpath'], "Discretization", "KokkosDrivers.cpp"), "r") as file:
        text = file.read()

    # Write the content to a new file
    with open(os.path.join(foldername, "KokkosDrivers.cpp"), "w") as fid:
        fid.write(text)

    xdg, udg, udg1, udg2, wdg, wdg1, wdg2, odg, odg1, odg2, uhg, nlg, tau, uinf, param, time = syminit(app)[0:16];

    pde = import_module(app['modelfile']);

    if app['modelnumber']==0:
        strn = "";
    else:
        strn = str(app['modelnumber']);

    nc = app['nc'];
    ncu = app['ncu'];
    u = udg[0:ncu];
    u1 = udg1[0:ncu];
    u2 = udg2[0:ncu];
    if nc>ncu:
        q = udg[ncu:];
        q1 = udg1[ncu:];
        q2 = udg2[ncu:];
    else:
        q = [];
        q1 = [];
        q2 = [];

    if hasattr(pde, 'flux'):
        #f = pde.flux(xdg, udg, odg, wdg, uinf, param, time);
        f = pde.flux(u, q, wdg, odg, xdg, time, param, uinf);
        f = f.flatten('F');
        gencodeelem("Flux" + strn, f, xdg, udg, odg, wdg, uinf, param, time, foldername);
        if app['hybrid'] == 1:
            hdggencodeelem("Flux" + str(strn), f, xdg, udg, odg, wdg, uinf, param, time, foldername)
        else:
            hdgnocodeelem("Flux" + str(strn), foldername)    
    else:
        sys.exit('pde.flux is not defined')
    if hasattr(pde, 'source'):
        #f = pde.source(xdg, udg, odg, wdg, uinf, param, time);
        f = pde.source(u, q, wdg, odg, xdg, time, param, uinf);
        gencodeelem("Source" + strn, f, xdg, udg, odg, wdg, uinf, param, time, foldername);
        if app['hybrid'] == 1:
            hdggencodeelem("Source" + str(strn), f, xdg, udg, odg, wdg, uinf, param, time, foldername)
        else:
            hdgnocodeelem("Source" + str(strn), foldername)    
    else:
        nocodeelem("Source" + strn, foldername);
        hdgnocodeelem("Source" + str(strn), foldername)    
    if hasattr(pde, 'eos'):
        f = pde.eos(u, q, wdg, odg, xdg, time, param, uinf);
        gencodeelem2("EoS" + strn, f, xdg, udg, odg, wdg, uinf, param, time, foldername);      

        nf = len(f);
        nu = len(u);
        nw = len(wdg);
        
        dfdu = array(symbols("dfdu1:" + str(nu*nf+1)));
        k = 0;
        for n in range(0, nu): 
          for m in range(0, nf):  
            dfdu[k] = diff(f[m],u[n]);      
            k = k + 1;
          end
        end        
        gencodeelem2("EoSdu" + strn, dfdu, xdg, udg, odg, wdg, uinf, param, time, foldername);      

        dfdw = array(symbols("dfdw1:" + str(nw*nf+1)));
        k = 0;
        for n in range(0, nw): 
          for m in range(0, nf):  
            dfdw[k] = diff(f[m],w[n]);      
            k = k + 1;
          end
        end        
        gencodeelem2("EoSdw" + strn, dfdw, xdg, udg, odg, wdg, uinf, param, time, foldername);      
        if app['hybrid'] == 1:
            hdggencodeelem("EoS" + str(strn), f, xdg, udg, odg, wdg, uinf, param, time, foldername)
        else:
            hdgnocodeelem("EoS" + str(strn), foldername)    
    else:
        nocodeelem2("EoS" + strn, foldername);
        nocodeelem2("EoSdu" + strn, foldername);
        nocodeelem2("EoSdw" + strn, foldername);
        hdgnocodeelem("EoS" + strn, foldername);
    if hasattr(pde, 'sourcew'):
        f = pde.sourcew(u, q, wdg, odg, xdg, time, param, uinf);
        gencodeelem2("Sourcew" + strn, f, xdg, udg, odg, wdg, uinf, param, time, foldername);
        if app['hybrid'] == 1:
            hdggencodeelem("Sourcew" + str(strn), f, xdg, udg, odg, wdg, uinf, param, time, foldername)
            hdggencodeelem2("Sourcewonly" + str(strn), f, xdg, udg, odg, wdg, uinf, param, time, foldername)
        else:
            hdgnocodeelem("Sourcew" + str(strn), foldername)
            hdgnocodeelem2("Sourcewonly" + str(strn), foldername)    
    else:
        nocodeelem2("Sourcew" + strn, foldername);
        hdgnocodeelem("Sourcew" + strn, foldername)
        hdgnocodeelem2("Sourcewonly" + strn, foldername)
    if hasattr(pde, 'mass'):
        #f = pde.mass(xdg, udg, odg, wdg, uinf, param, time);
        f = pde.mass(u, q, wdg, odg, xdg, time, param, uinf);
        gencodeelem("Tdfunc" + strn, f, xdg, udg, odg, wdg, uinf, param, time, foldername);
    else:
        if app['model']=="ModelW" or app['model']=="modelW" or app['tdep']==1:
            error("pde.inituq is not defined");
        else:
            nocodeelem("Tdfunc" + strn, foldername);
    if hasattr(pde, 'avfield'):
        #f = pde.avfield(xdg, udg, odg, wdg, uinf, param, time);
        f = pde.avfield(u, q, wdg, odg, xdg, time, param, uinf);
        gencodeelem2("Avfield" + strn, f, xdg, udg, odg, wdg, uinf, param, time, foldername);
    else:
        nocodeelem2("Avfield" + strn, foldername);
    if hasattr(pde, 'output'):
        #f = pde.output(xdg, udg, odg, wdg, uinf, param, time);
        f = pde.output(u, q, wdg, odg, xdg, time, param, uinf);
        gencodeelem2("Output" + strn, f, xdg, udg, odg, wdg, uinf, param, time, foldername);
    else:
        nocodeelem2("Output" + strn, foldername);
    if hasattr(pde, 'monitor'):
        f = pde.monitor(u, q, wdg, odg, xdg, time, param, uinf);
        gencodeelem2("Monitor" + strn, f, xdg, udg, odg, wdg, uinf, param, time, foldername);
    else:
        nocodeelem2("Monitor" + strn, foldername);
    if app['hybrid'] == 1:
        if hasattr(pde, 'fbouhdg'):
            f = pde.fbouhdg(u, q, wdg, odg, xdg, time, param, uinf, uhg, nlg, tau)
            f = numpy.reshape(f.flatten('F'),(app['ncu'],round(f.size/app['ncu'])),'F');
            hdggencodeface("Fbou" + str(strn), f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, foldername)
            hdggencodeface2("Fbouonly" + str(strn), f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, foldername)
        else:
            raise AttributeError("pde.fbouhdg is not defined")
        if hasattr(pde, 'fint'):
            f = pde.fint(u, q, wdg, odg, xdg, time, param, uinf, uhg, nlg, tau)
            ncu12 = pde['interfacefluxmap'].size
            if ncu<=0 :
              ncu12 = 1
            f = numpy.reshape(f.flatten('F'),(ncu12,round(f.size/ncu12)),'F');
            hdggencodeface("Fint" + str(strn), f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, foldername)
            hdggencodeface2("Fintonly" + str(strn), f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, foldername)
        else:
            hdgnocodeface("Fint" + str(strn), foldername)
            hdgnocodeface2("Fintonly" + str(strn), foldername)    
    else:
        hdgnocodeface("Fbou" + str(strn), foldername)
        hdgnocodeface2("Fbouonly" + str(strn), foldername)    
    if hasattr(pde, 'fbou'):
        #f = pde.fbou(xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time);
        f = pde.fbou(u, q, wdg, odg, xdg, time, param, uinf, uhg, nlg, tau);
        f = numpy.reshape(f.flatten('F'),(app['ncu'],round(f.size/app['ncu'])),'F');
        gencodeface("Fbou" + strn, f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, foldername);
    else:
        sys.exit("pde.fbou is not defined");
    if hasattr(pde, 'ubou'):
        #f = pde.ubou(xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time);
        f = pde.ubou(u, q, wdg, odg, xdg, time, param, uinf, uhg, nlg, tau);
        f = numpy.reshape(f.flatten('F'),(app['ncu'],round(f.size/app['ncu'])),'F');
        gencodeface("Ubou" + strn, f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, foldername);
    else:
        sys.exit("pde.ubou is not defined");
    if hasattr(pde, 'fhat'):
        #f = pde.fhat(xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
        f = pde.fhat(u1, q1, wdg1, odg1, xdg, time, param, uinf, uhg, nlg, tau, u2, q2, wdg2, odg2);
        gencodeface2("Fhat" + strn, f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, foldername);
    else:
        nocodeface2("Fhat" + strn, foldername);
    if hasattr(pde, 'uhat'):
        #f = pde.uhat(xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
        f = pde.uhat(u1, q1, wdg1, odg1, xdg, time, param, uinf, uhg, nlg, tau, u2, q2, wdg2, odg2);
        gencodeface2("Uhat" + strn, f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, foldername);
    else:
        nocodeface2("Uhat" + strn, foldername);
    if hasattr(pde, 'stab'):
        #f = pde.stab(xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
        f = pde.stab(u1, q1, wdg1, odg1, xdg, time, param, uinf, uhg, nlg, tau, u2, q2, wdg2, odg2);
        gencodeface3("Stab" + strn, f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, foldername);
    else:
        nocodeface2("Stab" + strn, foldername);
    if hasattr(pde, 'initu'):
        udg = pde.initu(xdg, param, uinf);
        gencodeelem3("Initu" + strn, udg, xdg, uinf, param, foldername);
        gencodeelem4("Initu" + strn, udg, xdg, uinf, param, foldername);
    else:
        error("pde.initu is not defined");
    if hasattr(pde, 'initw'):
        wdg = pde.initw(xdg, param, uinf);
        gencodeelem3("Initwdg" + strn, wdg, xdg, uinf, param, foldername);
        gencodeelem4("Initwdg" + strn, wdg, xdg, uinf, param, foldername);
    else:
        nocodeelem3("Initwdg" + strn, foldername);
        nocodeelem4("Initwdg" + strn, foldername);
    if hasattr(pde, 'initv'):
        odg = pde.initv(xdg, param, uinf);
        gencodeelem3("Initodg" + strn, odg, xdg, uinf, param, foldername);
        gencodeelem4("Initodg" + strn, odg, xdg, uinf, param, foldername);
    else:
        nocodeelem3("Initodg" + strn, foldername);
        nocodeelem4("Initodg" + strn, foldername);

    if hasattr(pde, 'initq'):
        q = pde.initq(xdg, param, uinf);
        q = q.flatten('F');
        gencodeelem3("Initq" + strn, q, xdg, uinf, param, foldername);
        gencodeelem4("Initq" + strn, q, xdg, uinf, param, foldername);

        u = pde.initu(xdg, param, uinf);

        udg = numpy.concatenate((u.flatten('F'), q.flatten('F')));
        gencodeelem3("Initudg" + strn, udg, xdg, uinf, param, foldername);
        gencodeelem4("Initudg" + strn, udg, xdg, uinf, param, foldername);
    else:
        if app['model']=="ModelW" or app['model']=="modelW":
            error("pde.initq is not defined");
        else:
            nocodeelem3("Initq" + strn, foldername);
            nocodeelem3("Initudg" + strn, foldername);
            nocodeelem4("Initq" + strn, foldername);
            nocodeelem4("Initudg" + strn, foldername);

    # if hasattr(pde, 'initq'):
    #     udg = pde.initq(xdg, param, uinf);
    #     gencodeelem3("Initq", udg, xdg, uinf, param);
    # else:
    #     nocodeelem3("Initq");
    # if hasattr(pde, 'inituq'):
    #     udg = pde.inituq(xdg, param, uinf);
    #     gencodeelem3("Initudg", udg, xdg, uinf, param);
    # else:
    #     if app['model']=="ModelW" or app['model']=="modelW":
    #         error("pde.inituq is not defined");
    #     else:
    #         nocodeelem3("Initudg");

    return 0;
