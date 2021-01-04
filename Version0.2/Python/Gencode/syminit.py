from sympy import *
from numpy import *

def syminit(app):

    nd = app['nd'];
    nc = app['nc'];
    ncx = app['ncx'];
    ncu = app['ncu'];
    ncq = app['ncq'];
    ncp = app['ncp'];
    nce = app['nce'];
    nco = app['nco'];
    ncw = app['ncw'];
    ntau = len(app['tau']);
    nuinf = len(app['uinf']);
    nparam = len(app['physicsparam']);

    time = array(symbols("time"));
    xdg = array(symbols("xdg1:" + str(ncx+1)));
    udg = array(symbols("udg1:" + str(nc+1)));
    udg1 = array(symbols("udgp1:" + str(nc+1)));
    udg2 = array(symbols("udgm1:" + str(nc+1)));
    wdg = array(symbols("wdg1:" + str(ncw+1)));
    wdg1 = array(symbols("wdgp1:" + str(ncw+1)));
    wdg2 = array(symbols("wdgm1:" + str(ncw+1)));
    odg = array(symbols("odg1:" + str(nco+1)));
    odg1 = array(symbols("odgp1:" + str(nco+1)));
    odg2 = array(symbols("odgm1:" + str(nco+1)));
    uhg = array(symbols("uhg1:" + str(ncu+1)));
    nlg = array(symbols("nlg1:" + str(nd+1)));
    tau = array(symbols("tau1:" + str(ntau+1)));
    uinf = array(symbols("uinf1:" + str(nuinf+1)));
    param = array(symbols("param1:" + str(nparam+1)));

    return xdg, udg, udg1, udg2, wdg, wdg1, wdg2, odg, odg1, odg2, uhg, nlg, tau, uinf, param, time
