import Gencode, Preprocessing, os
from numpy import *

def gmshcall(pde, filename, nd, elemtype):

    opts = "-format msh3";

    # find gmsh executable
    gmsh = Preprocessing.findexec(pde['gmsh'], pde['version']);

    print("Gmsh mesh generator...\n");
    mystr = gmsh + " " + filename + ".geo -" + str(nd) + " " + opts;
    os.system(mystr);

    file1 = open(filename + ".msh", 'r')
    tm = file1.readlines();
    tm = [x.strip() for x in tm];
    file1.close();

    nlines = len(tm);
    j = 0;
    for i in range(0,nlines):
        if tm[i] == "$Nodes":
            j = i;
            break;
    np = int(tm[j+1]);
    p = zeros((nd,np));
    i = j;
    for ii in range(0,np):
        j = i+2+ii;
        pii = float64(tm[j].split());
        p[:,ii] = pii[1:(1+nd)];

    i = j+1;
    for ii in range(i,nlines):
        if tm[ii] == "$Elements":
            j = ii;
            break;

    ne = int(tm[j+1]);
    if nd==1: # line
        nve = 2;
        wcase = 1;
    elif nd==2 and elemtype==0: # tri
        nve = 3;
        wcase = 2;
    elif nd==2 and elemtype==1: # quad
        nve = 4;
        wcase = 3;
    elif nd==3 and elemtype==0: # tet
        nve = 4;
        wcase = 4;
    elif nd==3 and elemtype==1: # hex
        nve = 8;
        wcase = 5;

    t = zeros((nve,ne)).astype(int64);
    m = 0;
    i = j;
    for ii in range(0,ne):
        j = i+2+ii;
        tii = int64(tm[j].split());
        if (tii[1] == wcase):
            t[:,m] = tii[4:];
            m = m + 1;
    t = t[:,0:m];
    t = t - 1;

    return p, t
