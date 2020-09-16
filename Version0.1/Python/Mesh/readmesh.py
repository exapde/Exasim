from numpy import *

def readmesh(filename, mode):

    #fileID = fopen(filename,'r');
    if mode==0: # binary
        tmp = fromfile(open(filename, "r"), dtype=float64);
        nd = int_(tmp[0]);
        np = int_(tmp[1]);
        nve = int_(tmp[2]);
        ne = int_(tmp[3]);
        p = reshape(tmp[4:(4+nd*np)], (nd, np), order='F');
        t = reshape(tmp[(4+nd*np):(4+nd*np+nve*ne)], (nve, ne), order='F');
        t = int_(t);
    else: # text
        file1 = open(filename, "r");
        tmp = file1.readline();
        tmp = int64(tmp.strip().split());
        nd = tmp[0];
        np = tmp[1];
        nve = tmp[2];
        ne = tmp[3];
        file1.close();
        p = genfromtxt(filename, skip_header = 1, skip_footer = ne);
        t = int64(genfromtxt(filename, skip_header = 1+np));
        p = p.T;
        t = t.T;
    return p,t
