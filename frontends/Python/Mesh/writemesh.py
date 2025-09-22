from numpy import *

def writemesh(p, t, filename, mode):

    nd,np = p.shape;
    nve,ne = t.shape;
    tmp = int64(array([nd, np, nve, ne]));

    if mode==0: # binary
        fileID = open(filename, "wb");
        tmp.astype('float64').tofile(fileID);
        p = p.flatten(order = 'F');
        p.astype('float64').tofile(fileID);
        t = t.flatten(order = 'F');
        t.astype('float64').tofile(fileID);
        fileID.close();
    else: # text
        fileID = open(filename, "w");
        savetxt(fileID, tmp.reshape((1,4)), fmt='%d');
        savetxt(fileID, p.T, fmt='%.16f');
        savetxt(fileID, t.T, fmt='%d');
        fileID.close();
    return 0;
