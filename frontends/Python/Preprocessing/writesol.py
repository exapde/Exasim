from numpy import *

def writesol(filename,ifile,xdg, udg, vdg, wdg):

    
    ncx = xdg.shape[1];

    if isinstance(udg, list) and udg == []:
        ncu = 0;
    else:
        ncu = udg.shape[1];
    
    if isinstance(vdg, list) and vdg == []:
        nco = 0;
    else:
        nco = vdg.shape[1];
    
    if isinstance(wdg, list) and wdg == []:
        ncw = 0;
    else:
        ncw = wdg.shape[1];

    ndims = zeros((12,1));
    ndims[1-1] = ncx;
    ndims[2-1] = ncu;
    ndims[3-1] = nco;
    ndims[4-1] = ncw;

    nsize = zeros((20,1));
    nsize[1-1] = size(ndims);
    nsize[2-1] = size(xdg);
    nsize[3-1] = size(udg);
    nsize[4-1] = size(vdg);
    nsize[5-1] = size(wdg);

    if (ifile>0):
        print("Writing solution into file " + str(ifile));
        fileID = open(filename + str(ifile) + ".bin","wb");
    else:
        print("Writing solution into file");
        fileID = open(filename + ".bin","wb");

    array(size(nsize), dtype=float64).tofile(fileID);
    nsize.astype('float64').tofile(fileID);
    ndims.astype('float64').tofile(fileID);
    if nsize[1]>0:
        xdg = array(xdg).flatten(order = 'F');
        xdg.astype('float64').tofile(fileID);
    if nsize[2]>0:
        udg = array(udg).flatten(order = 'F');
        udg.astype('float64').tofile(fileID);
    if nsize[3]>0:
        vdg = array(vdg).flatten(order = 'F');
        vdg.astype('float64').tofile(fileID);
    if nsize[4]>0:
        wdg = array(wdg).flatten(order = 'F');
        wdg.astype('float64').tofile(fileID);

    fileID.close();
