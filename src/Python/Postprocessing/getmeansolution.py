from numpy import *

def getmeansolution(filename,dmd,npe):

    nproc = len(dmd);
    if nproc==1:
        fn = filename + "_np0.bin";
        UDG = fromfile(open(fn, "r"), dtype=float64);
        ne = len(dmd[0]['elempart'].flatten('F'));
        nc = int(round((size(UDG)-1)/(npe*ne)));
        UDG = reshape(UDG[0:npe*nc*ne]/UDG[-1],[npe,nc,ne],'F');
    else:
        nei = zeros(nproc);
        for i in range (0,nproc):
            nei[i] = sum(dmd[i]['elempartpts'][0:2]);
        ne = sum(nei);

        fn = filename + "_np0.bin";
        tmp = fromfile(open(fn, "r"), dtype=float64);
        nc = int(round((size(tmp)-1)/(npe*nei[0])));

        UDG = zeros((npe,nc,ne));
        for i in range(0,nproc):
            elempart = dmd[i]['elempart'][0:nei[i]];
            fn = filename + "_np" + str(i-1) + ".bin"
            tmp = fromfile(open(fn, "r"), dtype=float64);
            UDG[:,:,elempart]  = reshape(tmp[1:npe*nc*nei[i]]/tmp[-1],[npe,nc,nei[i]],'F');

    return UDG
