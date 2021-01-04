from numpy import *

def getsolution(filename,dmd,npe):

    nproc = len(dmd);
    if nproc==1:
        fn = filename + "_np0.bin";
        UDG = fromfile(open(fn, "r"), dtype=float64);
        ne = len(dmd[0]['elempart'].flatten('F'));
        nc = int(round(size(UDG)/(npe*ne)));
        UDG = reshape(UDG,[npe,nc,ne],'F');
    else:
        nei = zeros(nproc).astype(int);
        for i in range (0,nproc):
            nei[i] = sum(dmd[i]['elempartpts'][0:2]);
        ne = sum(nei);

        fn = filename + "_np0.bin";
        tmp = fromfile(open(fn, "r"), dtype=float64);
        nc = int(round(size(tmp)/(npe*nei[0])));

        UDG = zeros((npe,nc,ne));
        for i in range(0,nproc):
            elempart = dmd[i]['elempart'][0:nei[i]].flatten('F')-1;
            fn = filename + "_np" + str(i) + ".bin"
            tmp = fromfile(open(fn, "r"), dtype=float64);
            UDG[:,:,elempart]  = reshape(tmp,(npe,nc,nei[i]),'F');

    return UDG
