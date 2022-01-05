from numpy import *
from partition import partition
from neighboringelements import neighboringelements
from ainb import ainb
from writebin import writebin
from sortrows import sortrows

def elementpartition2(dmd,t,t2t,nproc,metis):

    nve,ne = t.shape;

    writebin('t.bin',t);
    writebin('t2t.bin',t2t);

    if nproc==1:
        i = 0;
        dmd[i]['elempart'] = arange(0,ne).astype(int);
        dmd[i]['elempartpts'] = array([ne]).astype(int);
        dmd[i]['nbsd'] = [];
        #dmd[i].elem2cpu = reshape([0],1,1);
        dmd[i]['elemrecv'] = [];
        dmd[i]['elemsend'] = [];
        dmd[i]['elemrecvpts'] = [];
        dmd[i]['elemsendpts'] = [];
        return dmd;

    elem2cpu = partition(t+1,ne,nproc,metis)[0]
    elem2cpu = array(elem2cpu).flatten().astype(int64);

    for i in range(0,nproc):
        intelem = nonzero(elem2cpu == i)[0]; # elements in subdomain i
        intelem = array(intelem).astype(int).flatten();
        elem = neighboringelements(t2t, intelem); # all elements connected to elements in subdomain i
        extelem = sort(setdiff1d(elem,intelem)); # exterior elements

        intextelem = concatenate([intelem, extelem]);
        elem = neighboringelements(t2t, extelem); # all elements connected to exterior elements
        bndelem = intersect1d(elem,intelem);  # boundary elements in subdmain i
        outelem = sort(setdiff1d(elem,intextelem)); #  elements outside subdmain i

        tmp = concatenate([setdiff1d(intelem,bndelem), bndelem, extelem, outelem]);
        dmd[i]['elempart'] = tmp.flatten('F'); # partitions of elements
        dmd[i]['elempartpts'] = array([len(intelem)-len(bndelem), len(bndelem), len(extelem), len(outelem)]).flatten();
        #dmd[i].elem2cpu = elem2cpu[dmd[i].elempart];
        nelem = dmd[i]['elempartpts'];

        recvelem = concatenate([extelem, outelem]); # elements received from neighbors
        tmp = array(arange(sum(nelem[0:2])+1,sum(nelem)+1)).flatten().astype(int);
        elemrecv = zeros((len(recvelem),3)).astype(int);
        elemrecv[:,0] = elem2cpu[recvelem];
        elemrecv[:,1] = tmp-1;
        elemrecv[:,2] = recvelem;
        tmp,ind = sortrows(elemrecv)[0:2];
        dmd[i]['elemrecv'] = elemrecv[ind.flatten('F'),:];
        dmd[i]['nbsd'] = (sort(unique(dmd[i]['elemrecv'][:,0]))).astype(int); # neighboring subdomains

        #writebin("elemrecv" + str(i+1) +  ".bin",dmd[i]['elemrecv']);

    # store elements sent to neighboring subdomains to assemble the linear system
    for k in range(0,nproc):
        dmd[k]['elemsend'] = zeros(0);

    for i in range(0,nproc):
        for j in range(0,size(dmd[i]['nbsd'])):
            # cpu k sends information to cpu i
            k = dmd[i]['nbsd'][j];
            ii = nonzero(dmd[i]['elemrecv'][:,0] == k)[0];
            tm = dmd[i]['elemrecv'][ii,:];
            tm[:,0] = i;
            tm[:,1] = ainb(tm[:,2].flatten(), dmd[k]['elempart'].flatten());
            if size(dmd[k]['elemsend'])==0:
                dmd[k]['elemsend'] = tm;
            else:
                dmd[k]['elemsend'] = vstack((dmd[k]['elemsend'], tm));

    for i in range(0,nproc):
        dmd[i]['elemsendpts'] = zeros(size(dmd[i]['nbsd'])).astype(int);
        dmd[i]['elemrecvpts'] = zeros(size(dmd[i]['nbsd'])).astype(int);
        for j in range(0,size(dmd[i]['nbsd'])):
            dmd[i]['elemsendpts'][j] = size(nonzero(dmd[i]['elemsend'][:,0] == dmd[i]['nbsd'][j])[0]);
            dmd[i]['elemrecvpts'][j] = size(nonzero(dmd[i]['elemrecv'][:,0] == dmd[i]['nbsd'][j])[0]);

        dmd[i]['elemsend'] = dmd[i]['elemsend'][:,1];
        dmd[i]['elemrecv'] = dmd[i]['elemrecv'][:,1];

    return dmd
