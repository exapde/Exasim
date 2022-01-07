from numpy import *
from partition import partition
from getelemface import getelemface
from mkv2t import mkv2t
from node2elem import node2elem
from ainb import ainb

def elementpartition(dmd,t,elemtype,nproc,metis):

    nve,ne = t.shape;

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
    elem2cpu = array(elem2cpu).flatten();

    if nve==2:
        dim=1;
    else:
        if elemtype==0: # tri/tet elements
            dim=nve-1;
        elif elemtype==1: # quad/hex elements
            dim=log2(nve);

    face = getelemface(dim,elemtype);
    nvf,nfe = shape(face);

    print("run mkv2t...");
    re,ce = mkv2t(t,ne);

    for i in range(0,nproc):
        intelem = nonzero(elem2cpu == i)[0]; # elements in subdomain i
        intelem = array(intelem).astype(int).flatten();
        elem = node2elem(t[:,intelem],re,ce);
        extelem = sort(setdiff1d(elem,intelem)); # exterior elements

        # fix extelem
        n1 = size(intelem);
        n2 = size(extelem);
        t1 = reshape(t[ix_(face.flatten(order='F'),intelem.flatten(order='F'))], (nvf, nfe*n1), order='F');
        t2 = reshape(t[ix_(face.flatten(order='F'),extelem.flatten(order='F'))], (nvf, nfe, n2), order='F');
        t1 = sort(t1,axis=0);
        t2 = sort(t2,axis=0);
        match = zeros(n2);
        for j in range(0,n2):
            for k in range(0,nfe):
                if min(sum(abs(reshape(t2[:,k,j],(nvf,1)) - t1),axis=0))==0:
                    match[j] = 1;
                    break;
        ind = nonzero(match==1)[0];
        extelem = extelem[ind];

        intextelem = concatenate([intelem, extelem]);
        elem = node2elem(t[:,extelem],re,ce); # all elements connected to exterior elements
        bndelem = intersect1d(elem,intelem);  # boundary elements in subdmain i
        outelem = sort(setdiff1d(elem,intextelem)); #  elements outside subdmain i

        # fix outelem
        n1 = size(intextelem);
        n2 = size(outelem);
        t1 = reshape(t[ix_(face.flatten(order='F'),intextelem)], (nvf, nfe*n1), order='F');
        t2 = reshape(t[ix_(face.flatten(order='F'),outelem)], (nvf, nfe, n2), order='F');
        t1 = sort(t1,axis=0);
        t2 = sort(t2,axis=0);
        match = zeros(n2);
        for j in range(0,n2):
            for k in range(0,nfe):
                if min(sum(abs(reshape(t2[:,k,j],(nvf,1)) - t1),axis=0))==0:
                    match[j] = 1;
                    break;
        ind = nonzero(match==1)[0];
        outelem = outelem[ind];

        tmp = concatenate([setdiff1d(intelem,bndelem), bndelem, extelem, outelem]);
        dmd[i]['elempart'] = tmp.flatten('F'); # partitions of elements
        dmd[i]['elempartpts'] = array([len(intelem)-len(bndelem), len(bndelem), len(extelem), len(outelem)]).flatten();
        #dmd[i].elem2cpu = elem2cpu[dmd[i].elempart];
        nelem = dmd[i]['elempartpts'];

        recvelem = concatenate([extelem, outelem]); # elements received from neighbors
        tmp = array(arange(sum(nelem[0:2])+1,sum(nelem)+1)).flatten();
        dmd[i]['elemrecv'] = vstack((elem2cpu[recvelem], tmp-1, recvelem)).transpose();
        ind = argsort(dmd[i]['elemrecv'][:,0]);
        dmd[i]['elemrecv'] = dmd[i]['elemrecv'][ind,:];
        dmd[i]['nbsd'] = (sort(unique(dmd[i]['elemrecv'][:,0]))).astype(int); # neighboring subdomains

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
        dmd[i]['elemsendpts'] = zeros(size(dmd[i]['nbsd']));
        dmd[i]['elemrecvpts'] = zeros(size(dmd[i]['nbsd']));
        for j in range(0,size(dmd[i]['nbsd'])):
            dmd[i]['elemsendpts'][j] = size(nonzero(dmd[i]['elemsend'][:,0] == dmd[i]['nbsd'])[0]);
            dmd[i]['elemrecvpts'][j] = size(nonzero(dmd[i]['elemrecv'][:,0] == dmd[i]['nbsd'])[0]);
        dmd[i]['elemsend'] = dmd[i]['elemsend'][:,1];
        dmd[i]['elemrecv'] = dmd[i]['elemrecv'][:,1];

    # for i in range(0,nproc):
    #     print(dmd[i]['nbsd'])
    #     print(dmd[i]['elemrecv'])
    #     print(dmd[i]['elemsend'])
    #     print(dmd[k]['elempart'])

    return dmd
