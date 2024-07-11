import numpy as np
from ainb import ainb
from sortrows import sortrows
from partition import partition
from neighboringelements import neighboringelements

def elementpartitionhdg(dmd, t, t2t, nproc, metis):
    nve, ne = t.shape

    if nproc == 1:
        i = 0
        dmd[i]['elempart'] = np.arange(0,ne).astype(int);
        dmd[i]['elempartpts'] = np.array([ne]).astype(int);
        dmd[i]['nbsd'] = [];
        #dmd[i].elem2cpu = reshape([0],1,1);
        dmd[i]['elemrecv'] = [];
        dmd[i]['elemsend'] = [];
        dmd[i]['elemrecvpts'] = [];
        dmd[i]['elemsendpts'] = [];
        # #dmd[i]['elempart'] = np.reshape(np.arange(ne), (1, ne))
        # #dmd[i]['elempartpts'] = np.reshape([ne], (1, 1))
        # dmd[i]['nbsd'] = np.reshape([0], (1, 1))
        # # dmd[i]['elem2cpu'] = np.reshape([0], (1, 1))
        # dmd[i]['elemrecv'] = np.reshape([0], (1, 1))
        # dmd[i]['elemsend'] = np.reshape([0], (1, 1))
        # dmd[i]['elemrecvpts'] = np.reshape([0], (1, 1))
        # dmd[i]['elemsendpts'] = np.reshape([0], (1, 1))
        return dmd

    elem2cpu = partition(t+1,ne,nproc,metis)[0]
    elem2cpu = np.array(elem2cpu).flatten().astype(int);

    for i in range(nproc):
        intelem = np.where(elem2cpu == i)[0]  # elements in subdomain i
        intelem = np.array(intelem).astype(int).flatten();
        elem = neighboringelements(t2t, intelem)  # all elements connected to elements in subdomain i
        extelem = np.sort(np.setdiff1d(elem, intelem))  # exterior elements

        elem = neighboringelements(t2t, extelem)  # all elements connected to exterior elements
        bndelem = np.intersect1d(elem, intelem)  # boundary elements in subdomain i

        tmp = np.concatenate((np.setdiff1d(intelem, bndelem), bndelem, extelem))
        dmd[i]['elempart'] = tmp.flatten('F'); # partitions of elements        
        dmd[i]['elempartpts'] = np.array([len(intelem)-len(bndelem), len(bndelem), len(extelem)]).flatten();
        # dmd[i]['elem2cpu'] = elem2cpu[dmd[i]['elempart']]
        nelem = dmd[i]['elempartpts']

        recvelem = extelem
        tmp = np.array(np.arange(np.sum(nelem[0:2])+1,np.sum(nelem)+1)).flatten().astype(int);
        elemrecv = np.zeros((len(recvelem),3)).astype(int);
        elemrecv[:,0] = elem2cpu[recvelem];
        elemrecv[:,1] = tmp-1;
        elemrecv[:,2] = recvelem;
        tmp,ind = sortrows(elemrecv)[0:2];
        dmd[i]['elemrecv'] = elemrecv[ind.flatten('F'),:];
        dmd[i]['nbsd'] = (np.sort(np.unique(dmd[i]['elemrecv'][:,0]))).astype(int); # neighboring subdomains

    # store elements sent to neighboring subdomains to assemble the linear system
    for k in range(0,nproc):
        dmd[k]['elemsend'] = np.zeros(0);

    for i in range(0,nproc):
        for j in range(0,np.size(dmd[i]['nbsd'])):
            # cpu k sends information to cpu i
            k = dmd[i]['nbsd'][j];
            ii = np.nonzero(dmd[i]['elemrecv'][:,0] == k)[0];
            tm = dmd[i]['elemrecv'][ii,:];
            tm[:,0] = i;
            tm[:,1] = ainb(tm[:,2].flatten(), dmd[k]['elempart'].flatten());
            if np.size(dmd[k]['elemsend'])==0:
                dmd[k]['elemsend'] = tm;
            else:
                dmd[k]['elemsend'] = np.vstack((dmd[k]['elemsend'], tm));

    for i in range(0,nproc):
        dmd[i]['elemsendpts'] = np.zeros(np.size(dmd[i]['nbsd'])).astype(int);
        dmd[i]['elemrecvpts'] = np.zeros(np.size(dmd[i]['nbsd'])).astype(int);
        for j in range(0,np.size(dmd[i]['nbsd'])):
            dmd[i]['elemsendpts'][j] = np.size(np.nonzero(dmd[i]['elemsend'][:,0] == dmd[i]['nbsd'][j])[0]);
            dmd[i]['elemrecvpts'][j] = np.size(np.nonzero(dmd[i]['elemrecv'][:,0] == dmd[i]['nbsd'][j])[0]);

        dmd[i]['elemsend'] = dmd[i]['elemsend'][:,1];
        dmd[i]['elemrecv'] = dmd[i]['elemrecv'][:,1];

    return dmd  # elem2cpu; #elempart, elempartall

# Example usage:
# dmd = some_initial_value
# t, t2t, nproc, metis = some_other_values
# dmd = elementpartitionhdg(dmd, t, t2t, nproc, metis)
