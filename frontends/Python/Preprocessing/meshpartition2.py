import numpy as np
from elementpartition2 import elementpartition2
from facepartition2 import facepartition2

def meshpartition2(dmd,t,f,t2t,bcm,dim,elemtype,porder,nproc,metis,Cxxpreprocessing):

    print("run elementpartition...\n");
    dmd = elementpartition2(dmd,t,t2t,nproc,metis);

    if Cxxpreprocessing == 0:    
        print("run facepartition...\n");
        dmd = facepartition2(dmd,t,f,bcm,dim,elemtype,porder,nproc);
    else:
        for i in range(nproc):
            elempart = dmd[i]['elempart'].flatten();
            fi = f[:, elempart]  # Ensure correct indexing in Python
            dmd[i]['bf'] = np.zeros_like(fi)  # Initialize boundary flags    
            for j in range(len(bcm)):
                ind = np.where(fi == (j+1))  # Find indices where fi equals j
                dmd[i]['bf'][ind] = bcm[j]  # Assign bcm[j] to the corresponding indices

    return dmd;
