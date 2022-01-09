from elementpartition import elementpartition
from facepartition import facepartition

def meshpartition(dmd,p,t,f,tprd,elemtype,bcm,bndexpr,prdexpr,porder,nproc,metis):

    # display("run facenumbering...");
    # f, tprd = facenumbering(p,t,elemtype,bndexpr,prdexpr);

    print("run elementpartition...");
    if tprd.shape==t.shape:
        dmd = elementpartition(dmd,tprd,elemtype,nproc,metis);
    else:
        dmd = elementpartition(dmd,t,elemtype,nproc,metis);

    print("run facepartition...");
    dmd = facepartition(dmd,p,t,tprd,f,bcm,elemtype,prdexpr,porder,nproc);

    return dmd
