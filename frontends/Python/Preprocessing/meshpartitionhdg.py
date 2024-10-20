from elementpartitionhdg import elementpartitionhdg
from facepartitionhdg import facepartitionhdg

def meshpartitionhdg(dmd, t, f, t2t, bcm, dim, elemtype, porder, nproc, metis):
    print("run elementpartition...")
    dmd = elementpartitionhdg(dmd, t, t2t, nproc, metis)

    print("run facepartition...")
    dmd = facepartitionhdg(dmd, t, f, bcm, dim, elemtype, porder, nproc)

    return dmd

# Example usage:
# dmd = some_initial_value
# t, f, t2t, bcm, dim, elemtype, porder, nproc, metis = some_other_values
# dmd = meshpartitionhdg(dmd, t, f, t2t, bcm, dim, elemtype, porder, nproc, metis)
