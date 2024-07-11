import numpy as np
import sys
from faceconnectivity2 import faceconnectivity2
from mkf2e import mkf2e

def facepartitionhdg(dmd, t, f, bcm, dim, elemtype, porder, nproc):
    
    for i in range(nproc):
        print(f"face partition {i}")
        elempart = dmd[i]['elempart'].flatten();

        fi = f[:, elempart]  # Ensure correct indexing in Python
        dmd[i]['bf'] = np.zeros_like(fi)  # Initialize boundary flags    
        for j in range(len(bcm)):
            ind = np.where(fi == (j+1))  # Find indices where fi equals j
            dmd[i]['bf'][ind] = bcm[j]  # Assign bcm[j] to the corresponding indices
        
        f2t, _ = mkf2e(t[:, elempart], elemtype, dim)

        # Only on [intelem bndelem]
        nelem = dmd[i]['elempartpts']
        if len(nelem) >= 2:
            ne3 = np.sum(nelem[:2])
            ind = np.where(f2t[0, :] <= ne3)[0]
            f2t = f2t[:, ind]

        # Reorder so that boundary faces are last
        ina = np.where(f2t[2, :] > 0)[0]  # Interior faces
        inb = np.where(f2t[2, :] == 0)[0]  # Boundary faces
        inc = (f2t[0, inb] - 1) * fi.shape[0] + f2t[1, inb] - 1  # Global indices of boundary faces
        fb = fi.flatten('F')[inc.flatten()]  # List of boundary indices

        fa = np.sort(np.unique(fb))  # Boundary indices
        bcn = np.sort(np.unique(bcm[fa-1]))  # A list of boundary conditions
        nbc = len(bcn)

        dmd[i]['facepartpts'] = np.array([len(ina)]).reshape(1, 1)
        dmd[i]['facepartbnd'] = np.array([0]).reshape(1, 1)
        ind = np.zeros(len(fb), dtype=int)
        m = 0
        for j in range(nbc):  # For each boundary condition bcn[j]
            bj = np.where(bcm == bcn[j])[0] + 1 # Find all boundaries that have condition bcn[j]
            n = 0
            for k in range(len(bj)):  # For each boundary that has condition bcn[j]
                ii = np.where(fb == bj[k])[0]  # Indices of the boundary bj[k]
                l = len(ii)
                n += l
                ind[m:m + l] = ii
                m += l
            dmd[i]['facepartpts'] = np.append(dmd[i]['facepartpts'], [n])
            dmd[i]['facepartbnd'] = np.append(dmd[i]['facepartbnd'], [bcn[j]])
                
        # [interior faces, boundary faces]
        f2t = f2t[:, np.concatenate((ina, inb[ind[:]]))]
        dmd[i]['f2t'] = f2t    

        # Fix f2t to ensure DOF consistency across interior faces
        # The global ID of the FIRST face must be smaller than that of the SECOND face        
        N = dmd[i]['f2t'].shape[1]
        for face in range(N):  # Loop over each face
            e1 = dmd[i]['f2t'][0, face]
            e2 = dmd[i]['f2t'][2, face]
            if e2 > 0:  # If the face is an interior face                
                g1 = elempart[e1-1]
                g2 = elempart[e2-1]
                if g2 < g1:
                    # Ensure g2 > g1
                    tm = dmd[i]['f2t'][2:4, face].copy()
                    dmd[i]['f2t'][2:4, face] = dmd[i]['f2t'][0:2, face]
                    dmd[i]['f2t'][0:2, face] = tm
                e1 = dmd[i]['f2t'][0, face]
                e2 = dmd[i]['f2t'][2, face]
                g1 = elempart[e1-1]
                g2 = elempart[e2-1]
                if g2 <= g1:
                    raise ValueError("something wrong")

        dmd[i]['facecon'], dmd[i]['elemcon'] = faceconnectivity2(t[:, elempart], dmd[i]['f2t'], dim, elemtype, porder)

    return dmd
