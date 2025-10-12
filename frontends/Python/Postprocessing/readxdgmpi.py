import numpy as np
from pathlib import Path

def readbin(filename):
    """
    Reads a binary file of doubles, equivalent to MATLAB readbin().
    """
    with open(filename, "rb") as f:
        data = np.fromfile(f, dtype=np.float64)
    return data


def readxdg(filesol):
    """
    Python equivalent of MATLAB readxdg.m

    Parameters
    ----------
    filesol : str
        Path to the binary solution file.

    Returns
    -------
    xdg : ndarray
        Array of shape (npe, ncx, ne)
    """

    tm = readbin(filesol)

    # Step 1: read header sizes
    sz = int(tm[0])

    k1 = 1
    k2 = k1 + sz
    nsize = tm[k1:k2].astype(int)

    k1 = k2
    k2 = k1 + nsize[0]
    ndims = tm[k1:k2].astype(int)

    k1 = k2
    k2 = k1 + nsize[1]
    xdg = tm[k1:k2]

    # Step 2: reshape according to dimensions
    ne = ndims[0]
    npe = ndims[3]
    ncx = int(len(xdg) / (npe * ne))

    xdg = np.reshape(xdg, (npe, ncx, ne), order='F')
    return xdg

def readxdgmpi(base, nprocs, ne):
    """
    Python equivalent of MATLAB readxdgmpi.m

    Parameters
    ----------
    base : str
        Base filename (e.g., 'xdg' for files like 'xdg1.bin', 'xdg2.bin', ...).
    nprocs : int
        Number of processor files to read.

    Returns
    -------
    xdg : ndarray
        Concatenated array of shape (npe, ncx, total_ne)
    """

    if nprocs == 1:
        filesol = f"{base}.bin"
        if not Path(filesol).exists():
            raise FileNotFoundError(f"Cannot open file: {filesol}")
        xdg = readxdg(filesol)
        return xdg

    xdg = None
    for i in range(1, nprocs + 1):
        filesol = f"{base}{i}.bin"
        if not Path(filesol).exists():
            raise FileNotFoundError(f"Cannot open file: {filesol}")
        xdgi = readxdg(filesol)
        idx = np.array(ne[i - 1]) - 1  
        xdgi = xdgi[:, :, idx]

        if xdg is None:
            xdg = xdgi;
        else:
            xdg = np.concatenate((xdg, xdgi), axis=2)

    return xdg

