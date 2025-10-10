import os
import numpy as np

def readsol(base, nsteps=1, stepoffsets=0):
    """
    Python equivalent of MATLAB readsol function.

    Parameters
    ----------
    base : str
        Base filename (without processor suffix, e.g., "solution").
    nsteps : int, optional
        Number of time steps to read. Default is 1.
    stepoffsets : int, optional
        Number of timesteps to skip from the beginning. Default is 0.

    Returns
    -------
    n1, n2, n3 : int
        Dimensions of the solution.
    timesteps : int
        Total number of timesteps in the file.
    sol : ndarray, optional
        Solution array of shape (n1, n2, n3, nsteps).
    """

    fname = f"{base}_np0.bin"

    if not os.path.exists(fname):
        raise FileNotFoundError(f"Cannot open file: {fname}")

    with open(fname, "rb") as f:
        # Read header (3 doubles)
        hdr = np.fromfile(f, dtype=np.float64, count=3)
        if hdr.size != 3:
            raise ValueError(f"Header read failed for {fname}")

        n1, n2, n3 = map(int, hdr)
        N = n1 * n2 * n3

        # Determine total number of timesteps
        L = os.path.getsize(fname) / 8  # total doubles in file
        timesteps = int((L - 3) / N)

        sol = None
        if nsteps is not None and nsteps > 0:
            if stepoffsets > 0:
                f.seek(stepoffsets * N * 8, os.SEEK_CUR)

            sol = np.zeros((n1, n2, n3, nsteps), dtype=np.float64)
            for i in range(nsteps):
                tm = np.fromfile(f, dtype=np.float64, count=N)
                if tm.size != N:
                    raise ValueError(f"Unexpected end of file while reading step {i+1} in {fname}")
                sol[:, :, :, i] = tm.reshape((n1, n2, n3), order='F')

    return n1, n2, n3, timesteps, sol