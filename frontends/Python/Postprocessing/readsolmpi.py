import os
import numpy as np

def readsolmpi(base: str, nprocs: int, nsteps: int = 1, stepoffsets: int = 0):
    """
    Python equivalent of MATLAB readsolmpi.m

    Parameters
    ----------
    base : str
        Base filename (e.g., 'solution').
    nprocs : int
        Number of processor files to read.
    nsteps : int, optional
        Number of time steps to read (default: 1).
    stepoffsets : int, optional
        Number of time steps to skip before reading (default: 0).

    Returns
    -------
    k : ndarray
        Array of shape (nprocs, 4): [n1, n2, n3, timesteps] for each processor file.
    sol : ndarray, optional
        Combined solution array of shape (n1, n2, ne, nsteps),
        where ne = sum(k[:, 2]).
    """

    # --- Single-process shortcut
    if nprocs == 1:
        from pathlib import Path
        fname = f"{base}_np0.bin"
        if not Path(fname).exists():
            raise FileNotFoundError(f"Cannot open file: {fname}")

        with open(fname, "rb") as f:
            hdr = np.fromfile(f, dtype=np.float64, count=3)
            if hdr.size != 3:
                raise ValueError(f"Header read failed for {fname}")

            n1, n2, n3 = map(int, hdr)
            N = n1 * n2 * n3
            L = os.path.getsize(fname) / 8
            timesteps = int((L - 3) / N)

            sol = np.zeros((n1, n2, n3, nsteps), dtype=np.float64)
            if stepoffsets > 0:
                f.seek(stepoffsets * N * 8, os.SEEK_CUR)

            for i in range(nsteps):
                tm = np.fromfile(f, dtype=np.float64, count=N)
                sol[:, :, :, i] = tm.reshape((n1, n2, n3), order='F')

        k = np.array([[n1, n2, n3, timesteps]], dtype=float)
        return k, sol

    # --- Multi-process case
    k = np.zeros((nprocs, 4), dtype=float)

    # Pass 1: read headers and compute metadata
    for i in range(nprocs):
        fname = f"{base}_np{i}.bin"
        if not os.path.exists(fname):
            raise FileNotFoundError(f"Cannot open file: {fname}")

        with open(fname, "rb") as f:
            hdr = np.fromfile(f, dtype=np.float64, count=3)
            if hdr.size != 3:
                raise ValueError(f"Header read failed for {fname}")

            n1, n2, n3 = map(int, hdr)
            N = n1 * n2 * n3
            L = os.path.getsize(fname) / 8
            timesteps = int((L - 3) / N)
            k[i, :] = [n1, n2, n3, timesteps]

    # Consistency checks
    if not np.allclose(k[:, 0], k[0, 0]):
        raise ValueError("First dimension in binary files do not match.")
    if not np.allclose(k[:, 1], k[0, 1]):
        raise ValueError("Second dimension in binary files do not match.")
    if not np.allclose(k[:, 3], k[0, 3]):
        raise ValueError("Number of timesteps in binary files do not match.")

    # Pass 2: read and assemble full solution
    n1, n2 = int(k[0, 0]), int(k[0, 1])
    ne = int(np.sum(k[:, 2]))
    sol = np.zeros((n1, n2, ne, nsteps), dtype=np.float64)

    mn = 0
    for n in range(nprocs):
        fname = f"{base}_np{n}.bin"
        with open(fname, "rb") as f:
            hdr = np.fromfile(f, dtype=np.float64, count=3)
            if hdr.size != 3:
                raise ValueError(f"Header read failed for {fname}")

            n1, n2, n3 = map(int, hdr)
            N = n1 * n2 * n3

            if stepoffsets > 0:
                f.seek(stepoffsets * N * 8, os.SEEK_CUR)

            for i in range(nsteps):
                tm = np.fromfile(f, dtype=np.float64, count=N)
                if tm.size != N:
                    raise ValueError(f"Unexpected end of file while reading {fname}")
                sol[:, :, mn:mn+n3, i] = tm.reshape((n1, n2, n3), order='F')
            mn += n3

    return k, sol