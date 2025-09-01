import numpy as np

def read_rank(fname):
    """
    Read one rank's binary solution file.

    File format:
      - header: 3 float64 values [n1, n2, n3]
      - payload: nsteps blocks of size n1*n2*n3 (float64)
    """
    with open(fname, "rb") as f:
        # Read header (3 float64)
        hdr = np.fromfile(f, dtype=np.float64, count=3)
        if hdr.size != 3:
            raise ValueError(f"Header read failed for {fname}")
        n1, n2, n3 = hdr.astype(int)
        N = n1 * n2 * n3

        # Read the rest as doubles
        payload = np.fromfile(f, dtype=np.float64)

    if payload.size % N != 0:
        raise ValueError(
            f"Payload size {payload.size} is not a multiple of n1*n2*n3={N} in {fname}"
        )
    nsteps = payload.size // N

    # Reshape into [n1, n2, n3, nsteps]
    data4d = payload.reshape((n1, n2, n3, nsteps), order="F")  # MATLAB uses column-major

    return n1, n2, n3, nsteps, data4d