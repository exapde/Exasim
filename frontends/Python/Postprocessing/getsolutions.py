import numpy as np

def getsolutions(basename, dmd=None):
    """
    Assemble solution blocks from per-rank binary files.

    File format per rank:
      header (3 float64): [n1, n2, n3]
      followed by data for all timesteps, each block of size n1*n2*n3 (float64)

    Parameters
    ----------
    basename : str
        Base path/prefix of files. Files are "<basename>_np<rank>.bin".
    dmd : list
        Partition metadata per rank. Each entry must provide:
          - elempartpts: array-like, use sum(elempartpts[:2]) for element count
          - elempart:    array-like global element indices (1-based or 0-based; see note)
        If None, assumes single-rank file "<basename>_np0.bin".

    Returns
    -------
    sol : np.ndarray
        Shape [n1, n2, net, nsteps] for multi-rank, or [n1, n2, n3, nsteps] for single-rank.
        Fortran order (column-major) to mirror MATLAB reshape semantics.
    """
    def _f(x, name):
        return x[name] if isinstance(x, dict) else getattr(x, name)

    nproc = 1 if dmd is None else len(dmd)

    if nproc == 1:
        n1, n2, n3, nsteps, block = read_rank(rankfile(basename, 0))
        # block is [n1, n2, n3, nsteps]
        return block

    # ---- multi-rank assembly ----
    # how many elements each rank contributes
    nei = np.zeros(nproc, dtype=int)
    for r in range(nproc):
        epp = _f(dmd[r], 'elempartpts')
        nei[r] = int(np.sum(epp[:2]))

    net = int(np.sum(nei))

    # read rank 0 for sizing
    n1, n2, _, nsteps, block0 = read_rank(rankfile(basename, 0))
    sol = np.zeros((n1, n2, net, nsteps), dtype=np.float64, order='F')

    # place rank 0
    elempart0 = np.asarray(_f(dmd[0], 'elempart')[:nei[0]])
    # If your elempart is 1-based (MATLAB style), convert to 0-based:
    if elempart0.min() == 1:
        elempart0 = elempart0 - 1
    sol[:, :, elempart0, :] = block0

    # place ranks 1..nproc-1
    for r in range(1, nproc):
        fname = rankfile(basename, r)
        n1r, n2r, n3r, nsteps_r, blockr = read_rank(fname)
        # (Optional) consistency checks—comment out if layouts can differ:
        # if n1r != n1 or n2r != n2 or nsteps_r != nsteps:
        #     raise ValueError(f"Rank {r} has incompatible sizes.")

        elempart = np.asarray(_f(dmd[r], 'elempart')[:nei[r]])
        if elempart.min() == 1:
            elempart = elempart - 1
        sol[:, :, elempart, :] = blockr

    return sol


def rankfile(base, r):
    """Build '<base>_np<r>.bin'."""
    return f"{base}_np{r}.bin"


def read_rank(fname):
    """
    Read one rank's binary solution file.

    Format:
      - header: 3 float64 values [n1, n2, n3]
      - payload: nsteps blocks of size n1*n2*n3 (float64)
    """    
    with open(fname, "rb") as f:
        hdr = np.fromfile(f, dtype=np.float64, count=3)
        if hdr.size != 3:
            raise ValueError(f"Header read failed for {fname}")
        n1, n2, n3 = hdr.astype(int)
        N = n1 * n2 * n3

        payload = np.fromfile(f, dtype=np.float64)

    if payload.size % N != 0:
        raise ValueError(
            f"Payload size {payload.size} is not a multiple of n1*n2*n3={N} in {fname}"
        )
    nsteps = payload.size // N

    # Reshape into [n1, n2, n3, nsteps] with Fortran order to match MATLAB
    payload = payload.reshape((n1, n2, n3, nsteps), order="F")

    return n1, n2, n3, nsteps, payload

# import os
# import numpy as np

# def getsolutions(pde, dmd):
#     """
#     Python version of MATLAB getsolutions.m
#     - pde: dict with at least pde['buildpath']
#     - dmd: list (length = nproc) of dicts; each dict must have:
#         dmd[r]['elempart']      : 1D integer array of global element indices for this rank
#         dmd[r]['elempartpts']   : array-like with at least 2 entries; sum of first two = local ne
#     Returns:
#         UDG  : shape (npe, nc,  net, nsteps) if nproc>1 else (npe, nc,  ne, nsteps)
#         WDG  : shape (npe, ncw, net, nsteps) if nproc>1 else (npe, ncw, ne, nsteps)
#         UHAT : shape (ncu, npf, nft, nsteps)
#     Notes:
#       * Assumes the binary files are written in the same layout as your MATLAB code
#         (header of 10 doubles, then blocks in column-major order).
#       * Uses Fortran ('F') order reshapes to match MATLAB column-major semantics.
#       * If your dmd[r]['elempart'] is 1-based, subtract 1 before using it for NumPy indexing.
#     """
#     def binpath(rank):
#         return os.path.join(pde['buildpath'], "dataout", f"outsol_np{rank}.bin")

#     def read_all_doubles(path):
#         return np.fromfile(path, dtype=np.float64)

#     def read_header10(path):
#         return np.fromfile(path, dtype=np.float64, count=10)

#     nproc = len(dmd)

#     # --- Single-process case -------------------------------------------------
#     if nproc == 1:
#         tmp = read_all_doubles(binpath(0))
#         if tmp.size < 10:
#             raise RuntimeError("File too short: missing 10-double header in rank 0 file.")

#         # Header (MATLAB 1-based: tmp(1..10)) -> Python 0-based: tmp[0..9]
#         npe = int(tmp[0])
#         nc  = int(tmp[1])
#         ne  = int(tmp[2])
#         ncw = int(tmp[4])
#         ncu = int(tmp[6])
#         npf = int(tmp[7])
#         nf  = int(tmp[8])
#         N   = int(tmp[9])

#         # Steps
#         total_payload = tmp.size - 10
#         if total_payload % N != 0:
#             raise RuntimeError("File payload not divisible by N; inconsistent file layout.")
#         nsteps = total_payload // N

#         # Allocate
#         UDG  = np.zeros((npe, nc,  ne,  nsteps), dtype=np.float64)
#         WDG  = np.zeros((npe, ncw, ne,  nsteps), dtype=np.float64)
#         UHAT = np.zeros((ncu, npf, nf, nsteps), dtype=np.float64)

#         n1 = npe * nc  * ne
#         n2 = npe * ncw * ne
#         n3 = ncu * npf * nf

#         m = 10  # 0-based start after header (MATLAB had m=11 in 1-based)
#         for i in range(nsteps):
#             UDG[:, :, :, i]  = np.reshape(tmp[m : m+n1],      (npe, nc,  ne),  order='F'); m += n1
#             WDG[:, :, :, i]  = np.reshape(tmp[m : m+n2],      (npe, ncw, ne),  order='F'); m += n2
#             UHAT[:, :, :, i] = np.reshape(tmp[m : m+n3],      (ncu, npf, nf),  order='F'); m += n3

#         return UDG, WDG, UHAT

#     # --- Multi-process case --------------------------------------------------
#     # Per-rank element counts (from dmd[r]['elempartpts'][0:2])
#     nei = np.zeros(nproc, dtype=np.int64)
#     for r in range(nproc):
#         # Sum of first two entries
#         nei[r] = int(np.sum(np.asarray(dmd[r]['elempartpts'])[:2]))

#     net = int(np.sum(nei))

#     # Read rank 0 fully (we also need header values)
#     tmp0 = read_all_doubles(binpath(0))
#     if tmp0.size < 10:
#         raise RuntimeError("File too short: missing 10-double header in rank 0 file.")

#     npe = int(tmp0[0])
#     nc  = int(tmp0[1])
#     ne0 = int(tmp0[2])
#     ncw = int(tmp0[4])
#     ncu = int(tmp0[6])
#     npf = int(tmp0[7])
#     nf0 = int(tmp0[8])
#     N0  = int(tmp0[9])

#     total_payload0 = tmp0.size - 10
#     if total_payload0 % N0 != 0:
#         raise RuntimeError("Rank 0 file payload not divisible by N0; inconsistent layout.")
#     nsteps = total_payload0 // N0

#     # Collect nf for all ranks (including rank 0) to size UHAT properly
#     nf_list = [nf0]
#     for r in range(1, nproc):
#         hdr = read_header10(binpath(r))
#         if hdr.size != 10:
#             raise RuntimeError(f"File too short: missing header in rank {r} file.")
#         nf_list.append(int(hdr[8]))
#     nft = int(np.sum(nf_list))

#     # Allocate global arrays
#     UDG  = np.zeros((npe, nc,  net, nsteps), dtype=np.float64)
#     WDG  = np.zeros((npe, ncw, net, nsteps), dtype=np.float64)
#     UHAT = np.zeros((ncu, npf, nft, nsteps), dtype=np.float64)

#     # Helper to import one rank’s payload into global arrays
#     def fill_rank(rank, tmp, elem_offset_slice, face_offset_start):
#         ne_r = int(tmp[2])     # header slot 3 in MATLAB
#         nf_r = int(tmp[8])     # header slot 9 in MATLAB
#         N    = int(tmp[9])     # header slot 10 in MATLAB

#         n1 = npe * nc  * ne_r
#         n2 = npe * ncw * ne_r
#         n3 = ncu * npf * nf_r

#         m = 10
#         # Sanity check the per-step chunk size
#         if N != (n1 + n2 + n3):
#             # Some simulations may legitimately vary N across ranks; we do not error out,
#             # but we do verify we can step through the file.
#             pass

#         # Steps (length consistency)
#         payload = tmp.size - 10
#         if payload % N != 0:
#             raise RuntimeError(f"Rank {rank} payload not divisible by its N; inconsistent layout.")
#         steps_here = payload // N
#         if steps_here != nsteps:
#             raise RuntimeError(f"Rank {rank} has {steps_here} steps but rank 0 has {nsteps}.")

#         for i in range(nsteps):
#             # UDG block
#             UDG[:, :, elem_offset_slice, i] = np.reshape(
#                 tmp[m : m+n1], (npe, nc, ne_r), order='F'
#             )
#             m += n1

#             # WDG block
#             WDG[:, :, elem_offset_slice, i] = np.reshape(
#                 tmp[m : m+n2], (npe, ncw, ne_r), order='F'
#             )
#             m += n2

#             # UHAT block (faces append contiguously across ranks)
#             UHAT[:, :, face_offset_start:face_offset_start+nf_r, i] = np.reshape(
#                 tmp[m : m+n3], (ncu, npf, nf_r), order='F'
#             )
#             m += n3

#         return nf_r  # how many faces we filled

#     # Rank 0 fill
#     # NOTE: If your dmd[0]['elempart'] is 1-based, subtract 1: elempart0 = np.asarray(dmd[0]['elempart'][:nei[0]]) - 1
#     elempart0 = np.asarray(dmd[0]['elempart'][:nei[0]], dtype=np.int64)
#     k_face = 0
#     k_face += fill_rank(0, tmp0, elempart0, k_face)

#     # Remaining ranks
#     for r in range(1, nproc):
#         tmp_r = read_all_doubles(binpath(r))
#         elempart_r = np.asarray(dmd[r]['elempart'][:nei[r]], dtype=np.int64)
#         k_face += fill_rank(r, tmp_r, elempart_r, k_face)

#     return UDG, WDG, UHAT