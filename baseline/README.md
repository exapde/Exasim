# Numerical regression baseline (Phase 0.2)

Captured outputs from the existing `cput2cEXASIM` and `cpumpit2cEXASIM`
executables before the library port begins. Every later phase must
reproduce these (text outputs bit-for-bit; binary md5s identical for a
deterministic run, or numerical diff within a tight tolerance).

## How it was captured

```
# from repo root, after building Kokkos / METIS / text2code / Exasim
cd apps/poisson/poisson2d
build/text2code pdeapp.txt           # generates Model/* + libpdemodelserial.dylib + datain/*.bin
build/cput2cEXASIM     pdeapp.txt    # serial run
mpirun -np 2 build/cpumpit2cEXASIM pdeapp.txt   # 2-rank MPI run
```

## Reference numerics

The Poisson 2D test problem solves `-Δu = 2π² sin(πx) sin(πy)` on the
unit square with `u=0` on the boundary; exact solution
`u = sin(πx) sin(πy)`. The QoI text output is the same for serial and
MPI runs:

```
Time            Domain_QoI1     Domain_QoI2     Boundary_QoI1
0.000000e+00    4.992401e-13    4.052847e-01    8.000000e+00
```

`Domain_QoI1` is the L² error against the analytic solution
(`(uq[0]-uexact)²` integrated). `4.99e-13` confirms the discretization
converged to machine precision on this mesh + porder=3 setup.

## What's stored

```
baseline/
  poisson2d_serial/
    outqoi.txt        # the QoI table above
    md5.txt           # md5 hashes of every dataout/*.bin
  poisson2d_mpi2/
    outqoi.txt        # identical QoI table
    md5.txt           # md5 hashes of every dataout/*.bin (np0, np1)
  verify.sh           # re-runs both, diffs against this snapshot
```

## How to verify after a refactor

```
bash baseline/verify.sh
```

Exit code 0 = numerics preserved. Exit code != 0 = regression; inspect
the diff output.
