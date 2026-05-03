# Known test divergences

This file catalogs tests that fail on lanka but pass on Mac, with
diagnosis as either **platform numerics** (real solver, just
different last-bit precision than the recorded baseline) or **real
divergence** (a bug to investigate).

Now that the test harness uses element-L2 baselines, most of the
"platform numerics" failures from earlier runs no longer appear —
element-L2 is FP-order-tolerant, so per-platform drift in the
1e-7 to 1e-3 range passes cleanly. What's left here are bugs that
need fixing.

Updated: 2026-05-03.

## Open

### `orion:mpi`, `orion:mpi_gpu` — QoI ~62% off serial

**Fail signal**: `outqoi.txt` Domain_QoI1 = `2.030988e-01` (lanka MPI
np=2) vs baseline `1.255923e-01` (Mac serial). Identical for both
codegen + facade stems → not facade-specific.

**Smoking gun**: `orion:cpu` and `orion:gpu` pass byte-identical
against the same baseline. The divergence appears only under MPI.
Ratio `0.2030988 / 0.1255923 ≈ 1.617` is consistent with each rank
summing over `owned + ghost` elements (with ~62% ghost overlap for
np=2) instead of just owned. Plausible cause: `QoIvolume` evaluation
in the runtime iterates `[0, ne)` where `ne == ne_full` (owned +
ghosts) rather than `ne1` (owned only).

**Investigation pointer**: `backend/Solution/postsolution.hpp`
`SaveQoI` calls `qoiElement<M>(...)` which uses `disc.common.ne`.
Check whether that includes ghost rows.

**Status**: real bug, should be fixed before relying on QoI-based
MPI baselines.

## Closed (recently)

### `periodic:mpi`, `periodic:mpi_gpu` — was 32.6% off, now ~1e-7

**Was failing** with `relative element-L2 = 3.258e-01` over 64
elements (both codegen + facade stems lockstep), serial passing.

**Cause**: parallel `setperiodicfaces` in
`backend/Preprocessing/parmetisexasim.hpp` updated the element-to-
element neighbor array `t2t` with global element IDs of periodic
neighbors but never merged the global node IDs across the periodic
boundary. The runtime's `buildConn` matches faces by sorted vertex
IDs (via `mkf2e_hash`), so periodic-paired faces with non-merged
node IDs looked like physical boundaries.

**Fix**: added `mergePeriodicNodeIDs(mesh, comm)` in the same file.
Collects boundary-1 and boundary-2 nodes per rank, q-coord-matches
them via `xiny()` (same matching `setperiodicfaces` does on faces),
and rewrites `mesh.nodeGlobalID` so each pair shares a canonical
(smaller) global ID. Transitive resolution handles 3-way matches at
xy-periodic corners. Called from `initializeDMD` before
`setperiodicfaces`.

**Result**: `relative element-L2` now ~1.7e-7 at np=2, ~1.4e-7 at
np=4. Five orders of magnitude better. Both codegen + facade stems
pass.

## Platform-numerics divergences (now passing via element-L2)

The following examples used to fail on lanka with bin-md5 gating
but now pass under element-L2. Magnitudes for reference (from a
recent lanka run):

| Example | Variant | element-L2 vs Mac baseline |
|---|---|---|
| isoq3d | cpu, gpu | 2.6e-7 |
| isoq | cpu, gpu | 2.3e-7 (cpu), 8.6e-8 (gpu) |
| sharpb2 | cpu, gpu | (passes; magnitude similar) |
| nsmach8 | cpu, gpu | 3.8e-12 |
| cone | gpu | 9.4e-7 |
| naca0012steady | cpu, gpu | 2.6e-11 |
| naca0012unsteady | cpu, gpu | 5.4e-11 |

These are pure platform-numerics drift — both `codegen:X:Y` and
`facade:X:Y` produce bit-for-bit identical output on lanka, just
last-bit different from the Mac baseline (Apple Accelerate vs Intel
MKL + AVX-512, accumulated through GMRES + Newton). Element-L2's
1e-3 threshold cleanly admits this drift while still catching real
bugs at 1e-2 and above.

## Validation notes

The strongest validation signal is **codegen vs facade lockstep**:
bit-identical output between the two stems proves the facade is
correct relative to the legacy `<exasim/run.hpp>` path. Cross-
platform agreement is a separate question; element-L2 with a 1e-3
threshold catches real divergences while admitting platform drift.

Apps with analytical solutions (Shu/Yee isentropic vortex, Method
of Manufactured Solutions, etc.) would give us cross-platform AND
solver-correctness validation simultaneously. None of the current
`apps/library_example/*_codegen` tests are set up that way; adding
manufactured-solution variants is queued as future work.
