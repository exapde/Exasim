# Known test divergences

This file catalogs tests that fail on lanka but pass on Mac, with
diagnosis as either **platform numerics** (real solver, just
different last-bit precision than the recorded baseline) or **real
divergence** (a bug to investigate).

Updated: 2026-05-02 (HOT.7.16 lanka run).

## Real divergences (TODO)

### `periodic:mpi`, `periodic:mpi_gpu` — element-L2 ~32.6% off serial

**Fail signal**: `relative element-L2 = 3.258e-01 (>1e-03)` over 64
elements. Identical magnitude on both codegen + facade stems →
not facade-specific. Both serial variants (`periodic:cpu`,
`periodic:gpu`) pass with element-L2 ~ 3e-14.

**Smoking gun**: divergence appears only under MPI partitioning.
The `periodic` test exercises periodic boundary conditions
(periodicboundaries1/2 in pdeapp.txt). Suspected cause: the
periodic-face-pairing logic in `parmetisexasim.hpp` /
`connectivity.hpp` may not correctly pair faces that are split
across rank partitions, leading to local boundary conditions that
differ from serial.

**Status**: real bug, related family to orion:mpi (both are
partition-handling regressions, both lockstep). Investigate
together as HOT.7.17.

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

**Status**: real bug, not blocking HOT.7.x but should be fixed
before relying on QoI-based MPI baselines.

## Platform-numerics divergences (no action needed)

These all show **lockstep behavior** — `codegen:X:Y` and
`facade:X:Y` produce bit-for-bit identical output. Only the
recorded baseline (Mac, Apple Accelerate) differs from lanka's
output (Intel MKL + AVX-512), drift accumulated through GMRES
+ Newton iterations to a few parts in 1e-3.

| Example | Variant | rel-RMS vs Mac baseline | Severity |
|---|---|---|---|
| isoq3d | cpu, gpu | 4.796e-04 | small |
| isoq | cpu, gpu | 4.561e-04 | small |
| sharpb2 | cpu, gpu | 2.867e-04 | small |
| nsmach8 | cpu, gpu | 3.850e-03 | medium |
| cone | gpu | 8.947e-03 | medium |
| cone | cpu | (similar to gpu) | medium |

For these the validation strategy is `bin md5` (no QoI baseline),
which is sensitive to FP order. Two ways to make them
cross-platform-stable:

1. **Add `QoIvolume` / `QoIboundary` to each app's `pdemodel.txt`**
   and re-record `outqoi.txt` baselines. The validate harness
   already prefers `outqoi.txt` over bin md5 (and tolerates 1e-4
   relative diff on individual QoI components). Domain integrals
   of mass/momentum/energy are FP-stable to ~1e-8.
2. **Record per-platform baselines** (`baseline/<name>_serial_lanka/`,
   `baseline/<name>_serial_mac/`) and have validate select based
   on `uname` / hostname. More fragile; more state.

Recommended: option 1 for the apps above. Option 2 is a fallback
for apps where adding a QoI is non-trivial (e.g. moving boundaries).

## Validation notes

The strongest validation signal we have today is **codegen vs
facade lockstep**: bit-identical output between the two stems
proves the facade is correct relative to the legacy
`<exasim/run.hpp>` path. Cross-platform agreement is a separate
question that depends on having a partition- and FP-order-invariant
quantity to compare.

Apps with analytical solutions (Shu/Yee isentropic vortex, Method
of Manufactured Solutions, etc.) would give us cross-platform
*and* solver-correctness validation simultaneously. None of the
current `apps/library_example/*_codegen` tests are set up that
way.
