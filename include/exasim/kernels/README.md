# `<exasim/kernels/>` — templated model kernels (Phase 3)

Empty during Phases 1–2. Populated in Phase 3 with one header per kernel
family:

```
flux.hpp           kokkos_flux_kernel<M>,    hdg_flux_kernel<M>
source.hpp         kokkos_source_kernel<M>,  hdg_source_kernel<M>
sourcew.hpp
boundary.hpp       fbou, ubou, fhat, uhat, stab kernels
init.hpp           initu, initq, initwdg, initudg, initodg
tdfunc.hpp
eos.hpp            eos, eos_du, eos_dw, avfield
output.hpp         output, monitor
visualization.hpp  vis_scalars, vis_vectors, vis_tensors
qoi.hpp            qoi_volume, qoi_boundary
```

Each header contains 1–2 `template <class M> void *_kernel(...)` functions
that perform SoA gather → call `M::flux(...)` (and, for HDG, `M::flux_jac_uq`)
→ SoA scatter, inside a `Kokkos::parallel_for`. They are pure plumbing —
the user supplies the math (and Jacobians) on `M`. **No autodiff path.**

Reference template (will be generalized from `backend/Model/KokkosFlux1.cpp`):

```cpp
namespace exasim {
template <class M>
void kokkos_flux_kernel(dstype* f, const dstype* xdg, const dstype* udg,
                        const dstype* odg, const dstype* wdg,
                        const dstype* uinf, const dstype* param,
                        dstype t, int ng) {
    constexpr int nd = M::nd, ncu = M::ncu, ncw = M::ncw;
    constexpr int Nq = ncu * (1 + nd);
    Kokkos::parallel_for("Flux", ng, KOKKOS_LAMBDA(size_t i) {
        double x[nd], uq[Nq], w[ncw];
        for (int k = 0; k < nd;  ++k) x[k]  = xdg[k * ng + i];
        for (int k = 0; k < Nq;  ++k) uq[k] = udg[k * ng + i];
        for (int k = 0; k < ncw; ++k) w[k]  = wdg[k * ng + i];
        double fi[ncu * nd];
        M::flux(fi, x, uq, w, param, t);
        for (int k = 0; k < ncu * nd; ++k) f[k * ng + i] = fi[k];
    });
}
} // namespace exasim
```
