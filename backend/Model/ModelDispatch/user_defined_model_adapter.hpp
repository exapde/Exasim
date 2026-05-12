/**
 * @file user_defined_model_adapter.hpp
 * @brief Template helpers for source-compiled user-defined Exasim models.
 *
 * This header implements the model-specific source adapter path used for
 * user-defined models that are compiled with an application and linked against
 * the non-templated backend static library.
 */
#ifndef __EXASIM_USER_DEFINED_MODEL_ADAPTER_HPP__
#define __EXASIM_USER_DEFINED_MODEL_ADAPTER_HPP__

#include <type_traits>

#include "../../../include/exasim/kernels/init.hpp"
#include "../../../include/exasim/model.hpp"
#include "../ModelDispatch/model_binding.hpp"

namespace exasim {

namespace detail {

#define EXASIM_MODEL_OPTIONAL_INT_GETTER(name)                                      \
template <class, class = void> struct has_##name : std::false_type {};              \
template <class M>                                                                  \
struct has_##name<M, std::void_t<decltype(M::name)>> : std::true_type {};           \
template <class M>                                                                  \
inline constexpr Int get_##name() {                                                 \
    if constexpr (has_##name<M>::value)                                             \
        return static_cast<Int>(M::name);                                           \
    else                                                                            \
        return 0;                                                                   \
}

EXASIM_MODEL_OPTIONAL_INT_GETTER(nexternalparam);
EXASIM_MODEL_OPTIONAL_INT_GETTER(nce);
EXASIM_MODEL_OPTIONAL_INT_GETTER(ncuext);
EXASIM_MODEL_OPTIONAL_INT_GETTER(ncm);
EXASIM_MODEL_OPTIONAL_INT_GETTER(nsca);
EXASIM_MODEL_OPTIONAL_INT_GETTER(nvec);
EXASIM_MODEL_OPTIONAL_INT_GETTER(nten);
EXASIM_MODEL_OPTIONAL_INT_GETTER(nsurf);
EXASIM_MODEL_OPTIONAL_INT_GETTER(nvqoi);

#undef EXASIM_MODEL_OPTIONAL_INT_GETTER

template <class M>
inline void user_defined_initq_adapter(dstype* f, const dstype* xg,
                                       appstruct& app, Int ncx, Int nc,
                                       Int npe, Int ne, Int /*backend*/)
{
    initq_kernel<M>(f, xg, app.uinf, app.physicsparam,
                    /*modelnumber=*/0, npe * ne, ncx, nc, npe, ne);
}

template <class M>
inline void user_defined_cpuInitq_adapter(dstype* f, const dstype* xg,
                                          appstruct& app, Int ncx, Int nc,
                                          Int npe, Int ne, Int /*backend*/)
{
    for (Int elem = 0; elem < ne; ++elem) {
        for (Int j = 0; j < npe; ++j) {
            double x[M::nd];
            for (int k = 0; k < M::nd; ++k) {
                x[k] = xg[j + npe * k + npe * ncx * elem];
            }

            double q[M::ncu * M::nd];
            M::cpuinitq(q, x, app.uinf, app.physicsparam);

            for (int k = 0; k < M::ncu * M::nd; ++k) {
                f[j + npe * k + npe * nc * elem] = q[k];
            }
        }
    }
}

template <class M>
inline void user_defined_initu_adapter(dstype* f, const dstype* xg,
                                       appstruct& app, Int ncx, Int nc,
                                       Int npe, Int ne, Int /*backend*/)
{
    initu_kernel<M>(f, xg, app.uinf, app.physicsparam,
                    /*modelnumber=*/0, npe * ne, ncx, nc, npe, ne);
}

template <class M>
inline void user_defined_cpuInitu_adapter(dstype* f, const dstype* xg,
                                          appstruct& app, Int ncx, Int nc,
                                          Int npe, Int ne, Int /*backend*/)
{
    for (Int elem = 0; elem < ne; ++elem) {
        for (Int j = 0; j < npe; ++j) {
            double x[M::nd];
            for (int k = 0; k < M::nd; ++k) {
                x[k] = xg[j + npe * k + npe * ncx * elem];
            }

            double u[M::ncu];
            M::cpuinitu(u, x, app.uinf, app.physicsparam);

            for (int k = 0; k < M::ncu; ++k) {
                f[j + npe * k + npe * nc * elem] = u[k];
            }
        }
    }
}

template <class M>
inline void user_defined_initudg_adapter(dstype* f, const dstype* xg,
                                         appstruct& app, Int ncx, Int nc,
                                         Int npe, Int ne, Int /*backend*/)
{
    initudg_kernel<M>(f, xg, app.uinf, app.physicsparam,
                      /*modelnumber=*/0, npe * ne, ncx, nc, npe, ne);
}

template <class M>
inline void user_defined_cpuInitudg_adapter(dstype* f, const dstype* xg,
                                            appstruct& app, Int ncx, Int nc,
                                            Int npe, Int ne, Int /*backend*/)
{
    for (Int elem = 0; elem < ne; ++elem) {
        for (Int j = 0; j < npe; ++j) {
            double x[M::nd];
            for (int k = 0; k < M::nd; ++k) {
                x[k] = xg[j + npe * k + npe * ncx * elem];
            }

            double udg[M::Nq];
            M::cpuinitudg(udg, x, app.uinf, app.physicsparam);

            for (int k = 0; k < M::Nq; ++k) {
                f[j + npe * k + npe * nc * elem] = udg[k];
            }
        }
    }
}

template <class M>
inline void user_defined_initodg_adapter(dstype* f, const dstype* xg,
                                         appstruct& app, Int ncx, Int nco,
                                         Int npe, Int ne, Int /*backend*/)
{
    initodg_kernel<M>(f, xg, app.uinf, app.physicsparam,
                      /*modelnumber=*/0, npe * ne, ncx, nco, npe, ne);
}

template <class M>
inline void user_defined_cpuInitodg_adapter(dstype* f, const dstype* xg,
                                            appstruct& app, Int ncx, Int nco,
                                            Int npe, Int ne, Int /*backend*/)
{
    for (Int elem = 0; elem < ne; ++elem) {
        for (Int j = 0; j < npe; ++j) {
            for (Int k = 0; k < nco; ++k) {
                f[j + npe * k + npe * nco * elem] = 0.0;
            }

            if constexpr (M::nco > 0) {
                double x[M::nd];
                for (int k = 0; k < M::nd; ++k) {
                    x[k] = xg[j + npe * k + npe * ncx * elem];
                }

                double odg[M::nco];
                for (int k = 0; k < M::nco; ++k) {
                    odg[k] = 0.0;
                }

                M::cpuinitodg(odg, x, app.uinf, app.physicsparam);

                const Int ncopy = (nco < static_cast<Int>(M::nco)) ? nco : static_cast<Int>(M::nco);
                for (Int k = 0; k < ncopy; ++k) {
                    f[j + npe * k + npe * nco * elem] = odg[k];
                }
            }
        }
    }
}

template <class M>
inline void user_defined_initwdg_adapter(dstype* f, const dstype* xg,
                                         appstruct& app, Int ncx, Int ncw,
                                         Int npe, Int ne, Int /*backend*/)
{
    initwdg_kernel<M>(f, xg, app.uinf, app.physicsparam,
                      /*modelnumber=*/0, npe * ne, ncx, ncw, npe, ne);
}

template <class M>
inline void user_defined_cpuInitwdg_adapter(dstype* f, const dstype* xg,
                                            appstruct& app, Int ncx, Int ncw,
                                            Int npe, Int ne, Int /*backend*/)
{
    if constexpr (M::ncw > 0) {
        for (Int elem = 0; elem < ne; ++elem) {
            for (Int j = 0; j < npe; ++j) {
                double x[M::nd];
                for (int k = 0; k < M::nd; ++k) {
                    x[k] = xg[j + npe * k + npe * ncx * elem];
                }

                double w[M::ncw];
                M::cpuinitwdg(w, x, app.uinf, app.physicsparam);

                for (int k = 0; k < M::ncw; ++k) {
                    f[j + npe * k + npe * ncw * elem] = w[k];
                }
            }
        }
    } else {
        (void)f; (void)xg; (void)app; (void)ncx; (void)ncw; (void)npe; (void)ne;
    }
}

} // namespace detail

template <class M>
const ModelOps& getModelOps() {
    static const ModelOps ops = [] {
        ModelOps value;
        value.initodg = &detail::user_defined_initodg_adapter<M>;
        value.initq = &detail::user_defined_initq_adapter<M>;
        value.initudg = &detail::user_defined_initudg_adapter<M>;
        value.initu = &detail::user_defined_initu_adapter<M>;
        value.initwdg = &detail::user_defined_initwdg_adapter<M>;
        value.cpuInitodg = &detail::user_defined_cpuInitodg_adapter<M>;
        value.cpuInitq = &detail::user_defined_cpuInitq_adapter<M>;
        value.cpuInitu = &detail::user_defined_cpuInitu_adapter<M>;
        value.cpuInitudg = &detail::user_defined_cpuInitudg_adapter<M>;
        value.cpuInitwdg = &detail::user_defined_cpuInitwdg_adapter<M>;
        return value;
    }();
    return ops;
}

template <class M>
ModelBinding makeModelBinding(ModelProvider provider = ModelProvider::Unknown) {
    static_assert(is_model_v<M>,
                  "makeModelBinding<M>() requires M to satisfy the Exasim model contract.");

    ModelBinding binding;
    binding.ops = &getModelOps<M>();
    binding.provider = provider;
    binding.nd = static_cast<Int>(M::nd);
    binding.ncu = static_cast<Int>(M::ncu);
    binding.ncq = static_cast<Int>(M::ncu * M::nd);
    binding.ncw = static_cast<Int>(M::ncw);
    binding.nco = static_cast<Int>(M::nco);
    binding.nch = static_cast<Int>(M::ncu);
    binding.ncx = static_cast<Int>(M::nd);
    binding.nce = detail::get_nce<M>();
    binding.ncuext = detail::get_ncuext<M>();
    binding.ncm = detail::get_ncm<M>();
    binding.nsca = detail::get_nsca<M>();
    binding.nvec = detail::get_nvec<M>();
    binding.nten = detail::get_nten<M>();
    binding.nsurf = detail::get_nsurf<M>();
    binding.nvqoi = detail::get_nvqoi<M>();
    binding.nparam = static_cast<Int>(M::nparam);
    binding.nexternalparam = detail::get_nexternalparam<M>();
    binding.ncuq = static_cast<Int>(M::Nq);
    binding.nc = binding.ncuq;
    return binding;
}

template <class M>
ModelBinding makeUserDefinedSourceModelBinding() {
    return makeModelBinding<M>(ModelProvider::UserDefinedSource);
}

} // namespace exasim

#endif
