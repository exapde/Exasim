// SPDX-License-Identifier: see LICENSE
//
// <exasim/detail/driver_dispatch.hpp>
//
// EXASIM_DRIVER_CALL(Name, args...) — used inside the templated FEM
// classes (CDiscretization<M> & friends) and the templated chain
// free functions (residual.hpp, uequation.hpp, …) to dispatch driver
// calls based on the user Model type M:
//
//   - If `M == exasim::detail::AbiAdapter` (legacy build): call the
//     global non-templated `::Name(args...)`, which is defined by
//     backend/Model/{KokkosDrivers,ModelDrivers}.cpp and ultimately
//     dispatches through the libpdemodel.hpp ABI symbols.
//
//   - Otherwise (user code instantiating CSolution<MyModel>): call
//     `::exasim::Name<M>(args...)`, which is defined in
//     <exasim/drivers.hpp> and routes through the templated kernels
//     in <exasim/kernels/*.hpp>, ultimately invoking the user's
//     pointwise math (`M::flux`, `M::flux_jac_uq`, …).
//
// The macro must be expanded inside a context where `M` is in scope
// (a `template <class M>` class member or free function). `if constexpr`
// requires C++17.
//
// We use a macro rather than a wrapper function template because the
// driver argument lists are heterogeneous (mesh/master/app refs +
// scalars + variable counts) and there are many overload sets where
// argument-dependent template deduction would be brittle. The macro
// is a thin pass-through.

#pragma once

#include <type_traits>

#include "abi_adapter.hpp"

#define EXASIM_DRIVER_CALL(Name, ...)                                      \
    do {                                                                   \
        if constexpr (std::is_same_v<M, exasim::detail::AbiAdapter>) {     \
            Name(__VA_ARGS__);                                             \
        } else {                                                           \
            exasim::Name<M>(__VA_ARGS__);                                  \
        }                                                                  \
    } while (0)
