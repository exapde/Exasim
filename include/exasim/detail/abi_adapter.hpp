// SPDX-License-Identifier: see LICENSE
//
// <exasim/detail/abi_adapter.hpp> — AbiAdapter marker type.
//
// The four FEM class templates (CDiscretization, CSolver,
// CPreconditioner, CSolution) all default to `M = AbiAdapter` to
// preserve source-level compatibility with backend/Main/main.cpp,
// which constructs them without specifying a Model. With M defaulted
// to AbiAdapter, the templates pick up the legacy code path that
// calls into the libpdemodel.hpp ABI symbols (KokkosFlux, HdgFbou,
// …) — the same path Exasim has always used.
//
// User code that instantiates CSolution<MyModel> bypasses the ABI;
// the templated *Driver<M> wrappers in <exasim/drivers.hpp> route
// through exasim::*_kernel<M> and the user's pointwise math.
//
// This type is empty by design — it carries no compile-time shape.
// The legacy code path inside the FEM classes does not call any
// kernel<M> templates against it; its only role is to identify
// "use the libpdemodel.hpp ABI" at compile time.

#pragma once

namespace exasim::detail {

struct AbiAdapter {};

} // namespace exasim::detail
