/**
 * @file model_adapter.hpp
 * @brief Shared ModelDispatch factory interface.
 *
 * This header exposes the non-template provider factories used by the runtime
 * dispatch layer and re-exports the user-defined source-model template helpers
 * implemented under backend/Model/UserDefined/.
 */
#ifndef __EXASIM_MODEL_DISPATCH_MODEL_ADAPTER_HPP__
#define __EXASIM_MODEL_DISPATCH_MODEL_ADAPTER_HPP__

#include "../../../include/exasim/detail/abi_adapter.hpp"
#include "model_binding.hpp"

namespace exasim {

// Source-compiled model providers
// -------------------------------
// Generated-source models bind through a per-model adapter translation unit
// that is compiled with the generated backend/Model sources, not into the
// reusable core static library. That adapter should include:
//   backend/Model/KokkosDrivers.cpp
ModelBinding makeGeneratedSourceModelBinding();

// Shared-library model providers
// ------------------------------
// text2code / libpdemodel shared-library path:
//   backend/Model/ModelDrivers.cpp
ModelBinding makeSharedLibraryModelBinding();

// builtin PDE-model shared-library path:
//   backend/Model/BuiltIn/BuiltinModelDrivers.cpp
ModelBinding makeBuiltinSharedLibraryModelBinding();

// Compatibility hook retained during the transition; the implementation can
// forward to makeSharedLibraryModelBinding() once the legacy path is fully
// absorbed by ModelDispatch.
ModelBinding makeAbiModelBinding();

} // namespace exasim

#include "../UserDefined/user_defined_model_adapter.hpp"

#endif
