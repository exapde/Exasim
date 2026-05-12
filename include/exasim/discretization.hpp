// SPDX-License-Identifier: see LICENSE
//
// <exasim/discretization.hpp> — discretization tags shared by models and runtime bindings.

#pragma once

namespace exasim {

// Discretization tag — published as `Self::disc` by the user model so
// kernels can `if constexpr` on it. Defaults to LDG (Local DG) since
// it requires no Jacobians.
enum class Discretization {
    LDG,   // Local DG: only flux/source/initu values needed.
    HDG    // Hybridized DG: requires hand-written pointwise Jacobians.
};

} // namespace exasim
