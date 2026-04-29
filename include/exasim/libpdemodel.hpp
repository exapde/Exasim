// <exasim/libpdemodel.hpp> — the model ABI that Exasim's runtime calls into.
//
// Users of Exasim-as-a-library implement these ~30 free functions in their own
// glue.cpp (typically as one-liners that delegate to <exasim/kernels/*.hpp>
// templated kernels parameterized on a user-defined Model struct).
//
// text2code keeps using the same ABI: it generates one .so that fills these
// symbols. The runtime does not care which producer fulfilled the contract.
//
// Transitional shim. After Phase 1.2 the canonical declarations move here.
#pragma once
#include "../../backend/Model/libpdemodel.hpp"
