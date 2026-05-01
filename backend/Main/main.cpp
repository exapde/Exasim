/**
 * @file main.cpp
 * @brief Legacy `cput2cEXASIM` / `cpumpit2cEXASIM` / `cpuEXASIM`
 *        entry point.
 *
 * Calls `exasim::run<exasim::detail::AbiAdapter>(argc, argv)`, which
 * dispatches every kernel call through the legacy `libpdemodel.hpp`
 * ABI (the `if constexpr` branch in `EXASIM_DRIVER_CALL` that chooses
 * the global non-templated symbol).
 *
 * The body that used to live here moved to `<exasim/run.hpp>` so that
 * user code can do
 *
 *     exasim::run<MyModel>(argc, argv)
 *
 * with the same MPI / Kokkos / DIRK / multi-domain logic.
 *
 * Author: Exasim Team
 * Date: 2025
 */

#include "../../include/exasim/run.hpp"

int main(int argc, char** argv) {
    return exasim::run<exasim::detail::AbiAdapter>(argc, argv);
}
