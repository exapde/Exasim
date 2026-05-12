#pragma once

#include <cstdint>

namespace staticabi {

inline constexpr std::uint32_t kDriverAbiVersion = 1;

struct DriverABI {
    std::uint32_t abi_version;
    double (*add)(double a, double b);
};

using CreateDriverABIFn = const DriverABI* (*)();

} // namespace staticabi
