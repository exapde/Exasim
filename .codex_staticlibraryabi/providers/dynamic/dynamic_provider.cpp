#include "driver_abi.hpp"

namespace {

double Add(double a, double b)
{
    return a + b;
}

} // namespace

extern "C" const staticabi::DriverABI* GetDriverABI()
{
    static const staticabi::DriverABI abi{
        staticabi::kDriverAbiVersion,
        &Add
    };
    return &abi;
}
