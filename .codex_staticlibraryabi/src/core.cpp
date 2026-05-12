#include "core.hpp"

#include <stdexcept>

namespace staticabi {

namespace {

void ValidateDriver(const DriverABI* driver)
{
    if (driver == nullptr) {
        throw std::runtime_error("Provider returned a null DriverABI pointer.");
    }
    if (driver->abi_version != kDriverAbiVersion) {
        throw std::runtime_error("Provider ABI version does not match the core library.");
    }
    if (driver->add == nullptr) {
        throw std::runtime_error("Provider did not populate DriverABI::add.");
    }
}

} // namespace

double RunAdd(const DriverABI& driver, double a, double b)
{
    ValidateDriver(&driver);
    return driver.add(a, b);
}

} // namespace staticabi
