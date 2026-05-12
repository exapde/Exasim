#include "core.hpp"

#include <iostream>

extern "C" const staticabi::DriverABI* GetDriverABI();

int main()
{
    const double value = staticabi::RunAdd(*GetDriverABI(), 3.0, 4.0);

    std::cout << "dynamic provider: 3 + 4 = " << value << '\n';
    return 0;
}
