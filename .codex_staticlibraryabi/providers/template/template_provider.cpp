#include "driver_abi.hpp"

struct MyModel {
    static double Add(double a, double b)
    {
        return a + b;
    }
};

template <class Model>
double AddTemplate(double a, double b)
{
    return Model::Add(a, b);
}

namespace {

double Add(double a, double b)
{
    return AddTemplate<MyModel>(a, b);
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
