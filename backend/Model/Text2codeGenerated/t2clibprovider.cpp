#include "../ModelDispatch/driver_abi.h"

#include <stdexcept>

extern "C" const ExasimDriverABI* GetText2CodeExasimDriverABI();

const ExasimDriverABI& getText2codeLibraryExasimDriverABI()
{
    const ExasimDriverABI* abi = GetText2CodeExasimDriverABI();

    if (!abi)
        throw std::runtime_error("GetText2CodeExasimDriverABI returned a null ABI pointer");
    if (abi->abi_version != kExasimDriverABIVersion)
        throw std::runtime_error("Text2Code model library ABI version does not match Exasim");
    if (abi->struct_size != sizeof(ExasimDriverABI))
        throw std::runtime_error("Text2Code model library ABI struct size does not match Exasim");

    return *abi;
}
