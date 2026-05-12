/**
 * @file model_binding.hpp
 * @brief Runtime model binding metadata for the hybrid Exasim architecture.
 *
 * ModelBinding is the constructor-time object passed into CSolution and
 * CDiscretization. It couples a model's dispatch table (ModelOps) with the
 * compile-time metadata the non-templated backend still needs at runtime.
 */
#ifndef __EXASIM_MODEL_DISPATCH_MODEL_BINDING_HPP__
#define __EXASIM_MODEL_DISPATCH_MODEL_BINDING_HPP__

#include "model_ops.hpp"

namespace exasim {

enum class ModelProvider {
    Unknown = 0,
    GeneratedSource,
    SharedLibrary,
    BuiltinSharedLibrary,
    UserDefinedSource
};

struct ModelBinding {
    const ModelOps* ops = nullptr;
    ModelProvider provider = ModelProvider::Unknown;

    // Standard Exasim component counts. Keep these aligned with the
    // backend `commonstruct` vocabulary so constructor-time binding can
    // seed the non-templated runtime without ad hoc renaming.
    Int nc = 0;
    Int nd = 0;
    Int ncu = 0;
    Int ncq = 0;
    Int ncw = 0;
    Int nco = 0;
    Int nch = 0;
    Int ncx = 0;
    Int nce = 0;
    Int ncuext = 0;
    Int ncm = 0;
    Int nsca = 0;
    Int nvec = 0;
    Int nten = 0;
    Int nsurf = 0;
    Int nvqoi = 0;
    Int nparam = 0;
    Int nexternalparam = 0;
    Int ncuq = 0;

    // Runtime discretization / mesh counts. These are not model constants;
    // they must be populated after app/master/mesh are read, before the
    // corresponding dispatch entries are used.
    Int npe = 0;
    Int nge = 0;
    Int npf = 0;
    Int ngf = 0;
    Int ne = 0;
    Int N = 0;

    [[nodiscard]] bool isBound() const {
        return ops != nullptr;
    }
};

} // namespace exasim

#endif
