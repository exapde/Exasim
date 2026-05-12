#include "model_binding.hpp"

namespace exasim::modeldispatch::builtin_shared_library_legacy {

#include "../BuiltIn/BuiltinModelDrivers.cpp"

} // namespace exasim::modeldispatch::builtin_shared_library_legacy

namespace exasim {

ModelBinding makeBuiltinSharedLibraryModelBinding()
{
    static const ModelOps ops = [] {
        ModelOps value;
        value.flux = static_cast<ModelOps::ElementDriverFn>(&modeldispatch::builtin_shared_library_legacy::FluxDriver);
        value.source = static_cast<ModelOps::ElementDriverFn>(&modeldispatch::builtin_shared_library_legacy::SourceDriver);
        value.sourcew = static_cast<ModelOps::ElementDriverFn>(&modeldispatch::builtin_shared_library_legacy::SourcewDriver);
        value.output = static_cast<ModelOps::GlobalElementDriverFn>(&modeldispatch::builtin_shared_library_legacy::OutputDriver);
        value.monitor = static_cast<ModelOps::MonitorDriverFn>(&modeldispatch::builtin_shared_library_legacy::MonitorDriver);
        value.avfield = static_cast<ModelOps::GlobalElementDriverFn>(&modeldispatch::builtin_shared_library_legacy::AvfieldDriver);
        value.eos = static_cast<ModelOps::ElementDriverFn>(&modeldispatch::builtin_shared_library_legacy::EosDriver);
        value.eosdu = static_cast<ModelOps::ElementDriverFn>(&modeldispatch::builtin_shared_library_legacy::EosduDriver);
        value.eosdw = static_cast<ModelOps::ElementDriverFn>(&modeldispatch::builtin_shared_library_legacy::EosdwDriver);
        value.tdfunc = static_cast<ModelOps::ElementDriverFn>(&modeldispatch::builtin_shared_library_legacy::TdfuncDriver);
        value.fhatLDG = static_cast<ModelOps::FaceCoupledDriverFn>(&modeldispatch::builtin_shared_library_legacy::FhatDriver);
        value.fbouLDG = static_cast<ModelOps::BoundaryDriverFn>(&modeldispatch::builtin_shared_library_legacy::FbouDriver);
        value.uhat = static_cast<ModelOps::UhatDriverFn>(&modeldispatch::builtin_shared_library_legacy::UhatDriver);
        value.ubou = static_cast<ModelOps::BoundaryDriverFn>(&modeldispatch::builtin_shared_library_legacy::UbouDriver);
        value.initodg = static_cast<ModelOps::InitDriverFn>(&modeldispatch::builtin_shared_library_legacy::InitodgDriver);
        value.initq = static_cast<ModelOps::InitDriverFn>(&modeldispatch::builtin_shared_library_legacy::InitqDriver);
        value.initudg = static_cast<ModelOps::InitDriverFn>(&modeldispatch::builtin_shared_library_legacy::InitudgDriver);
        value.initu = static_cast<ModelOps::InitDriverFn>(&modeldispatch::builtin_shared_library_legacy::InituDriver);
        value.initwdg = static_cast<ModelOps::InitDriverFn>(&modeldispatch::builtin_shared_library_legacy::InitwdgDriver);
        value.cpuInitodg = static_cast<ModelOps::InitDriverFn>(&modeldispatch::builtin_shared_library_legacy::cpuInitodgDriver);
        value.cpuInitq = static_cast<ModelOps::InitDriverFn>(&modeldispatch::builtin_shared_library_legacy::cpuInitqDriver);
        value.cpuInitudg = static_cast<ModelOps::InitDriverFn>(&modeldispatch::builtin_shared_library_legacy::cpuInitudgDriver);
        value.cpuInitu = static_cast<ModelOps::InitDriverFn>(&modeldispatch::builtin_shared_library_legacy::cpuInituDriver);
        value.cpuInitwdg = static_cast<ModelOps::InitDriverFn>(&modeldispatch::builtin_shared_library_legacy::cpuInitwdgDriver);
        value.fluxJac = static_cast<ModelOps::FluxJacDriverFn>(&modeldispatch::builtin_shared_library_legacy::FluxDriver);
        value.sourceJac = static_cast<ModelOps::ElementJacDriverFn>(&modeldispatch::builtin_shared_library_legacy::SourceDriver);
        value.sourcewJac = static_cast<ModelOps::ElementJacDriverFn>(&modeldispatch::builtin_shared_library_legacy::SourcewDriver);
        value.sourcewWdg = static_cast<ModelOps::SourcewWdgDriverFn>(&modeldispatch::builtin_shared_library_legacy::SourcewDriver);
        value.eosJac = static_cast<ModelOps::ElementJacDriverFn>(&modeldispatch::builtin_shared_library_legacy::EosDriver);
        value.fbouHDGJac = static_cast<ModelOps::BoundaryJacDriverFn>(&modeldispatch::builtin_shared_library_legacy::FbouDriver);
        value.fbouHDG = static_cast<ModelOps::BoundaryStateDriverFn>(&modeldispatch::builtin_shared_library_legacy::FbouDriver);
        value.fintJac = static_cast<ModelOps::BoundaryJacDriverFn>(&modeldispatch::builtin_shared_library_legacy::FintDriver);
        value.fintState = static_cast<ModelOps::BoundaryStateDriverFn>(&modeldispatch::builtin_shared_library_legacy::FintDriver);
        value.fextJac = static_cast<ModelOps::BoundaryExternalJacDriverFn>(&modeldispatch::builtin_shared_library_legacy::FextDriver);
        value.fextState = static_cast<ModelOps::BoundaryExternalStateDriverFn>(&modeldispatch::builtin_shared_library_legacy::FextDriver);
        value.fhatHDGJac = static_cast<ModelOps::FaceJacDriverFn>(&modeldispatch::builtin_shared_library_legacy::FhatDriver);
        value.fhatHDG = static_cast<ModelOps::FaceStateDriverFn>(&modeldispatch::builtin_shared_library_legacy::FhatDriver);
        value.visScalars = static_cast<ModelOps::ElementDriverFn>(&modeldispatch::builtin_shared_library_legacy::VisScalarsDriver);
        value.visVectors = static_cast<ModelOps::ElementDriverFn>(&modeldispatch::builtin_shared_library_legacy::VisVectorsDriver);
        value.visTensors = static_cast<ModelOps::ElementDriverFn>(&modeldispatch::builtin_shared_library_legacy::VisTensorsDriver);
        value.qoiVolume = static_cast<ModelOps::ElementDriverFn>(&modeldispatch::builtin_shared_library_legacy::QoIvolumeDriver);
        value.qoiBoundary = static_cast<ModelOps::BoundaryDriverFn>(&modeldispatch::builtin_shared_library_legacy::QoIboundaryDriver);
        return value;
    }();

    ModelBinding binding;
    binding.ops = &ops;
    binding.provider = ModelProvider::BuiltinSharedLibrary;
    return binding;
}

} // namespace exasim
