#include "model_binding.hpp"

namespace exasim::frontend_generated_legacy {

#include "../FrontendGenerated/KokkosDrivers.cpp"

} // namespace exasim::frontend_generated_legacy

namespace exasim {

ModelBinding makeGeneratedSourceModelBinding()
{
    static const ModelOps ops = [] {
        ModelOps value;
        value.flux = static_cast<ModelOps::ElementDriverFn>(&frontend_generated_legacy::FluxDriver);
        value.source = static_cast<ModelOps::ElementDriverFn>(&frontend_generated_legacy::SourceDriver);
        value.sourcew = static_cast<ModelOps::ElementDriverFn>(&frontend_generated_legacy::SourcewDriver);
        value.output = static_cast<ModelOps::GlobalElementDriverFn>(&frontend_generated_legacy::OutputDriver);
        value.monitor = static_cast<ModelOps::MonitorDriverFn>(&frontend_generated_legacy::MonitorDriver);
        value.avfield = static_cast<ModelOps::GlobalElementDriverFn>(&frontend_generated_legacy::AvfieldDriver);
        value.eos = static_cast<ModelOps::ElementDriverFn>(&frontend_generated_legacy::EosDriver);
        value.eosdu = static_cast<ModelOps::ElementDriverFn>(&frontend_generated_legacy::EosduDriver);
        value.eosdw = static_cast<ModelOps::ElementDriverFn>(&frontend_generated_legacy::EosdwDriver);
        value.tdfunc = static_cast<ModelOps::ElementDriverFn>(&frontend_generated_legacy::TdfuncDriver);
        value.fhatLDG = static_cast<ModelOps::FaceCoupledDriverFn>(&frontend_generated_legacy::FhatDriver);
        value.fbouLDG = static_cast<ModelOps::BoundaryDriverFn>(&frontend_generated_legacy::FbouDriver);
        value.uhat = static_cast<ModelOps::UhatDriverFn>(&frontend_generated_legacy::UhatDriver);
        value.ubou = static_cast<ModelOps::BoundaryDriverFn>(&frontend_generated_legacy::UbouDriver);
        value.initodg = static_cast<ModelOps::InitDriverFn>(&frontend_generated_legacy::InitodgDriver);
        value.initq = static_cast<ModelOps::InitDriverFn>(&frontend_generated_legacy::InitqDriver);
        value.initudg = static_cast<ModelOps::InitDriverFn>(&frontend_generated_legacy::InitudgDriver);
        value.initu = static_cast<ModelOps::InitDriverFn>(&frontend_generated_legacy::InituDriver);
        value.initwdg = static_cast<ModelOps::InitDriverFn>(&frontend_generated_legacy::InitwdgDriver);
        value.cpuInitodg = static_cast<ModelOps::InitDriverFn>(&frontend_generated_legacy::cpuInitodgDriver);
        value.cpuInitq = static_cast<ModelOps::InitDriverFn>(&frontend_generated_legacy::cpuInitqDriver);
        value.cpuInitudg = static_cast<ModelOps::InitDriverFn>(&frontend_generated_legacy::cpuInitudgDriver);
        value.cpuInitu = static_cast<ModelOps::InitDriverFn>(&frontend_generated_legacy::cpuInituDriver);
        value.cpuInitwdg = static_cast<ModelOps::InitDriverFn>(&frontend_generated_legacy::cpuInitwdgDriver);
        value.fluxJac = static_cast<ModelOps::FluxJacDriverFn>(&frontend_generated_legacy::FluxDriver);
        value.sourceJac = static_cast<ModelOps::ElementJacDriverFn>(&frontend_generated_legacy::SourceDriver);
        value.sourcewJac = static_cast<ModelOps::ElementJacDriverFn>(&frontend_generated_legacy::SourcewDriver);
        value.sourcewWdg = static_cast<ModelOps::SourcewWdgDriverFn>(&frontend_generated_legacy::SourcewDriver);
        value.eosJac = static_cast<ModelOps::ElementJacDriverFn>(&frontend_generated_legacy::EosDriver);
        value.fbouHDGJac = static_cast<ModelOps::BoundaryJacDriverFn>(&frontend_generated_legacy::FbouDriver);
        value.fbouHDG = static_cast<ModelOps::BoundaryStateDriverFn>(&frontend_generated_legacy::FbouDriver);
        value.fintJac = static_cast<ModelOps::BoundaryJacDriverFn>(&frontend_generated_legacy::FintDriver);
        value.fintState = static_cast<ModelOps::BoundaryStateDriverFn>(&frontend_generated_legacy::FintDriver);
        value.fextJac = static_cast<ModelOps::BoundaryExternalJacDriverFn>(&frontend_generated_legacy::FextDriver);
        value.fextState = static_cast<ModelOps::BoundaryExternalStateDriverFn>(&frontend_generated_legacy::FextDriver);
        value.fhatHDGJac = static_cast<ModelOps::FaceJacDriverFn>(&frontend_generated_legacy::FhatDriver);
        value.fhatHDG = static_cast<ModelOps::FaceStateDriverFn>(&frontend_generated_legacy::FhatDriver);
        value.visScalars = static_cast<ModelOps::ElementDriverFn>(&frontend_generated_legacy::VisScalarsDriver);
        value.visVectors = static_cast<ModelOps::ElementDriverFn>(&frontend_generated_legacy::VisVectorsDriver);
        value.visTensors = static_cast<ModelOps::ElementDriverFn>(&frontend_generated_legacy::VisTensorsDriver);
        value.qoiVolume = static_cast<ModelOps::ElementDriverFn>(&frontend_generated_legacy::QoIvolumeDriver);
        value.qoiBoundary = static_cast<ModelOps::BoundaryDriverFn>(&frontend_generated_legacy::QoIboundaryDriver);
        return value;
    }();

    ModelBinding binding;
    binding.ops = &ops;
    binding.provider = ModelProvider::GeneratedSource;
    return binding;
}

} // namespace exasim
