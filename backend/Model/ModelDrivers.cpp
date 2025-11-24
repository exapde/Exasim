#include "modeltemplate.hpp"
#include "libpdemodel.hpp"
#include "modeltemplate.hpp"
#include "libpdemodel.hpp"

// ======================================================
// Volume flux / source (no derivatives, non-HDG)
// ======================================================

template <typename Model = DefaultModel>
void FluxDriver(dstype* f, const dstype* xg, const dstype* udg,
                const dstype* odg, const dstype* wdg,
                meshstruct &mesh, masterstruct &master,
                appstruct &app, solstruct &sol, tempstruct &temp,
                commonstruct &common, Int nge, Int e1, Int e2, Int backend)
{
  Int nc  = common.nc;
  Int ncu = common.ncu;
  Int ncw = common.ncw;
  Int nco = common.nco;
  Int ncx = common.ncx;
  Int nd  = common.nd;
  Int numPoints = nge * (e2 - e1);
  dstype time = common.time;

  typename Model::FluxFn flux{};
  flux(f, xg, udg, odg, wdg,
       app.uinf, app.physicsparam, time,
       common.modelnumber, numPoints,
       nc, ncu, nd, ncx, nco, ncw);
}

template <typename Model = DefaultModel>
void SourceDriver(dstype* f, const dstype* xg, const dstype* udg,
                  const dstype* odg, const dstype* wdg,
                  meshstruct &mesh, masterstruct &master,
                  appstruct &app, solstruct &sol, tempstruct &temp,
                  commonstruct &common, Int nge, Int e1, Int e2, Int backend)
{
  Int nc  = common.nc;
  Int ncu = common.ncu;
  Int ncw = common.ncw;
  Int nco = common.nco;
  Int ncx = common.ncx;
  Int nd  = common.nd;
  Int numPoints = nge * (e2 - e1);
  dstype time = common.time;

  typename Model::SourceFn source{};
  source(f, xg, udg, odg, wdg,
         app.uinf, app.physicsparam, time,
         common.modelnumber, numPoints,
         nc, ncu, nd, ncx, nco, ncw);
}

template <typename Model = DefaultModel>
void SourcewDriver(dstype* f, const dstype* xg, const dstype* udg,
                   const dstype* odg, const dstype* wdg,
                   meshstruct &mesh, masterstruct &master,
                   appstruct &app, solstruct &sol, tempstruct &temp,
                   commonstruct &common, Int npe, Int e1, Int e2, Int backend)
{
  Int nc  = common.nc;
  Int ncu = common.ncu;
  Int ncw = common.ncw;
  Int nco = common.nco;
  Int ncx = common.ncx;
  Int nd  = common.nd;
  Int ne  = e2 - e1;
  Int numPoints = npe * ne;
  dstype time = common.time;

  typename Model::SourcewFn sourcew{};
  sourcew(f, xg, udg, odg, wdg,
          app.uinf, app.physicsparam, time, common.modelnumber,
          numPoints, nc, ncu, nd, ncx, nco, ncw,
          ncw, npe, ne);
}

// ======================================================
// Output / monitor / Avfield / EoS (non-HDG, no derivatives)
// ======================================================

template <typename Model = DefaultModel>
void OutputDriver(dstype* f, const dstype* xg, const dstype* udg,
                  const dstype* odg, const dstype* wdg,
                  meshstruct &mesh, masterstruct &master,
                  appstruct &app, solstruct &sol, tempstruct &temp,
                  commonstruct &common, Int backend)
{
  Int nc  = common.nc;
  Int ncu = common.ncu;
  Int ncw = common.ncw;
  Int nco = common.nco;
  Int ncx = common.ncx;
  Int nce = common.nce;
  Int nd  = common.nd;
  Int npe = common.npe;
  Int ne  = common.ne;
  Int numPoints = npe * ne;
  dstype time = common.time;

  typename Model::OutputFn output{};
  output(f, xg, udg, odg, wdg,
         app.uinf, app.physicsparam, time,
         common.modelnumber, numPoints,
         nc, ncu, nd, ncx, nco, ncw, nce, npe, ne);
}

template <typename Model = DefaultModel>
void MonitorDriver(dstype* f, Int nc_sol,
                   const dstype* xg, const dstype* udg,
                   const dstype* odg, const dstype* wdg,
                   meshstruct &mesh, masterstruct &master,
                   appstruct &app, solstruct &sol, tempstruct &temp,
                   commonstruct &common, Int backend)
{
  Int ncu = common.ncu;
  Int ncw = common.ncw;
  Int nco = common.nco;
  Int ncx = common.ncx;
  Int ncm = common.ncm;
  Int nd  = common.nd;
  Int npe = common.npe;
  Int ne  = common.ne;
  Int numPoints = npe * ne;
  dstype time = common.time;

  typename Model::MonitorFn monitor{};
  monitor(f, xg, udg, odg, wdg,
          app.uinf, app.physicsparam, time,
          common.modelnumber, numPoints,
          nc_sol, ncu, nd, ncx, nco, ncw, ncm, npe, ne);
}

template <typename Model = DefaultModel>
void AvfieldDriver(dstype* f, const dstype* xg, const dstype* udg,
                   const dstype* odg, const dstype* wdg,
                   meshstruct &mesh, masterstruct &master,
                   appstruct &app, solstruct &sol, tempstruct &temp,
                   commonstruct &common, Int backend)
{
  Int nc  = common.nc;
  Int ncu = common.ncu;
  Int ncw = common.ncw;
  Int nco = common.nco;
  Int ncx = common.ncx;
  Int nd  = common.nd;
  Int npe = common.npe;
  Int ne  = common.ne;
  Int numPoints = npe * ne;
  dstype time = common.time;

  typename Model::AvfieldFn avfield{};
  avfield(f, xg, udg, odg, wdg,
          app.uinf, app.physicsparam, time, common.modelnumber,
          numPoints, nc, ncu, nd, ncx, nco, ncw,
          nco, npe, ne);
}

template <typename Model = DefaultModel>
void EosDriver(dstype* f, const dstype* xg, const dstype* udg,
               const dstype* odg, const dstype* wdg,
               meshstruct &mesh, masterstruct &master,
               appstruct &app, solstruct &sol, tempstruct &temp,
               commonstruct &common, Int npe, Int e1, Int e2, Int backend)
{
  Int nc  = common.nc;
  Int ncu = common.ncu;
  Int ncw = common.ncw;
  Int nco = common.nco;
  Int ncx = common.ncx;
  Int nd  = common.nd;
  Int ne  = e2 - e1;
  Int numPoints = npe * ne;
  dstype time = common.time;

  typename Model::EoSFn eos{};
  eos(f, xg, udg, odg, wdg,
      app.uinf, app.physicsparam, time, common.modelnumber,
      numPoints, nc, ncu, nd, ncx, nco, ncw,
      ncw, npe, ne);
}

template <typename Model = DefaultModel>
void EosduDriver(dstype* f, const dstype* xg, const dstype* udg,
                 const dstype* odg, const dstype* wdg,
                 meshstruct &mesh, masterstruct &master,
                 appstruct &app, solstruct &sol, tempstruct &temp,
                 commonstruct &common, Int npe, Int e1, Int e2, Int backend)
{
  Int nc  = common.nc;
  Int ncu = common.ncu;
  Int ncw = common.ncw;
  Int nco = common.nco;
  Int ncx = common.ncx;
  Int nd  = common.nd;
  Int ne  = e2 - e1;
  Int numPoints = npe * ne;
  dstype time = common.time;

  typename Model::EoSduFn eosdu{};
  eosdu(f, xg, udg, odg, wdg,
        app.uinf, app.physicsparam, time, common.modelnumber,
        numPoints, nc, ncu, nd, ncx, nco, ncw,
        ncw * ncu, npe, ne);
}

template <typename Model = DefaultModel>
void EosdwDriver(dstype* f, const dstype* xg, const dstype* udg,
                 const dstype* odg, const dstype* wdg,
                 meshstruct &mesh, masterstruct &master,
                 appstruct &app, solstruct &sol, tempstruct &temp,
                 commonstruct &common, Int npe, Int e1, Int e2, Int backend)
{
  Int nc  = common.nc;
  Int ncu = common.ncu;
  Int ncw = common.ncw;
  Int nco = common.nco;
  Int ncx = common.ncx;
  Int nd  = common.nd;
  Int ne  = e2 - e1;
  Int numPoints = npe * ne;
  dstype time = common.time;

  typename Model::EoSdwFn eosdw{};
  eosdw(f, xg, udg, odg, wdg,
        app.uinf, app.physicsparam, time, common.modelnumber,
        numPoints, nc, ncu, nd, ncx, nco, ncw,
        ncw * ncw, npe, ne);
}

// ======================================================
// Time-dependent forcing
// ======================================================

template <typename Model = DefaultModel>
void TdfuncDriver(dstype* f, const dstype* xg, const dstype* udg,
                  const dstype* odg, const dstype* wdg,
                  meshstruct &mesh, masterstruct &master,
                  appstruct &app, solstruct &sol, tempstruct &temp,
                  commonstruct &common, Int nge, Int e1, Int e2, Int backend)
{
  Int nc  = common.nc;
  Int ncu = common.ncu;
  Int ncw = common.ncw;
  Int nco = common.nco;
  Int ncx = common.ncx;
  Int nd  = common.nd;
  Int numPoints = nge * (e2 - e1);
  dstype time = common.time;

  typename Model::TdfuncFn tdfunc{};
  tdfunc(f, xg, udg, odg, wdg,
         app.uinf, app.physicsparam, time,
         common.modelnumber, numPoints,
         nc, ncu, nd, ncx, nco, ncw);
}

// ======================================================
// Fhat / Fbou / Uhat / Ubou (non-HDG boundary/interface)
// ======================================================

template <typename Model = DefaultModel>
void FhatDriver(dstype* fg, const dstype* xg,
                const dstype* ug1, const dstype* ug2,
                const dstype* og1, const dstype* og2,
                const dstype* wg1, const dstype* wg2,
                const dstype* uh, const dstype* nl,
                meshstruct &mesh, masterstruct &master,
                appstruct &app, solstruct &sol, tempstruct &tmp,
                commonstruct &common, Int ngf, Int f1, Int f2, Int backend)
{
  Int nc  = common.nc;
  Int ncu = common.ncu;
  Int ncw = common.ncw;
  Int nco = common.nco;
  Int ncx = common.ncx;
  Int nd  = common.nd;
  Int numPoints = ngf * (f2 - f1);
  Int M = numPoints * ncu;
  Int N = numPoints * ncu * nd;
  Int ntau = common.ntau;
  dstype time = common.time;

  if (common.extFhat == 1) {
    typename Model::FhatFn fhat{};
    fhat(fg, xg, ug1, ug2, og1, og2, wg1, wg2, uh, nl,
         app.tau, app.uinf, app.physicsparam,
         time, common.modelnumber,
         numPoints, nc, ncu, nd, ncx, nco, ncw);
  } else {
    // left flux
    FluxDriver<Model>(fg, xg, ug1, og1, wg1,
                      mesh, master, app, sol, tmp, common,
                      ngf, f1, f2, backend);
    // right flux
    dstype* fg2 = &fg[N];
    FluxDriver<Model>(fg2, xg, ug2, og2, wg2,
                      mesh, master, app, sol, tmp, common,
                      ngf, f1, f2, backend);

    // Part 1: fh = fg dot nl
    AverageFlux(fg, N);
    AverageFluxDotNormal(fg, nl, N, M, numPoints, nd);

    // Part 2: Contribution due to tau*(U-UH)
    if (common.extStab >= 1) {
      typename Model::StabFn stab{};
      stab(fg, xg, ug1, ug2, og1, og2, wg1, wg2, uh, nl,
           app.tau, app.uinf, app.physicsparam,
           time, common.modelnumber,
           numPoints, nc, ncu, nd, ncx, nco, ncw);
    } else if (ntau == 0) {
      AddStabilization1(fg, ug1, ug2, app.tau, M);
    } else if (ntau == 1) { // constant scalar
      AddStabilization1(fg, ug1, ug2, app.tau, M);
    } else if (ntau == ncu) { // constant diagonal tensor
      AddStabilization2(fg, ug1, ug2, app.tau, M, numPoints);
    } else if (ntau == ncu * ncu) { // constant full tensor
      AddStabilization3(fg, ug1, ug2, app.tau, M, numPoints, ncu);
    } else {
      printf("Stabilization option is not implemented");
      exit(-1);
    }
  }
}

template <typename Model = DefaultModel>
void FbouDriver(dstype* fb, const dstype* xg,
                const dstype* udg, const dstype* odg,
                const dstype* wdg, const dstype* uhg,
                const dstype* nl,
                meshstruct &mesh, masterstruct &master,
                appstruct &app, solstruct &sol, tempstruct &temp,
                commonstruct &common, Int ngf, Int f1, Int f2,
                Int ib, Int backend)
{
  Int nc  = common.nc;
  Int ncu = common.ncu;
  Int ncw = common.ncw;
  Int nco = common.nco;
  Int ncx = common.ncx;
  Int nd  = common.nd;
  Int numPoints = ngf * (f2 - f1);
  dstype time = common.time;

  typename Model::FbouFn fbou{};
  fbou(fb, xg, udg, odg, wdg, uhg, nl,
       app.tau, app.uinf, app.physicsparam, time,
       common.modelnumber, ib, numPoints,
       nc, ncu, nd, ncx, nco, ncw);
}

template <typename Model = DefaultModel>
void UhatDriver(dstype* fg, dstype* xg,
                dstype* ug1, dstype* ug2,
                const dstype* og1, const dstype* og2,
                const dstype* wg1, const dstype* wg2,
                const dstype* uh, const dstype* nl,
                meshstruct &mesh, masterstruct &master,
                appstruct &app, solstruct &sol, tempstruct &tmp,
                commonstruct &common, Int ngf, Int f1, Int f2, Int backend)
{
  Int nc  = common.nc;
  Int ncu = common.ncu;
  Int ncw = common.ncw;
  Int nco = common.nco;
  Int ncx = common.ncx;
  Int nd  = common.nd;
  Int numPoints = ngf * (f2 - f1);
  dstype time = common.time;

  if (common.extUhat == 1) {
    typename Model::UhatFn uhat{};
    uhat(fg, xg, ug1, ug2, og1, og2, wg1, wg2, uh, nl,
         app.tau, app.uinf, app.physicsparam,
         time, common.modelnumber,
         numPoints, nc, ncu, nd, ncx, nco, ncw);
  } else {
    ArrayAXPBY(fg, ug1, ug2,
               (dstype)0.5, (dstype)0.5,
               ngf * common.ncu * (f2 - f1));
  }
}

template <typename Model = DefaultModel>
void UbouDriver(dstype* ub, const dstype* xg,
                const dstype* udg, const dstype* odg,
                const dstype* wdg, const dstype* uhg,
                const dstype* nl,
                meshstruct &mesh, masterstruct &master,
                appstruct &app, solstruct &sol, tempstruct &temp,
                commonstruct &common, Int ngf, Int f1, Int f2,
                Int ib, Int backend)
{
  Int nc  = common.nc;
  Int ncu = common.ncu;
  Int ncw = common.ncw;
  Int nco = common.nco;
  Int ncx = common.ncx;
  Int nd  = common.nd;
  Int numPoints = ngf * (f2 - f1);
  dstype time = common.time;

  typename Model::UbouFn ubou{};
  ubou(ub, xg, udg, odg, wdg, uhg, nl,
       app.tau, app.uinf, app.physicsparam, time,
       common.modelnumber, ib, numPoints,
       nc, ncu, nd, ncx, nco, ncw);
}

// ======================================================
// Init (Kokkos) drivers
// ======================================================

template <typename Model = DefaultModel>
void InitodgDriver(dstype* f, const dstype* xg,
                   appstruct &app, Int ncx, Int nco,
                   Int npe, Int ne, Int backend)
{
  Int numPoints   = npe * ne;
  Int modelnumber = app.flag[12];

  typename Model::InitodgFn initodg{};
  initodg(f, xg, app.uinf, app.physicsparam,
          modelnumber, numPoints, ncx, nco, npe, ne);
}

template <typename Model = DefaultModel>
void InitqDriver(dstype* f, const dstype* xg,
                 appstruct &app, Int ncx, Int nc,
                 Int npe, Int ne, Int backend)
{
  Int numPoints   = npe * ne;
  Int modelnumber = app.flag[12];

  typename Model::InitqFn initq{};
  initq(f, xg, app.uinf, app.physicsparam,
        modelnumber, numPoints, ncx, nc, npe, ne);
}

template <typename Model = DefaultModel>
void InitudgDriver(dstype* f, const dstype* xg,
                   appstruct &app, Int ncx, Int nc,
                   Int npe, Int ne, Int backend)
{
  Int numPoints   = npe * ne;
  Int modelnumber = app.flag[12];

  typename Model::InitudgFn initudg{};
  initudg(f, xg, app.uinf, app.physicsparam,
          modelnumber, numPoints, ncx, nc, npe, ne);
}

template <typename Model = DefaultModel>
void InituDriver(dstype* f, const dstype* xg,
                 appstruct &app, Int ncx, Int nc,
                 Int npe, Int ne, Int backend)
{
  Int numPoints   = npe * ne;
  Int modelnumber = app.flag[12];

  typename Model::InituFn initu{};
  initu(f, xg, app.uinf, app.physicsparam,
        modelnumber, numPoints, ncx, nc, npe, ne);
}

template <typename Model = DefaultModel>
void InitwdgDriver(dstype* f, const dstype* xg,
                   appstruct &app, Int ncx, Int ncw,
                   Int npe, Int ne, Int backend)
{
  Int numPoints   = npe * ne;
  Int modelnumber = app.flag[12];

  typename Model::InitwdgFn initwdg{};
  initwdg(f, xg, app.uinf, app.physicsparam,
          modelnumber, numPoints, ncx, ncw, npe, ne);
}

// ======================================================
// Init (CPU) drivers
// ======================================================

template <typename Model = DefaultModel>
void cpuInitodgDriver(dstype* f, const dstype* xg,
                      appstruct &app, Int ncx, Int nco,
                      Int npe, Int ne, Int backend)
{
  Int numPoints   = npe * ne;
  Int modelnumber = app.flag[12];

  typename Model::CpuInitodgFn initodg{};
  initodg(f, xg, app.uinf, app.physicsparam,
          modelnumber, numPoints, ncx, nco, npe, ne);
}

template <typename Model = DefaultModel>
void cpuInitqDriver(dstype* f, const dstype* xg,
                    appstruct &app, Int ncx, Int nc,
                    Int npe, Int ne, Int backend)
{
  Int numPoints   = npe * ne;
  Int modelnumber = app.flag[12];

  typename Model::CpuInitqFn initq{};
  initq(f, xg, app.uinf, app.physicsparam,
        modelnumber, numPoints, ncx, nc, npe, ne);
}

template <typename Model = DefaultModel>
void cpuInitudgDriver(dstype* f, const dstype* xg,
                      appstruct &app, Int ncx, Int nc,
                      Int npe, Int ne, Int backend)
{
  Int numPoints   = npe * ne;
  Int modelnumber = app.flag[12];

  typename Model::CpuInitudgFn initudg{};
  initudg(f, xg, app.uinf, app.physicsparam,
          modelnumber, numPoints, ncx, nc, npe, ne);
}

template <typename Model = DefaultModel>
void cpuInituDriver(dstype* f, const dstype* xg,
                    appstruct &app, Int ncx, Int nc,
                    Int npe, Int ne, Int backend)
{
  Int numPoints   = npe * ne;
  Int modelnumber = app.flag[12];

  typename Model::CpuInituFn initu{};
  initu(f, xg, app.uinf, app.physicsparam,
        modelnumber, numPoints, ncx, nc, npe, ne);
}

template <typename Model = DefaultModel>
void cpuInitwdgDriver(dstype* f, const dstype* xg,
                      appstruct &app, Int ncx, Int ncw,
                      Int npe, Int ne, Int backend)
{
  Int numPoints   = npe * ne;
  Int modelnumber = app.flag[12];

  typename Model::CpuInitwdgFn initwdg{};
  initwdg(f, xg, app.uinf, app.physicsparam,
          modelnumber, numPoints, ncx, ncw, npe, ne);
}

// ======================================================
// HDG volume flux (already templated, just uses HdgFluxFn)
// ======================================================

template <typename Model = DefaultModel>
void FluxDriver(dstype* f, dstype* f_udg, dstype* f_wdg,
                const dstype* xg, dstype* udg,
                const dstype* odg, const dstype* wdg,
                meshstruct &mesh, masterstruct &master,
                appstruct &app, solstruct &sol, tempstruct &temp,
                commonstruct &common, Int nge, Int e1, Int e2, Int backend)
{
  Int nc  = common.nc;
  Int ncu = common.ncu;
  Int ncw = common.ncw;
  Int nco = common.nco;
  Int ncx = common.ncx;
  Int nd  = common.nd;
  Int numPoints = nge * (e2 - e1);
  dstype time = common.time;

  typename Model::HdgFluxFn flux{};
  flux(f, f_udg, f_wdg, xg, udg, odg, wdg,
       app.uinf, app.physicsparam, time,
       common.modelnumber, numPoints,
       nc, ncu, nd, ncx, nco, ncw);
}

// ======================================================
// HDG volume/source/EoS
// ======================================================

template <typename Model = DefaultModel>
void SourceDriver(dstype* f, dstype* f_udg, dstype* f_wdg,
                  const dstype* xg, const dstype* udg,
                  const dstype* odg, const dstype* wdg,
                  meshstruct &mesh, masterstruct &master,
                  appstruct &app, solstruct &sol, tempstruct &temp,
                  commonstruct &common, Int nge, Int e1, Int e2, Int backend)
{
  Int nc  = common.nc;
  Int ncu = common.ncu;
  Int ncw = common.ncw;
  Int nco = common.nco;
  Int ncx = common.ncx;
  Int nd  = common.nd;
  Int numPoints = nge * (e2 - e1);
  dstype time = common.time;

  typename Model::HdgSourceFn src{};
  src(f, f_udg, f_wdg, xg, udg, odg, wdg,
      app.uinf, app.physicsparam, time,
      common.modelnumber, numPoints,
      nc, ncu, nd, ncx, nco, ncw);
}

template <typename Model = DefaultModel>
void SourcewDriver(dstype* f, dstype* f_udg, dstype* f_wdg,
                   const dstype* xg, const dstype* udg,
                   const dstype* odg, const dstype* wdg,
                   meshstruct &mesh, masterstruct &master,
                   appstruct &app, solstruct &sol, tempstruct &temp,
                   commonstruct &common, Int nge, Int e1, Int e2, Int backend)
{
  Int nc  = common.nc;
  Int ncu = common.ncu;
  Int ncw = common.ncw;
  Int nco = common.nco;
  Int ncx = common.ncx;
  Int nd  = common.nd;
  Int numPoints = nge * (e2 - e1);
  dstype time = common.time;

  typename Model::HdgSourcewFn srcw{};
  srcw(f, f_udg, f_wdg, xg, udg, odg, wdg,
       app.uinf, app.physicsparam, time,
       common.modelnumber, numPoints,
       nc, ncu, nd, ncx, nco, ncw);
}

template <typename Model = DefaultModel>
void SourcewDriver(dstype* f, dstype* f_wdg,
                   const dstype* xg, const dstype* udg,
                   const dstype* odg, const dstype* wdg,
                   meshstruct &mesh, masterstruct &master,
                   appstruct &app, solstruct &sol, tempstruct &temp,
                   commonstruct &common, Int nge, Int e1, Int e2, Int backend)
{
  Int nc  = common.nc;
  Int ncu = common.ncu;
  Int ncw = common.ncw;
  Int nco = common.nco;
  Int ncx = common.ncx;
  Int nd  = common.nd;
  Int numPoints = nge * (e2 - e1);
  dstype time = common.time;

  typename Model::HdgSourcewonlyFn srcwonly{};
  srcwonly(f, f_wdg, xg, udg, odg, wdg,
           app.uinf, app.physicsparam, time,
           common.modelnumber, numPoints,
           nc, ncu, nd, ncx, nco, ncw);
}

template <typename Model = DefaultModel>
void EosDriver(dstype* f, dstype* f_udg, dstype* f_wdg,
               const dstype* xg, const dstype* udg,
               const dstype* odg, const dstype* wdg,
               meshstruct &mesh, masterstruct &master,
               appstruct &app, solstruct &sol, tempstruct &temp,
               commonstruct &common, Int nge, Int e1, Int e2, Int backend)
{
  Int nc  = common.nc;
  Int ncu = common.ncu;
  Int ncw = common.ncw;
  Int nco = common.nco;
  Int ncx = common.ncx;
  Int nd  = common.nd;
  Int numPoints = nge * (e2 - e1);
  dstype time = common.time;

  typename Model::HdgEoSFn eos{};
  eos(f, f_udg, f_wdg, xg, udg, odg, wdg,
      app.uinf, app.physicsparam, time,
      common.modelnumber, numPoints,
      nc, ncu, nd, ncx, nco, ncw);
}

// ======================================================
// HDG Fbou / Fint (with and without derivatives)
// ======================================================

template <typename Model = DefaultModel>
void FbouDriver(dstype* f,  dstype* f_udg, dstype* f_wdg, dstype* f_uhg,
                dstype* xg, const dstype* udg,
                const dstype* odg, const dstype* wdg,
                dstype* uhg, const dstype* nl,
                meshstruct &mesh, masterstruct &master,
                appstruct &app, solstruct &sol, tempstruct &temp,
                commonstruct &common, Int nga, Int ib, Int backend)
{
  Int nc  = common.nc;
  Int ncu = common.ncu;
  Int ncw = common.ncw;
  Int nco = common.nco;
  Int ncx = common.ncx;
  Int nd  = common.nd;
  Int numPoints = nga;
  dstype time = common.time;

  typename Model::HdgFbouFn fbou{};
  fbou(f, f_udg, f_wdg, f_uhg,
       xg, udg, odg, wdg,
       uhg, nl, app.tau, app.uinf, app.physicsparam,
       time, common.modelnumber, ib,
       numPoints, nc, ncu, nd, ncx, nco, ncw);
}

template <typename Model = DefaultModel>
void FbouDriver(dstype* f, dstype* xg,
                const dstype* udg, const dstype* odg,
                const dstype* wdg, dstype* uhg,
                const dstype* nl,
                meshstruct &mesh, masterstruct &master,
                appstruct &app, solstruct &sol, tempstruct &temp,
                commonstruct &common, Int nga, Int ib, Int backend)
{
  Int nc  = common.nc;
  Int ncu = common.ncu;
  Int ncw = common.ncw;
  Int nco = common.nco;
  Int ncx = common.ncx;
  Int nd  = common.nd;
  Int numPoints = nga;
  dstype time = common.time;

  typename Model::HdgFbouonlyFn fbouonly{};
  fbouonly(f, xg, udg, odg, wdg,
           uhg, nl, app.tau, app.uinf, app.physicsparam,
           time, common.modelnumber, ib,
           numPoints, nc, ncu, nd, ncx, nco, ncw);
}

template <typename Model = DefaultModel>
void FintDriver(dstype* f,  dstype* f_udg, dstype* f_wdg, dstype* f_uhg,
                dstype* xg, const dstype* udg,
                const dstype* odg, const dstype* wdg,
                dstype* uhg, const dstype* nl,
                meshstruct &mesh, masterstruct &master,
                appstruct &app, solstruct &sol, tempstruct &temp,
                commonstruct &common, Int nga, Int ib, Int backend)
{
  Int nc  = common.nc;
  Int ncu = common.ncu;
  Int ncw = common.ncw;
  Int nco = common.nco;
  Int ncx = common.ncx;
  Int nd  = common.nd;
  Int numPoints = nga;
  dstype time = common.time;

  typename Model::HdgFintFn fint{};
  fint(f, f_udg, f_wdg, f_uhg,
       xg, udg, odg, wdg,
       uhg, nl, app.tau, app.uinf, app.physicsparam,
       time, common.modelnumber, ib,
       numPoints, nc, ncu, nd, ncx, nco, ncw);
}

template <typename Model = DefaultModel>
void FintDriver(dstype* f, dstype* xg,
                const dstype* udg, const dstype* odg,
                const dstype* wdg, dstype* uhg,
                const dstype* nl,
                meshstruct &mesh, masterstruct &master,
                appstruct &app, solstruct &sol, tempstruct &temp,
                commonstruct &common, Int nga, Int ib, Int backend)
{
  Int nc  = common.nc;
  Int ncu = common.ncu;
  Int ncw = common.ncw;
  Int nco = common.nco;
  Int ncx = common.ncx;
  Int nd  = common.nd;
  Int numPoints = nga;
  dstype time = common.time;

  typename Model::HdgFintonlyFn fintonly{};
  fintonly(f, xg, udg, odg, wdg,
           uhg, nl, app.tau, app.uinf, app.physicsparam,
           time, common.modelnumber, ib,
           numPoints, nc, ncu, nd, ncx, nco, ncw);
}

// ======================================================
// HDG Fhat (derivatives & non-derivative variant)
// ======================================================

template <typename Model = DefaultModel>
void FhatDriver(dstype* f,  dstype* f_udg, dstype* f_wdg, dstype* f_uhg,
                const dstype* xg, dstype* udg,
                const dstype* odg, const dstype* wdg,
                const dstype* uhg, dstype* nl,
                meshstruct &mesh, masterstruct &master,
                appstruct &app, solstruct &sol, tempstruct &temp,
                commonstruct &common, Int nga, Int backend)
{
  Int nc  = common.nc;
  Int ncu = common.ncu;
  Int ncw = common.ncw;
  Int nco = common.nco;
  Int ncx = common.ncx;
  Int nd  = common.nd;
  Int numPoints = nga;
  Int M = numPoints * ncu;
  Int N = numPoints * ncu * nd;
  dstype time = common.time;

  // copy u-compoment of udg to f_uhg
  ArrayCopy(f_uhg, udg, numPoints * ncu);

  // copy uhg to udg, so udg = (uh, q)
  ArrayCopy(udg, uhg, numPoints * ncu);

  // HdgFlux: f, f_udg, f_wdg
  typename Model::HdgFluxFn hdgflux{};
  hdgflux(f, f_udg, f_wdg, xg, udg, odg, wdg,
          app.uinf, app.physicsparam, time,
          common.modelnumber, numPoints,
          nc, ncu, nd, ncx, nco, ncw);

  // copy f_uhg back into u-component of udg
  ArrayCopy(udg, f_uhg, numPoints * ncu);

  // f = ng * ncu
  FluxDotNormal(f, f, nl, M, numPoints, nd);

  // f_udg = ng * ncu * nc
  for (Int n = 0; n < nc; n++) {
    FluxDotNormal(&f_udg[M * n], &f_udg[N * n],
                  nl, M, numPoints, nd);
  }

  // f_wdg = ng * ncu * ncw
  if ((ncw > 0) & (common.wave == 0)) {
    for (Int n = 0; n < ncw; n++) {
      FluxDotNormal(&f_wdg[M * n], &f_wdg[N * n],
                    nl, M, numPoints, nd);
    }
  }

  // f_uhg = ng * ncu * ncu
  ArrayCopy(f_uhg, f_udg, numPoints * ncu * ncu); // copy f_u to f_uhg

  // reset f_u to zero
  ArraySetValue(f_udg, zero, numPoints * ncu * ncu);

  // add tau*(u - uh) to f
  AddStabilization1(f, f_udg, f_uhg, udg, uhg, app.tau, M, numPoints);
}

template <typename Model = DefaultModel>
void FhatDriver(dstype* f, dstype* u,
                const dstype* xg, dstype* udg,
                const dstype* odg, const dstype* wdg,
                const dstype* uhg, dstype* nl,
                meshstruct &mesh, masterstruct &master,
                appstruct &app, solstruct &sol, tempstruct &temp,
                commonstruct &common, Int nga, Int backend)
{
  Int nc  = common.nc;
  Int ncu = common.ncu;
  Int ncw = common.ncw;
  Int nco = common.nco;
  Int ncx = common.ncx;
  Int nd  = common.nd;
  Int numPoints = nga;
  Int M = numPoints * ncu;
  dstype time = common.time;

  // copy u-compoment of udg to u
  ArrayCopy(u, udg, numPoints * ncu);

  // copy uhg to udg, so udg = (uh, q)
  ArrayCopy(udg, uhg, numPoints * ncu);

  // f = ng * ncu * nd
  typename Model::FluxFn flux{};
  flux(f, xg, udg, odg, wdg,
       app.uinf, app.physicsparam, time,
       common.modelnumber, numPoints,
       nc, ncu, nd, ncx, nco, ncw);

  // restore u into udg
  ArrayCopy(udg, u, numPoints * ncu);

  // f = ng * ncu
  FluxDotNormal(f, f, nl, M, numPoints, nd);

  // add tau*(u - uh) to f
  AddStabilization1(f, udg, uhg, app.tau, M);
}

// ======================================================
// Visualization and QoI
// ======================================================

template <typename Model = DefaultModel>
void VisScalarsDriver(dstype* f, const dstype* xg,
                      const dstype* udg, const dstype* odg,
                      const dstype* wdg,
                      meshstruct &mesh, masterstruct &master,
                      appstruct &app, solstruct &sol, tempstruct &temp,
                      commonstruct &common, Int nge, Int e1, Int e2,
                      Int backend)
{
  Int nc  = common.nc;
  Int ncu = common.ncu;
  Int ncw = common.ncw;
  Int nco = common.nco;
  Int ncx = common.ncx;
  Int nd  = common.nd;
  Int numPoints = nge * (e2 - e1);
  dstype time = common.time;

  typename Model::VisScalarsFn vis{};
  vis(f, xg, udg, odg, wdg,
      app.uinf, app.physicsparam, time,
      common.modelnumber, numPoints,
      nc, ncu, nd, ncx, nco, ncw);
}

template <typename Model = DefaultModel>
void VisVectorsDriver(dstype* f, const dstype* xg,
                      const dstype* udg, const dstype* odg,
                      const dstype* wdg,
                      meshstruct &mesh, masterstruct &master,
                      appstruct &app, solstruct &sol, tempstruct &temp,
                      commonstruct &common, Int nge, Int e1, Int e2,
                      Int backend)
{
  Int nc  = common.nc;
  Int ncu = common.ncu;
  Int ncw = common.ncw;
  Int nco = common.nco;
  Int ncx = common.ncx;
  Int nd  = common.nd;
  Int numPoints = nge * (e2 - e1);
  dstype time = common.time;

  typename Model::VisVectorsFn vis{};
  vis(f, xg, udg, odg, wdg,
      app.uinf, app.physicsparam, time,
      common.modelnumber, numPoints,
      nc, ncu, nd, ncx, nco, ncw);
}

template <typename Model = DefaultModel>
void VisTensorsDriver(dstype* f, const dstype* xg,
                      const dstype* udg, const dstype* odg,
                      const dstype* wdg,
                      meshstruct &mesh, masterstruct &master,
                      appstruct &app, solstruct &sol, tempstruct &temp,
                      commonstruct &common, Int nge, Int e1, Int e2,
                      Int backend)
{
  Int nc  = common.nc;
  Int ncu = common.ncu;
  Int ncw = common.ncw;
  Int nco = common.nco;
  Int ncx = common.ncx;
  Int nd  = common.nd;
  Int numPoints = nge * (e2 - e1);
  dstype time = common.time;

  typename Model::VisTensorsFn vis{};
  vis(f, xg, udg, odg, wdg,
      app.uinf, app.physicsparam, time,
      common.modelnumber, numPoints,
      nc, ncu, nd, ncx, nco, ncw);
}

template <typename Model = DefaultModel>
void QoIvolumeDriver(dstype* f, const dstype* xg,
                     const dstype* udg, const dstype* odg,
                     const dstype* wdg,
                     meshstruct &mesh, masterstruct &master,
                     appstruct &app, solstruct &sol, tempstruct &temp,
                     commonstruct &common, Int nge, Int e1, Int e2,
                     Int backend)
{
  Int nc  = common.nc;
  Int ncu = common.ncu;
  Int ncw = common.ncw;
  Int nco = common.nco;
  Int ncx = common.ncx;
  Int nd  = common.nd;
  Int numPoints = nge * (e2 - e1);
  dstype time = common.time;

  typename Model::QoIvolumeFn qvol{};
  qvol(f, xg, udg, odg, wdg,
       app.uinf, app.physicsparam, time,
       common.modelnumber, numPoints,
       nc, ncu, nd, ncx, nco, ncw);
}

template <typename Model = DefaultModel>
void QoIboundaryDriver(dstype* fb, const dstype* xg,
                       const dstype* udg, const dstype* odg,
                       const dstype* wdg, const dstype* uhg,
                       const dstype* nl,
                       meshstruct &mesh, masterstruct &master,
                       appstruct &app, solstruct &sol, tempstruct &temp,
                       commonstruct &common, Int ngf, Int f1, Int f2,
                       Int ib, Int backend)
{
  Int nc  = common.nc;
  Int ncu = common.ncu;
  Int ncw = common.ncw;
  Int nco = common.nco;
  Int ncx = common.ncx;
  Int nd  = common.nd;
  Int numPoints = ngf * (f2 - f1);
  dstype time = common.time;

  typename Model::QoIboundaryFn qbd{};
  qbd(fb, xg, udg, odg, wdg, uhg, nl,
      app.tau, app.uinf, app.physicsparam, time,
      common.modelnumber, ib, numPoints,
      nc, ncu, nd, ncx, nco, ncw);
}
