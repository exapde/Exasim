/**
 * @file model_drivers_abi.cpp
 * @brief Core-side driver wrappers built on top of ExasimDriverabi.
 *
 * These wrappers preserve the legacy Exasim driver argument mapping while
 * dispatching the low-level kernel calls through a selected provider abi.
 * They intentionally keep the existing LDG/HDG data flow unchanged.
 */

#ifndef __EXASIM_MODEL_DRIVERS_ABI
#define __EXASIM_MODEL_DRIVERS_ABI

#include "driver_abi.h"

void FluxDriver(dstype* f, const dstype* xg, const dstype* udg,
                const dstype* odg, const dstype* wdg,
                ExasimDriverABI& abi, meshstruct& mesh, masterstruct& master, appstruct& app,
                solstruct& sol, tempstruct& temp, commonstruct& common,
                Int nge, Int e1, Int e2, Int backend)
{
    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nd = common.nd;
    Int numPoints = nge * (e2 - e1);
    dstype time = common.time;

    abi.KokkosFlux(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time,
                   common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);
}

void SourceDriver(dstype* f, const dstype* xg, const dstype* udg,
                  const dstype* odg, const dstype* wdg,
                  ExasimDriverABI& abi, meshstruct& mesh, masterstruct& master, appstruct& app,
                  solstruct& sol, tempstruct& temp, commonstruct& common,
                  Int nge, Int e1, Int e2, Int backend)
{
    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nd = common.nd;
    Int numPoints = nge * (e2 - e1);
    dstype time = common.time;

    abi.KokkosSource(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time,
                     common.modelnumber, numPoints, nc, ncu, nd, ncx, nco,
                     ncw);
}

void SourcewDriver(dstype* f, const dstype* xg, const dstype* udg,
                   const dstype* odg, const dstype* wdg,
                   ExasimDriverABI& abi, meshstruct& mesh, masterstruct& master, appstruct& app,
                   solstruct& sol, tempstruct& temp, commonstruct& common,
                   Int npe, Int e1, Int e2, Int backend)
{
    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nd = common.nd;
    Int ne = e2 - e1;
    Int numPoints = npe * ne;
    dstype time = common.time;

    abi.KokkosSourcew(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time,
                      common.modelnumber, numPoints, nc, ncu, nd, ncx, nco,
                      ncw, ncw, npe, ne);
}

void OutputDriver(dstype* f, const dstype* xg, const dstype* udg,
                  const dstype* odg, const dstype* wdg,
                  ExasimDriverABI& abi, meshstruct& mesh, masterstruct& master, appstruct& app,
                  solstruct& sol, tempstruct& temp, commonstruct& common,
                  Int backend)
{
    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nce = common.nce;
    Int nd = common.nd;
    Int npe = common.npe;
    Int ne = common.ne;
    Int numPoints = npe * ne;
    dstype time = common.time;

    abi.KokkosOutput(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time,
                     common.modelnumber, numPoints, nc, ncu, nd, ncx, nco,
                     ncw, nce, npe, ne);
}

void MonitorDriver(dstype* f, Int nc_sol, const dstype* xg,
                   const dstype* udg, const dstype* odg, const dstype* wdg,
                   ExasimDriverABI& abi, meshstruct& mesh,
                   masterstruct& master, appstruct& app, solstruct& sol,
                   tempstruct& temp, commonstruct& common, Int backend)
{
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int ncm = common.ncm;
    Int nd = common.nd;
    Int npe = common.npe;
    Int ne = common.ne;
    Int numPoints = npe * ne;
    dstype time = common.time;

    abi.KokkosMonitor(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time,
                      common.modelnumber, numPoints, nc_sol, ncu, nd, ncx,
                      nco, ncw, ncm, npe, ne);
}

void AvfieldDriver(dstype* f, const dstype* xg, const dstype* udg,
                   const dstype* odg, const dstype* wdg,
                   ExasimDriverABI& abi, meshstruct& mesh, masterstruct& master, appstruct& app,
                   solstruct& sol, tempstruct& temp, commonstruct& common,
                   Int backend)
{
    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nd = common.nd;
    Int npe = common.npe;
    Int ne = common.ne;
    Int numPoints = npe * ne;
    dstype time = common.time;

    abi.KokkosAvfield(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time,
                      common.modelnumber, numPoints, nc, ncu, nd, ncx, nco,
                      ncw, nco, npe, ne);
}

void EosDriver(dstype* f, const dstype* xg, const dstype* udg,
               const dstype* odg, const dstype* wdg,
               ExasimDriverABI& abi, meshstruct& mesh, masterstruct& master, appstruct& app,
               solstruct& sol, tempstruct& temp, commonstruct& common,
               Int npe, Int e1, Int e2, Int backend)
{
    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nd = common.nd;
    Int ne = e2 - e1;
    Int numPoints = npe * ne;
    dstype time = common.time;

    abi.KokkosEoS(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time,
                  common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw,
                  ncw, npe, ne);
}

void EosduDriver(dstype* f, const dstype* xg, const dstype* udg,
                 const dstype* odg, const dstype* wdg,
                 ExasimDriverABI& abi, meshstruct& mesh, masterstruct& master, appstruct& app,
                 solstruct& sol, tempstruct& temp, commonstruct& common,
                 Int npe, Int e1, Int e2, Int backend)
{
    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nd = common.nd;
    Int ne = e2 - e1;
    Int numPoints = npe * ne;
    dstype time = common.time;

    abi.KokkosEoSdu(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time,
                    common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw,
                    ncw * ncu, npe, ne);
}

void EosdwDriver(dstype* f, const dstype* xg, const dstype* udg,
                 const dstype* odg, const dstype* wdg,
                 ExasimDriverABI& abi, meshstruct& mesh, masterstruct& master, appstruct& app,
                 solstruct& sol, tempstruct& temp, commonstruct& common,
                 Int npe, Int e1, Int e2, Int backend)
{
    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nd = common.nd;
    Int ne = e2 - e1;
    Int numPoints = npe * ne;
    dstype time = common.time;

    abi.KokkosEoSdw(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time,
                    common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw,
                    ncw * ncw, npe, ne);
}

void TdfuncDriver(dstype* f, const dstype* xg, const dstype* udg,
                  const dstype* odg, const dstype* wdg,
                  ExasimDriverABI& abi, meshstruct& mesh, masterstruct& master, appstruct& app,
                  solstruct& sol, tempstruct& temp, commonstruct& common,
                  Int nge, Int e1, Int e2, Int backend)
{
    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nd = common.nd;
    Int numPoints = nge * (e2 - e1);
    dstype time = common.time;

    abi.KokkosTdfunc(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time,
                     common.modelnumber, numPoints, nc, ncu, nd, ncx, nco,
                     ncw);
}

void FhatDriver(dstype* fg, const dstype* xg, const dstype* ug1,
                const dstype* ug2, const dstype* og1, const dstype* og2,
                const dstype* wg1, const dstype* wg2, const dstype* uh,
                const dstype* nl, ExasimDriverABI& abi, meshstruct& mesh,
                masterstruct& master, appstruct& app, solstruct& sol,
                tempstruct& tmp, commonstruct& common, Int ngf, Int f1, Int f2,
                Int backend)
{
    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nd = common.nd;
    Int numPoints = ngf * (f2 - f1);
    Int M = numPoints * ncu;
    Int N = numPoints * ncu * nd;
    Int ntau = common.ntau;
    dstype time = common.time;

    if (common.extFhat == 1) {
        abi.KokkosFhat(fg, xg, ug1, ug2, og1, og2, wg1, wg2, uh, nl, app.tau,
                       app.uinf, app.physicsparam, time, common.modelnumber,
                       numPoints, nc, ncu, nd, ncx, nco, ncw);
    }
    else {
        FluxDriver(fg, xg, ug1, og1, wg1, abi, mesh, master, app, sol, tmp,
                   common, ngf, f1, f2, backend);
        dstype* fg2 = &fg[N];
        FluxDriver(fg2, xg, ug2, og2, wg2, abi, mesh, master, app, sol, tmp,
                   common, ngf, f1, f2, backend);

        AverageFlux(fg, N);
        AverageFluxDotNormal(fg, nl, N, M, numPoints, nd);

        if (common.extStab >= 1) {
            abi.KokkosStab(fg, xg, ug1, ug2, og1, og2, wg1, wg2, uh, nl,
                           app.tau, app.uinf, app.physicsparam, time,
                           common.modelnumber, numPoints, nc, ncu, nd, ncx,
                           nco, ncw);
        }
        else if (ntau == 0) {
            AddStabilization1(fg, ug1, ug2, app.tau, M);
        }
        else if (ntau == 1) {
            AddStabilization1(fg, ug1, ug2, app.tau, M);
        }
        else if (ntau == ncu) {
            AddStabilization2(fg, ug1, ug2, app.tau, M, numPoints);
        }
        else if (ntau == ncu * ncu) {
            AddStabilization3(fg, ug1, ug2, app.tau, M, numPoints, ncu);
        }
        else {
            printf("Stabilization option is not implemented");
            exit(-1);
        }
    }
}

void FbouDriver(dstype* fb, const dstype* xg, const dstype* udg,
                const dstype* odg, const dstype* wdg, const dstype* uhg,
                const dstype* nl, ExasimDriverABI& abi, meshstruct& mesh,
                masterstruct& master, appstruct& app, solstruct& sol,
                tempstruct& temp, commonstruct& common, Int ngf, Int f1,
                Int f2, Int ib, Int backend)
{
    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nd = common.nd;
    Int numPoints = ngf * (f2 - f1);
    dstype time = common.time;

    abi.KokkosFbou(fb, xg, udg, odg, wdg, uhg, nl, app.tau, app.uinf,
                   app.physicsparam, time, common.modelnumber, ib, numPoints,
                   nc, ncu, nd, ncx, nco, ncw);
}

void UhatDriver(dstype* fg, dstype* xg, dstype* ug1, dstype* ug2,
                const dstype* og1, const dstype* og2, const dstype* wg1,
                const dstype* wg2, const dstype* uh, const dstype* nl,
                ExasimDriverABI& abi, meshstruct& mesh, masterstruct& master,
                appstruct& app, solstruct& sol, tempstruct& tmp,
                commonstruct& common, Int ngf, Int f1, Int f2, Int backend)
{
    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nd = common.nd;
    Int numPoints = ngf * (f2 - f1);
    dstype time = common.time;

    if (common.extUhat == 1) {
        abi.KokkosUhat(fg, xg, ug1, ug2, og1, og2, wg1, wg2, uh, nl, app.tau,
                       app.uinf, app.physicsparam, time, common.modelnumber,
                       numPoints, nc, ncu, nd, ncx, nco, ncw);
    }
    else {
        ArrayAXPBY(fg, ug1, ug2, (dstype)0.5, (dstype)0.5,
                   ngf * common.ncu * (f2 - f1));
    }
}

void UbouDriver(dstype* ub, const dstype* xg, const dstype* udg,
                const dstype* odg, const dstype* wdg, const dstype* uhg,
                const dstype* nl, ExasimDriverABI& abi, meshstruct& mesh,
                masterstruct& master, appstruct& app, solstruct& sol,
                tempstruct& temp, commonstruct& common, Int ngf, Int f1,
                Int f2, Int ib, Int backend)
{
    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nd = common.nd;
    Int numPoints = ngf * (f2 - f1);
    dstype time = common.time;

    abi.KokkosUbou(ub, xg, udg, odg, wdg, uhg, nl, app.tau, app.uinf,
                   app.physicsparam, time, common.modelnumber, ib, numPoints,
                   nc, ncu, nd, ncx, nco, ncw);
}

void InitodgDriver(dstype* f, const dstype* xg,
                   ExasimDriverABI& abi, appstruct& app, Int ncx, Int nco, Int npe, Int ne,
                   Int backend)
{
    Int numPoints = npe * ne;
    Int modelnumber = app.modelnumber;

    abi.KokkosInitodg(f, xg, app.uinf, app.physicsparam, modelnumber,
                      numPoints, ncx, nco, npe, ne);
}

void InitqDriver(dstype* f, const dstype* xg,
                 ExasimDriverABI& abi, appstruct& app, Int ncx, Int nc, Int npe, Int ne,
                 Int backend)
{
    Int numPoints = npe * ne;
    Int modelnumber = app.modelnumber;

    abi.KokkosInitq(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints,
                    ncx, nc, npe, ne);
}

void InitudgDriver(dstype* f, const dstype* xg,
                   ExasimDriverABI& abi, appstruct& app, Int ncx, Int nc, Int npe, Int ne,
                   Int backend)
{
    Int numPoints = npe * ne;
    Int modelnumber = app.modelnumber;

    abi.KokkosInitudg(f, xg, app.uinf, app.physicsparam, modelnumber,
                      numPoints, ncx, nc, npe, ne);
}

void InituDriver(dstype* f, const dstype* xg,
                 ExasimDriverABI& abi, appstruct& app, Int ncx, Int nc, Int npe, Int ne,
                 Int backend)
{
    Int numPoints = npe * ne;
    Int modelnumber = app.modelnumber;

    abi.KokkosInitu(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints,
                    ncx, nc, npe, ne);
}

void InitwdgDriver(dstype* f, const dstype* xg,
                   ExasimDriverABI& abi, appstruct& app, Int ncx, Int ncw, Int npe, Int ne,
                   Int backend)
{
    Int numPoints = npe * ne;
    Int modelnumber = app.modelnumber;

    abi.KokkosInitwdg(f, xg, app.uinf, app.physicsparam, modelnumber,
                      numPoints, ncx, ncw, npe, ne);
}

void cpuInitodgDriver(dstype* f, const dstype* xg,
                      ExasimDriverABI& abi, appstruct& app, Int ncx, Int nco, Int npe, Int ne,
                      Int backend)
{
    Int numPoints = npe * ne;
    Int modelnumber = app.modelnumber;

    abi.cpuInitodg(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints,
                   ncx, nco, npe, ne);
}

void cpuInitqDriver(dstype* f, const dstype* xg,
                    ExasimDriverABI& abi, appstruct& app, Int ncx, Int nc, Int npe, Int ne,
                    Int backend)
{
    Int numPoints = npe * ne;
    Int modelnumber = app.modelnumber;

    abi.cpuInitq(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints,
                 ncx, nc, npe, ne);
}

void cpuInitudgDriver(dstype* f, const dstype* xg,
                      ExasimDriverABI& abi, appstruct& app, Int ncx, Int nc, Int npe, Int ne,
                      Int backend)
{
    Int numPoints = npe * ne;
    Int modelnumber = app.modelnumber;

    abi.cpuInitudg(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints,
                   ncx, nc, npe, ne);
}

void cpuInituDriver(dstype* f, const dstype* xg,
                    ExasimDriverABI& abi, appstruct& app, Int ncx, Int nc, Int npe, Int ne,
                    Int backend)
{
    Int numPoints = npe * ne;
    Int modelnumber = app.modelnumber;

    abi.cpuInitu(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints,
                 ncx, nc, npe, ne);
}

void cpuInitwdgDriver(dstype* f, const dstype* xg,
                      ExasimDriverABI& abi, appstruct& app, Int ncx, Int ncw, Int npe, Int ne,
                      Int backend)
{
    Int numPoints = npe * ne;
    Int modelnumber = app.modelnumber;

    abi.cpuInitwdg(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints,
                   ncx, ncw, npe, ne);
}

void FluxDriver(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xg,
                dstype* udg, const dstype* odg, const dstype* wdg,
                ExasimDriverABI& abi, meshstruct& mesh, masterstruct& master,
                appstruct& app, solstruct& sol, tempstruct& temp,
                commonstruct& common, Int nge, Int e1, Int e2, Int backend)
{
    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nd = common.nd;
    Int numPoints = nge * (e2 - e1);
    dstype time = common.time;

    abi.HdgFlux(f, f_udg, f_wdg, xg, udg, odg, wdg, app.uinf,
                app.physicsparam, time, common.modelnumber, numPoints, nc,
                ncu, nd, ncx, nco, ncw);
}

void SourceDriver(dstype* f, dstype* f_udg, dstype* f_wdg,
                  const dstype* xg, const dstype* udg, const dstype* odg,
                  const dstype* wdg, ExasimDriverABI& abi, meshstruct& mesh,
                  masterstruct& master, appstruct& app, solstruct& sol,
                  tempstruct& temp, commonstruct& common, Int nge, Int e1,
                  Int e2, Int backend)
{
    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nd = common.nd;
    Int numPoints = nge * (e2 - e1);
    dstype time = common.time;

    abi.HdgSource(f, f_udg, f_wdg, xg, udg, odg, wdg, app.uinf,
                  app.physicsparam, time, common.modelnumber, numPoints, nc,
                  ncu, nd, ncx, nco, ncw);
}

void SourcewDriver(dstype* f, dstype* f_udg, dstype* f_wdg,
                   const dstype* xg, const dstype* udg, const dstype* odg,
                   const dstype* wdg, ExasimDriverABI& abi, meshstruct& mesh,
                   masterstruct& master, appstruct& app, solstruct& sol,
                   tempstruct& temp, commonstruct& common, Int nge, Int e1,
                   Int e2, Int backend)
{
    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nd = common.nd;
    Int numPoints = nge * (e2 - e1);
    dstype time = common.time;

    abi.HdgSourcew(f, f_udg, f_wdg, xg, udg, odg, wdg, app.uinf,
                   app.physicsparam, time, common.modelnumber, numPoints, nc,
                   ncu, nd, ncx, nco, ncw);
}

void SourcewDriver(dstype* f, dstype* f_wdg, const dstype* xg,
                   const dstype* udg, const dstype* odg, const dstype* wdg,
                   ExasimDriverABI& abi, meshstruct& mesh,
                   masterstruct& master, appstruct& app, solstruct& sol,
                   tempstruct& temp, commonstruct& common, Int nge, Int e1,
                   Int e2, Int backend)
{
    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nd = common.nd;
    Int numPoints = nge * (e2 - e1);
    dstype time = common.time;

    abi.HdgSourcewonly(f, f_wdg, xg, udg, odg, wdg, app.uinf,
                       app.physicsparam, time, common.modelnumber, numPoints,
                       nc, ncu, nd, ncx, nco, ncw);
}

void EosDriver(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xg,
               const dstype* udg, const dstype* odg, const dstype* wdg,
               ExasimDriverABI& abi, meshstruct& mesh, masterstruct& master,
               appstruct& app, solstruct& sol, tempstruct& temp,
               commonstruct& common, Int nge, Int e1, Int e2, Int backend)
{
    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nd = common.nd;
    Int numPoints = nge * (e2 - e1);
    dstype time = common.time;

    abi.HdgEoS(f, f_udg, f_wdg, xg, udg, odg, wdg, app.uinf,
               app.physicsparam, time, common.modelnumber, numPoints, nc, ncu,
               nd, ncx, nco, ncw);
}

void FbouDriver(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg,
                dstype* xg, const dstype* udg, const dstype* odg,
                const dstype* wdg, dstype* uhg, const dstype* nl,
                ExasimDriverABI& abi, meshstruct& mesh, masterstruct& master,
                appstruct& app, solstruct& sol, tempstruct& temp,
                commonstruct& common, Int nga, Int ib, Int backend)
{
    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nd = common.nd;
    Int numPoints = nga;
    dstype time = common.time;

    abi.HdgFbou(f, f_udg, f_wdg, f_uhg, xg, udg, odg, wdg, uhg, nl, app.tau,
                app.uinf, app.physicsparam, time, common.modelnumber, ib,
                numPoints, nc, ncu, nd, ncx, nco, ncw);
}

void FbouDriver(dstype* f, dstype* xg, const dstype* udg,
                const dstype* odg, const dstype* wdg, dstype* uhg,
                const dstype* nl, ExasimDriverABI& abi, meshstruct& mesh,
                masterstruct& master, appstruct& app, solstruct& sol,
                tempstruct& temp, commonstruct& common, Int nga, Int ib,
                Int backend)
{
    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nd = common.nd;
    Int numPoints = nga;
    dstype time = common.time;

    abi.HdgFbouonly(f, xg, udg, odg, wdg, uhg, nl, app.tau, app.uinf,
                    app.physicsparam, time, common.modelnumber, ib, numPoints,
                    nc, ncu, nd, ncx, nco, ncw);
}

void FintDriver(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg,
                dstype* xg, const dstype* udg, const dstype* odg,
                const dstype* wdg, dstype* uhg, const dstype* nl,
                ExasimDriverABI& abi, meshstruct& mesh, masterstruct& master,
                appstruct& app, solstruct& sol, tempstruct& temp,
                commonstruct& common, Int nga, Int ib, Int backend)
{
    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nd = common.nd;
    Int numPoints = nga;
    dstype time = common.time;

    abi.HdgFint(f, f_udg, f_wdg, f_uhg, xg, udg, odg, wdg, uhg, nl, app.tau,
                app.uinf, app.physicsparam, time, common.modelnumber, ib,
                numPoints, nc, ncu, nd, ncx, nco, ncw);
}

void FintDriver(dstype* f, dstype* xg, const dstype* udg,
                const dstype* odg, const dstype* wdg, dstype* uhg,
                const dstype* nl, ExasimDriverABI& abi, meshstruct& mesh,
                masterstruct& master, appstruct& app, solstruct& sol,
                tempstruct& temp, commonstruct& common, Int nga, Int ib,
                Int backend)
{
    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nd = common.nd;
    Int numPoints = nga;
    dstype time = common.time;

    abi.HdgFintonly(f, xg, udg, odg, wdg, uhg, nl, app.tau, app.uinf,
                    app.physicsparam, time, common.modelnumber, ib, numPoints,
                    nc, ncu, nd, ncx, nco, ncw);
}

void FextDriver(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg,
                dstype* xg, const dstype* udg, const dstype* odg,
                const dstype* wdg, dstype* uhg, const dstype* nl,
                const dstype* uext, ExasimDriverABI& abi, meshstruct& mesh,
                masterstruct& master, appstruct& app, solstruct& sol,
                tempstruct& temp, commonstruct& common, Int nga, Int ib,
                Int backend)
{
    abi.HdgFext(f, f_udg, f_wdg, f_uhg, xg, udg, odg, wdg, uhg, nl, uext,
                app.tau, app.uinf, app.physicsparam, common.time,
                common.modelnumber, ib, nga, common.nc, common.ncu, common.nd,
                common.ncx, common.nco, common.ncw);
}

void FextDriver(dstype* f, dstype* xg, const dstype* udg,
                const dstype* odg, const dstype* wdg, dstype* uhg,
                const dstype* nl, const dstype* uext,
                ExasimDriverABI& abi, meshstruct& mesh, masterstruct& master, appstruct& app,
                solstruct& sol, tempstruct& temp, commonstruct& common,
                Int nga, Int ib, Int backend)
{
    abi.HdgFextonly(f, xg, udg, odg, wdg, uhg, nl, uext, app.tau, app.uinf,
                    app.physicsparam, common.time, common.modelnumber, ib, nga,
                    common.nc, common.ncu, common.nd, common.ncx, common.nco,
                    common.ncw);
}

void FhatDriver(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg,
                const dstype* xg, dstype* udg, const dstype* odg,
                const dstype* wdg, const dstype* uhg, dstype* nl,
                ExasimDriverABI& abi, meshstruct& mesh, masterstruct& master,
                appstruct& app, solstruct& sol, tempstruct& temp,
                commonstruct& common, Int nga, Int backend)
{
    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nd = common.nd;
    Int numPoints = nga;
    Int M = numPoints * ncu;
    Int N = numPoints * ncu * nd;
    dstype time = common.time;

    ArrayCopy(f_uhg, udg, numPoints * ncu);
    ArrayCopy(udg, uhg, numPoints * ncu);

    abi.HdgFlux(f, f_udg, f_wdg, xg, udg, odg, wdg, app.uinf,
                app.physicsparam, time, common.modelnumber, numPoints, nc,
                ncu, nd, ncx, nco, ncw);

    ArrayCopy(udg, f_uhg, numPoints * ncu);

    FluxDotNormal(f, f, nl, M, numPoints, nd);

    for (int n = 0; n < nc; n++) {
        FluxDotNormal(&f_udg[M * n], &f_udg[N * n], nl, M, numPoints, nd);
    }

    if ((ncw > 0) & (common.wave == 0)) {
        for (int n = 0; n < ncw; n++) {
            FluxDotNormal(&f_wdg[M * n], &f_wdg[N * n], nl, M, numPoints, nd);
        }
    }

    ArrayCopy(f_uhg, f_udg, numPoints * ncu * ncu);
    ArraySetValue(f_udg, zero, numPoints * ncu * ncu);
    AddStabilization1(f, f_udg, f_uhg, udg, uhg, app.tau, M, numPoints);
}

void FhatDriver(dstype* f, dstype* u, const dstype* xg, dstype* udg,
                const dstype* odg, const dstype* wdg, const dstype* uhg,
                dstype* nl, ExasimDriverABI& abi, meshstruct& mesh,
                masterstruct& master, appstruct& app, solstruct& sol,
                tempstruct& temp, commonstruct& common, Int nga, Int backend)
{
    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nd = common.nd;
    Int numPoints = nga;
    Int M = numPoints * ncu;
    dstype time = common.time;

    ArrayCopy(u, udg, numPoints * ncu);
    ArrayCopy(udg, uhg, numPoints * ncu);

    abi.KokkosFlux(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time,
                   common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);

    ArrayCopy(udg, u, numPoints * ncu);

    FluxDotNormal(f, f, nl, M, numPoints, nd);
    AddStabilization1(f, udg, uhg, app.tau, M);
}

void VisScalarsDriver(dstype* f, const dstype* xg, const dstype* udg,
                      const dstype* odg, const dstype* wdg,
                      ExasimDriverABI& abi, meshstruct& mesh,
                      masterstruct& master, appstruct& app, solstruct& sol,
                      tempstruct& temp, commonstruct& common, Int nge, Int e1,
                      Int e2, Int backend)
{
    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nd = common.nd;
    Int numPoints = nge * (e2 - e1);
    dstype time = common.time;

    abi.KokkosVisScalars(f, xg, udg, odg, wdg, app.uinf, app.physicsparam,
                         time, common.modelnumber, numPoints, nc, ncu, nd, ncx,
                         nco, ncw);
}

void VisVectorsDriver(dstype* f, const dstype* xg, const dstype* udg,
                      const dstype* odg, const dstype* wdg,
                      ExasimDriverABI& abi, meshstruct& mesh,
                      masterstruct& master, appstruct& app, solstruct& sol,
                      tempstruct& temp, commonstruct& common, Int nge, Int e1,
                      Int e2, Int backend)
{
    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nd = common.nd;
    Int numPoints = nge * (e2 - e1);
    dstype time = common.time;

    abi.KokkosVisVectors(f, xg, udg, odg, wdg, app.uinf, app.physicsparam,
                         time, common.modelnumber, numPoints, nc, ncu, nd, ncx,
                         nco, ncw);
}

void VisTensorsDriver(dstype* f, const dstype* xg, const dstype* udg,
                      const dstype* odg, const dstype* wdg,
                      ExasimDriverABI& abi, meshstruct& mesh,
                      masterstruct& master, appstruct& app, solstruct& sol,
                      tempstruct& temp, commonstruct& common, Int nge, Int e1,
                      Int e2, Int backend)
{
    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nd = common.nd;
    Int numPoints = nge * (e2 - e1);
    dstype time = common.time;

    abi.KokkosVisTensors(f, xg, udg, odg, wdg, app.uinf, app.physicsparam,
                         time, common.modelnumber, numPoints, nc, ncu, nd, ncx,
                         nco, ncw);
}

void QoIvolumeDriver(dstype* f, const dstype* xg, const dstype* udg,
                     const dstype* odg, const dstype* wdg,
                     ExasimDriverABI& abi, meshstruct& mesh,
                     masterstruct& master, appstruct& app, solstruct& sol,
                     tempstruct& temp, commonstruct& common, Int nge, Int e1,
                     Int e2, Int backend)
{
    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nd = common.nd;
    Int numPoints = nge * (e2 - e1);
    dstype time = common.time;

    abi.KokkosQoIvolume(f, xg, udg, odg, wdg, app.uinf, app.physicsparam,
                        time, common.modelnumber, numPoints, nc, ncu, nd, ncx,
                        nco, ncw);
}

void QoIboundaryDriver(dstype* fb, const dstype* xg, const dstype* udg,
                       const dstype* odg, const dstype* wdg,
                       const dstype* uhg, const dstype* nl,
                       ExasimDriverABI& abi, meshstruct& mesh,
                       masterstruct& master, appstruct& app, solstruct& sol,
                       tempstruct& temp, commonstruct& common, Int ngf, Int f1,
                       Int f2, Int ib, Int backend)
{
    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nd = common.nd;
    Int numPoints = ngf * (f2 - f1);
    dstype time = common.time;

    abi.KokkosQoIboundary(fb, xg, udg, odg, wdg, uhg, nl, app.tau, app.uinf,
                          app.physicsparam, time, common.modelnumber, ib,
                          numPoints, nc, ncu, nd, ncx, nco, ncw);
}

#endif
