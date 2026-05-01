// SPDX-License-Identifier: see LICENSE
//
// backend/Preprocessing/buildstructs.hpp
//
// Direct in-memory constructors for the runtime data structs
// (`appstruct`, `masterstruct`, `meshstruct`, `solstruct`).
//
// Today these are populated by `readappstruct(filename, ...)` etc.,
// which deserialize from `app.bin`/`master.bin`/`mesh.bin`/`sol.bin`
// written by `writepde`/`writemaster`/`writemesh`/`writesol` in
// CPreprocessing. The binaries exist solely to cross the
// CPreprocessing → CSolution boundary.
//
// `buildAppStruct(const PDE&)` etc. populate the same runtime structs
// by reading typed values straight off `PDE` / `Mesh` / `Master` /
// `DMD` — no serialization roundtrip. Layout-equivalent to the read
// functions: every `app.X = readY(...)` is replaced with
// `app.X = malloc + copy from pde.Y`.
//
// Used by the in-memory `Preprocessed` handoff (HOT.7.3); the
// file-based path (`readappstruct(filename, ...)`) is kept and
// exercised by every legacy `_codegen` / `_template` binary.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "../Common/common.h"            // appstruct, masterstruct, meshstruct, solstruct
#include "structs.hpp"                   // PDE, Mesh, Master, DMD

namespace exasim {

// ---------------------------------------------------------------
// Helpers — parallel the (Int*) malloc + copy and (dstype*) malloc
// + copy patterns used by readappstruct etc.
// ---------------------------------------------------------------
inline Int* mallocIntArray(const std::vector<double>& src)
{
    Int n = (Int)src.size();
    if (n <= 0) return nullptr;
    Int* a = (Int*)std::malloc(sizeof(Int) * n);
    for (Int i = 0; i < n; ++i) a[i] = (Int)std::round(src[i]);
    return a;
}

inline Int* mallocIntArrayShifted(const std::vector<double>& src, double shift)
{
    Int n = (Int)src.size();
    if (n <= 0) return nullptr;
    Int* a = (Int*)std::malloc(sizeof(Int) * n);
    for (Int i = 0; i < n; ++i) a[i] = (Int)std::round(src[i] + shift);
    return a;
}

inline dstype* mallocDoubleArray(const std::vector<double>& src)
{
    Int n = (Int)src.size();
    if (n <= 0) return nullptr;
    dstype* a = (dstype*)std::malloc(sizeof(dstype) * n);
    for (Int i = 0; i < n; ++i) a[i] = (dstype)src[i];
    return a;
}

inline Int* mallocIntArrayN(const Int* src, Int n)
{
    if (n <= 0) return nullptr;
    Int* a = (Int*)std::malloc(sizeof(Int) * n);
    for (Int i = 0; i < n; ++i) a[i] = src[i];
    return a;
}

inline dstype* mallocDoubleArrayN(const dstype* src, Int n)
{
    if (n <= 0) return nullptr;
    dstype* a = (dstype*)std::malloc(sizeof(dstype) * n);
    for (Int i = 0; i < n; ++i) a[i] = src[i];
    return a;
}

// ---------------------------------------------------------------
// buildAppStruct(pde) — mirrors readappstruct() / writepde() pair.
//
// Source mapping (from writepde at backend/Preprocessing/readpdeapp.hpp):
//   ndims[0]   = pde.mpiprocs
//   ndims[1]   = pde.nd
//   ndims[5]   = pde.nc            ndims[6]   = pde.ncu
//   ndims[7]   = pde.ncq           ndims[8]   = pde.ncp
//   ndims[9]   = pde.ncv           ndims[10]  = pde.nch
//   ndims[11]  = pde.ncx           ndims[12]  = pde.nce
//   ndims[13]  = pde.ncw           ndims[14]  = pde.nsca
//   ndims[15]  = pde.nvec          ndims[16]  = pde.nten
//   ndims[17]  = pde.nsurf         ndims[18]  = pde.nvqoi
//   nsize[0]   = ndims.size() (=40)
//   nsize[k]   = .size() of k-th vector below
//   payload    = flag, problem, externalparam, dt, factor,
//                physicsparam, solversparam, tau, stgdata, stgparam,
//                stgib, vindx (shifted by -1), dae_dt, interfaceFluxmap,
//                avparam (= avparam1 ++ avparam2)
//
// app.uinf is sourced from pde.externalparam (this is the legacy ABI;
// see readappstruct line 81: app.uinf <- nsize[3] = externalparam.size()).
// ---------------------------------------------------------------
inline appstruct buildAppStruct(const PDE& pde)
{
    appstruct app{};

    // ---- ndims (40 entries, dimensions) ----
    constexpr Int kNDims = 40;
    app.ndims = (Int*)std::calloc(kNDims, sizeof(Int));
    app.ndims[0]  = pde.mpiprocs;
    app.ndims[1]  = pde.nd;
    app.ndims[5]  = pde.nc;
    app.ndims[6]  = pde.ncu;
    app.ndims[7]  = pde.ncq;
    app.ndims[8]  = pde.ncp;
    app.ndims[9]  = pde.ncv;
    app.ndims[10] = pde.nch;
    app.ndims[11] = pde.ncx;
    app.ndims[12] = pde.nce;
    app.ndims[13] = pde.ncw;
    app.ndims[14] = pde.nsca;
    app.ndims[15] = pde.nvec;
    app.ndims[16] = pde.nten;
    app.ndims[17] = pde.nsurf;
    app.ndims[18] = pde.nvqoi;

    // ---- nsize (30 entries, sub-array sizes) ----
    constexpr Int kNSize = 30;
    std::vector<double> avparam;
    avparam.insert(avparam.end(), pde.avparam1.begin(), pde.avparam1.end());
    avparam.insert(avparam.end(), pde.avparam2.begin(), pde.avparam2.end());

    app.nsize = (Int*)std::calloc(kNSize, sizeof(Int));
    app.nsize[0]  = kNDims;
    app.nsize[1]  = (Int)pde.flag.size();
    app.nsize[2]  = (Int)pde.problem.size();
    app.nsize[3]  = (Int)pde.externalparam.size();
    app.nsize[4]  = (Int)pde.dt.size();
    app.nsize[5]  = (Int)pde.factor.size();
    app.nsize[6]  = (Int)pde.physicsparam.size();
    app.nsize[7]  = (Int)pde.solversparam.size();
    app.nsize[8]  = (Int)pde.tau.size();
    app.nsize[9]  = (Int)pde.stgdata.size();
    app.nsize[10] = (Int)pde.stgparam.size();
    app.nsize[11] = (Int)pde.stgib.size();
    app.nsize[12] = (Int)pde.vindx.size();
    app.nsize[13] = (Int)pde.dae_dt.size();
    app.nsize[14] = (Int)pde.interfaceFluxmap.size();
    app.nsize[15] = (Int)avparam.size();

    app.lsize = (Int*)std::malloc(sizeof(Int));
    app.lsize[0] = kNSize;

    // ---- typed payloads ----
    app.flag             = mallocIntArray(pde.flag);
    app.problem          = mallocIntArray(pde.problem);
    app.uinf             = mallocDoubleArray(pde.externalparam);
    app.dt               = mallocDoubleArray(pde.dt);
    app.factor           = mallocDoubleArray(pde.factor);
    app.physicsparam     = mallocDoubleArray(pde.physicsparam);
    app.solversparam     = mallocDoubleArray(pde.solversparam);
    app.tau              = mallocDoubleArray(pde.tau);
    app.stgdata          = mallocDoubleArray(pde.stgdata);
    app.stgparam         = mallocDoubleArray(pde.stgparam);
    app.stgib            = mallocIntArray(pde.stgib);
    app.vindx            = mallocIntArrayShifted(pde.vindx, -1.0);   // legacy `-1` on write, restored on read
    app.dae_dt           = mallocDoubleArray(pde.dae_dt);
    {
        // pde.interfaceFluxmap is std::vector<int>; the on-disk write
        // serializes it as doubles via writeiarraytodouble, then read
        // converts back. Direct-build skips both steps.
        Int n = (Int)pde.interfaceFluxmap.size();
        if (n > 0) {
            app.interfacefluxmap = (Int*)std::malloc(sizeof(Int) * n);
            for (Int i = 0; i < n; ++i)
                app.interfacefluxmap[i] = (Int)pde.interfaceFluxmap[i];
        }
    }
    app.avparam          = mallocDoubleArray(avparam);

    // ---- size fields (mirrors readappstruct lines 95-109) ----
    app.szflag             = app.nsize[1];
    app.szproblem          = app.nsize[2];
    app.szuinf             = app.nsize[3];
    app.szdt               = app.nsize[4];
    app.szfactor           = app.nsize[5];
    app.szphysicsparam     = app.nsize[6];
    app.szsolversparam     = app.nsize[7];
    app.sztau              = app.nsize[8];
    app.szstgdata          = app.nsize[9];
    app.szstgparam         = app.nsize[10];
    app.szstgib            = app.nsize[11];
    app.szvindx            = app.nsize[12];
    app.szdae_dt           = app.nsize[13];
    app.szinterfacefluxmap = app.nsize[14];
    app.szavparam          = app.nsize[15];

    // ---- derived fc_u/fc_q/fc_w (mirrors readappstruct lines 134-168) ----
    Int ncu = app.ndims[6];
    Int ncq = app.ndims[7];
    Int ncw = app.ndims[13];
    if (ncu > 0) {
        app.fc_u     = (dstype*)std::malloc(sizeof(dstype) * ncu);
        app.dtcoef_u = (dstype*)std::malloc(sizeof(dstype) * ncu);
        for (Int i = 0; i < ncu; ++i) { app.fc_u[i] = 1.0; app.dtcoef_u[i] = 1.0; }
        app.szfc_u = ncu; app.szdtcoef_u = ncu;
    }
    if (ncq > 0) {
        app.fc_q     = (dstype*)std::malloc(sizeof(dstype) * ncq);
        app.dtcoef_q = (dstype*)std::malloc(sizeof(dstype) * ncq);
        for (Int i = 0; i < ncq; ++i) { app.fc_q[i] = 1.0; app.dtcoef_q[i] = 1.0; }
        app.szfc_q = ncu; app.szdtcoef_q = ncu;   // note: legacy uses ncu (see readappstruct line 156)
    }
    if (ncw > 0) {
        app.fc_w     = (dstype*)std::malloc(sizeof(dstype) * ncw);
        app.dtcoef_w = (dstype*)std::malloc(sizeof(dstype) * ncw);
        for (Int i = 0; i < ncw; ++i) { app.fc_w[i] = 1.0; app.dtcoef_w[i] = 1.0; }
        app.szfc_w = ncw; app.szdtcoef_w = ncw;
    }

    // The runtime sets app.porder + app.comm in readInput AFTER the
    // app struct is built; we leave them nullptr here for parity with
    // readappstruct's behavior (porder is set in readInput).

    return app;
}

// ---------------------------------------------------------------
// buildMasterStruct(master) — mirrors readmasterstruct() / writemaster() pair.
//
// `Master` (preprocessing struct) holds 1D/2D/3D shape/quadrature
// arrays as `vector<double>`. `masterstruct` (runtime) holds the
// same data as `dstype*` arrays for solver consumption. This is a
// straight copy.
// ---------------------------------------------------------------
inline masterstruct buildMasterStruct(const Master& master)
{
    masterstruct ms{};

    // ---- ndims (20 entries) ----
    constexpr Int kNDims = 20;
    ms.ndims = (Int*)std::calloc(kNDims, sizeof(Int));
    ms.ndims[0]  = master.nd;
    ms.ndims[1]  = master.elemtype;
    ms.ndims[2]  = master.nodetype;
    ms.ndims[3]  = master.porder;
    ms.ndims[4]  = master.pgauss;
    ms.ndims[5]  = master.npe;
    ms.ndims[6]  = master.npf;
    ms.ndims[7]  = master.nge;
    ms.ndims[8]  = master.ngf;
    ms.ndims[9]  = master.np1d;
    ms.ndims[10] = master.ng1d;

    // ---- nsize (22 entries) ----
    constexpr Int kNSize = 22;
    ms.nsize = (Int*)std::calloc(kNSize, sizeof(Int));
    ms.nsize[0]  = kNDims;
    ms.nsize[1]  = (Int)master.shapegt.size();
    ms.nsize[2]  = (Int)master.shapegw.size();
    ms.nsize[3]  = (Int)master.shapfgt.size();
    ms.nsize[4]  = (Int)master.shapfgw.size();
    ms.nsize[5]  = (Int)master.shapent.size();
    ms.nsize[6]  = (Int)master.shapen.size();
    ms.nsize[7]  = (Int)master.shapfnt.size();
    ms.nsize[8]  = (Int)master.shapfn.size();
    ms.nsize[9]  = (Int)master.xpe.size();
    ms.nsize[10] = (Int)master.gpe.size();
    ms.nsize[11] = (Int)master.gwe.size();
    ms.nsize[12] = (Int)master.xpf.size();
    ms.nsize[13] = (Int)master.gpf.size();
    ms.nsize[14] = (Int)master.gwf.size();
    ms.nsize[15] = (Int)master.shap1dgt.size();
    ms.nsize[16] = (Int)master.shap1dgw.size();
    ms.nsize[17] = (Int)master.shap1dnt.size();
    ms.nsize[18] = (Int)master.shap1dn.size();   // → masterstruct::shap1dnl (legacy name)
    ms.nsize[19] = (Int)master.xp1d.size();
    ms.nsize[20] = (Int)master.gp1d.size();
    ms.nsize[21] = (Int)master.gw1d.size();

    ms.lsize = (Int*)std::malloc(sizeof(Int));
    ms.lsize[0] = kNSize;

    // ---- typed payloads ----
    ms.shapegt   = mallocDoubleArray(master.shapegt);
    ms.shapegw   = mallocDoubleArray(master.shapegw);
    ms.shapfgt   = mallocDoubleArray(master.shapfgt);
    ms.shapfgw   = mallocDoubleArray(master.shapfgw);
    ms.shapent   = mallocDoubleArray(master.shapent);
    ms.shapen    = mallocDoubleArray(master.shapen);
    ms.shapfnt   = mallocDoubleArray(master.shapfnt);
    ms.shapfn    = mallocDoubleArray(master.shapfn);
    ms.xpe       = mallocDoubleArray(master.xpe);
    ms.gpe       = mallocDoubleArray(master.gpe);
    ms.gwe       = mallocDoubleArray(master.gwe);
    ms.xpf       = mallocDoubleArray(master.xpf);
    ms.gpf       = mallocDoubleArray(master.gpf);
    ms.gwf       = mallocDoubleArray(master.gwf);
    ms.shap1dgt  = mallocDoubleArray(master.shap1dgt);
    ms.shap1dgw  = mallocDoubleArray(master.shap1dgw);
    ms.shap1dnt  = mallocDoubleArray(master.shap1dnt);
    ms.shap1dnl  = mallocDoubleArray(master.shap1dn);   // name change: shap1dn → shap1dnl
    ms.xp1d      = mallocDoubleArray(master.xp1d);
    ms.gp1d      = mallocDoubleArray(master.gp1d);
    ms.gw1d      = mallocDoubleArray(master.gw1d);

    // ---- size fields ----
    ms.szshapegt   = ms.nsize[1];
    ms.szshapegw   = ms.nsize[2];
    ms.szshapfgt   = ms.nsize[3];
    ms.szshapfgw   = ms.nsize[4];
    ms.szshapent   = ms.nsize[5];
    ms.szshapen    = ms.nsize[6];
    ms.szshapfnt   = ms.nsize[7];
    ms.szshapfn    = ms.nsize[8];
    ms.szxpe       = ms.nsize[9];
    ms.szgpe       = ms.nsize[10];
    ms.szgwe       = ms.nsize[11];
    ms.szxpf       = ms.nsize[12];
    ms.szgpf       = ms.nsize[13];
    ms.szgwf       = ms.nsize[14];
    ms.szshap1dgt  = ms.nsize[15];
    ms.szshap1dgw  = ms.nsize[16];
    ms.szshap1dnt  = ms.nsize[17];
    ms.szshap1dnl  = ms.nsize[18];
    ms.szxp1d      = ms.nsize[19];
    ms.szgp1d      = ms.nsize[20];
    ms.szgw1d      = ms.nsize[21];

    return ms;
}

// ---------------------------------------------------------------
// buildMeshStruct(mesh, master, dmd, bf, hybrid, mpiprocs)
//   — mirrors readmeshstruct() / writemesh() pair.
//
// Sets only the fields actually written by writemesh(). The runtime
// `meshstruct` has many more fields (facecon, eblks, fblks, e2f,
// f2e, f2f, …) that are populated later by `buildConn(mesh, sol,
// app, master, …)`; this builder is the equivalent of the
// readmeshstruct() step that loads ti / boundaryConditions /
// intepartpts and stops there. Caller (`finalizePreprocessed`)
// invokes buildConn afterward.
//
// MPI fields (nbsd / elemsend / elemrecv / …) are populated when
// `mpiprocs > 1` from the DMD struct, mirroring writemesh's
// `if (pde.mpiprocs>1)` branch.
//
// Outputs `ti_out` for the caller — needed by buildConn().
// ---------------------------------------------------------------
inline meshstruct buildMeshStruct(const Mesh& mesh, const Master& master,
                                  const DMD& dmd, const std::vector<int>& bf,
                                  int hybrid, int mpiprocs,
                                  std::vector<int>& ti_out)
{
    meshstruct ms{};

    const int ne   = (int)dmd.elempart.size();
    const int nve  = mesh.nve;
    const int nfe  = mesh.nfe;

    // ---- ndims (20 entries; only a handful are actually used) ----
    constexpr Int kNDims = 20;
    ms.ndims = (Int*)std::calloc(kNDims, sizeof(Int));
    ms.ndims[0] = mesh.dim;
    ms.ndims[1] = ne;
    {
        int maxt = 0;
        for (int v : mesh.t) if (v > maxt) maxt = v;
        ms.ndims[3] = maxt + 1;
    }
    ms.ndims[4] = nfe;

    // ---- nsize (50 entries) ----
    constexpr Int kNSize = 50;
    ms.nsize = (Int*)std::calloc(kNSize, sizeof(Int));
    ms.nsize[0] = kNDims;
    if (mpiprocs > 1) {
        ms.nsize[4] = (Int)dmd.nbsd.size();
        ms.nsize[5] = (Int)dmd.elemsend.size();
        ms.nsize[6] = (Int)dmd.elemrecv.size();
        ms.nsize[7] = (Int)dmd.elemsendpts.size();
        ms.nsize[8] = (Int)dmd.elemrecvpts.size();
    }
    ms.nsize[9]  = (Int)dmd.elempart.size();
    ms.nsize[10] = (Int)dmd.elempartpts.size();
    ms.nsize[23] = (Int)master.perm.size();
    ms.nsize[24] = (Int)bf.size();
    ms.nsize[25] = (Int)mesh.cartGridPart.size();
    ms.nsize[26] = (Int)(ne * nve);
    ms.nsize[27] = (Int)mesh.boundaryConditions.size();
    ms.nsize[28] = (Int)dmd.intepartpts.size();

    ms.lsize = (Int*)std::malloc(sizeof(Int));
    ms.lsize[0] = kNSize;

    // ---- MPI-only payloads (mirrors writemesh's pde.mpiprocs>1 branch) ----
    if (mpiprocs > 1) {
        // elemsend / elemrecv each store rows [sender, recv_local_idx, sender_global_idx]
        // — writemesh extracts column [1] only.
        std::vector<int> elemsend(dmd.elemsend.size());
        std::vector<int> elemrecv(dmd.elemrecv.size());
        for (size_t i = 0; i < dmd.elemsend.size(); ++i) elemsend[i] = dmd.elemsend[i][1];
        for (size_t i = 0; i < dmd.elemrecv.size(); ++i) elemrecv[i] = dmd.elemrecv[i][1];
        ms.nbsd        = mallocIntArrayN(dmd.nbsd.data(),       (Int)dmd.nbsd.size());
        ms.elemsend    = mallocIntArrayN(elemsend.data(),       (Int)elemsend.size());
        ms.elemrecv    = mallocIntArrayN(elemrecv.data(),       (Int)elemrecv.size());
        ms.elemsendpts = mallocIntArrayN(dmd.elemsendpts.data(),(Int)dmd.elemsendpts.size());
        ms.elemrecvpts = mallocIntArrayN(dmd.elemrecvpts.data(),(Int)dmd.elemrecvpts.size());
        ms.sznbsd        = ms.nsize[4];
        ms.szelemsend    = ms.nsize[5];
        ms.szelemrecv    = ms.nsize[6];
        ms.szelemsendpts = ms.nsize[7];
        ms.szelemrecvpts = ms.nsize[8];
    }
    ms.elempart    = mallocIntArrayN(dmd.elempart.data(),    (Int)dmd.elempart.size());
    ms.elempartpts = mallocIntArrayN(dmd.elempartpts.data(), (Int)dmd.elempartpts.size());
    ms.szelempart    = ms.nsize[9];
    ms.szelempartpts = ms.nsize[10];

    // ---- HDG-only payloads ----
    if (hybrid > 0) {
        ms.perm         = mallocIntArrayN(master.perm.data(),       (Int)master.perm.size());
        ms.bf           = mallocIntArrayN(bf.data(),                (Int)bf.size());
        ms.cartgridpart = mallocIntArrayN(mesh.cartGridPart.data(), (Int)mesh.cartGridPart.size());
        ms.szperm        = ms.nsize[23];
        ms.szbf          = ms.nsize[24];
        ms.szcartgridpart= ms.nsize[25];
    }

    // ---- ti (mesh.t reordered by elempart) — emitted to caller ----
    ti_out.assign((size_t)nve * ne, 0);
    for (int j = 0; j < ne; ++j) {
        int col = dmd.elempart[j];
        for (int i = 0; i < nve; ++i)
            ti_out[i + j*nve] = mesh.t[i + col*nve];
    }

    // ---- boundaryConditions, intepartpts ----
    ms.boundaryConditions = mallocIntArrayN(mesh.boundaryConditions.data(),
                                            (Int)mesh.boundaryConditions.size());
    ms.intepartpts        = mallocIntArrayN(dmd.intepartpts.data(),
                                            (Int)dmd.intepartpts.size());

    return ms;
}

// ---------------------------------------------------------------
// buildSolStruct(mesh, master, pde, ne, ncudg) — mirrors
// readsolstruct() / writesol() pair, *without* the udg / odg / wdg
// initial-condition computation. The legacy readsolstruct calls
// cpuInituDriver / cpuInitodgDriver / cpuInitwdgDriver when the
// corresponding nsize[k] is zero. That's hard to do here because
// app.flag isn't yet built — the orchestrator (see
// `buildPreprocessed` below) calls those drivers after solstruct is
// built but before buildConn().
//
// `ne_local` is `dmd.elempart.size()` (per-rank element count).
// `ncudg`   is the legacy heuristic from writesol: `mesh.udg.size()
// / (mesh.ne * master.npe)` if udg is nonempty, else `pde.nc`.
// ---------------------------------------------------------------
inline solstruct buildSolStruct(const Mesh& mesh, const Master& master,
                                const PDE& pde, const DMD& dmd)
{
    solstruct sol{};

    const int ne_local = (int)dmd.elempart.size();
    const int ne_global= mesh.ne;
    const int npe      = master.npe;

    // ---- ndims (12 entries) ----
    constexpr Int kNDims = 12;
    sol.ndims = (Int*)std::calloc(kNDims, sizeof(Int));
    sol.ndims[0]  = ne_local;
    sol.ndims[2]  = mesh.nfe;
    sol.ndims[3]  = master.npe;
    sol.ndims[4]  = master.npf;
    sol.ndims[5]  = pde.nc;
    sol.ndims[6]  = pde.ncu;
    sol.ndims[7]  = pde.ncq;
    sol.ndims[8]  = pde.ncw;
    sol.ndims[9]  = pde.ncv;
    sol.ndims[10] = pde.nch;
    sol.ndims[11] = pde.ncx;

    // ---- nsize (20 entries; xdg=1, udg=2, vdg=3, wdg=4, uhat=5, xcg=6) ----
    constexpr Int kNSize = 20;
    sol.nsize = (Int*)std::calloc(kNSize, sizeof(Int));
    sol.nsize[0] = kNDims;
    if (!mesh.xdg.empty())
        sol.nsize[1] = (Int)(npe * mesh.dim * ne_local);
    int ncudg = 0;
    if (!mesh.udg.empty()) {
        ncudg = (int)mesh.udg.size() / (ne_global * npe);
        if (ncudg == pde.ncu || ncudg == pde.nc)
            sol.nsize[2] = (Int)(npe * ncudg * ne_local);
    }
    if (!mesh.vdg.empty()) sol.nsize[3] = (Int)(npe * pde.ncv * ne_local);
    if (!mesh.wdg.empty()) sol.nsize[4] = (Int)(npe * pde.ncw * ne_local);
    if (!mesh.uhat.empty()) sol.nsize[5] = (Int)mesh.uhat.size();

    sol.lsize = (Int*)std::malloc(sizeof(Int));
    sol.lsize[0] = kNSize;

    // ---- xdg (always written if mesh.xdg present, selected by elempart) ----
    if (!mesh.xdg.empty()) {
        const int row = npe * mesh.dim;
        sol.xdg = (dstype*)std::malloc(sizeof(dstype) * row * ne_local);
        for (int j = 0; j < ne_local; ++j) {
            int col = dmd.elempart[j];
            for (int i = 0; i < row; ++i)
                sol.xdg[i + j*row] = (dstype)mesh.xdg[i + col*row];
        }
        sol.szxdg = sol.nsize[1];
    }

    // ---- udg, vdg, wdg, uhat: only the *file-present* paths.
    //     The legacy "compute initial via cpuInituDriver" branch
    //     happens in `finalizePreprocessed` after app/master/sol
    //     are all in hand. ----
    if (!mesh.udg.empty()) {
        const int row = npe * ncudg;
        sol.udg = (dstype*)std::malloc(sizeof(dstype) * row * ne_local);
        for (int j = 0; j < ne_local; ++j) {
            int col = dmd.elempart[j];
            for (int i = 0; i < row; ++i)
                sol.udg[i + j*row] = (dstype)mesh.udg[i + col*row];
        }
        sol.szudg = sol.nsize[2];
    }
    // The legacy ABI stores `vdg` (auxiliary scalar fields) under
    // `solstruct::odg` — see common.h: odg = auxilary term. Mesh
    // .vdg → solstruct.odg.
    if (!mesh.vdg.empty()) {
        const int row = npe * pde.ncv;
        sol.odg = (dstype*)std::malloc(sizeof(dstype) * row * ne_local);
        for (int j = 0; j < ne_local; ++j) {
            int col = dmd.elempart[j];
            for (int i = 0; i < row; ++i)
                sol.odg[i + j*row] = (dstype)mesh.vdg[i + col*row];
        }
        sol.szodg = sol.nsize[3];
    }
    if (!mesh.wdg.empty()) {
        const int row = npe * pde.ncw;
        sol.wdg = (dstype*)std::malloc(sizeof(dstype) * row * ne_local);
        for (int j = 0; j < ne_local; ++j) {
            int col = dmd.elempart[j];
            for (int i = 0; i < row; ++i)
                sol.wdg[i + j*row] = (dstype)mesh.wdg[i + col*row];
        }
        sol.szwdg = sol.nsize[4];
    }
    if (!mesh.uhat.empty()) {
        sol.uh = (dstype*)std::malloc(sizeof(dstype) * mesh.uhat.size());
        for (size_t i = 0; i < mesh.uhat.size(); ++i)
            sol.uh[i] = (dstype)mesh.uhat[i];
        sol.szuh = sol.nsize[5];
    }

    return sol;
}

// ---------------------------------------------------------------
// Bundle of the four runtime structs ready for CDiscretization /
// CSolution to consume directly. See HOT.7.3 plan.
// ---------------------------------------------------------------
struct Preprocessed {
    appstruct    app;
    masterstruct master;
    meshstruct   mesh;
    solstruct    sol;
    // ti retained for buildConn() inside CDiscretization's setup.
    std::vector<int> ti;
    // HOT.7.4 — when false, CSolution<M> skips opening output
    // bin/qoi files. The data still lives in disc.sol after the
    // solve and can be pulled via host_udg() / host_uhat() / etc.
    bool save_outputs = true;
};

} // namespace exasim
