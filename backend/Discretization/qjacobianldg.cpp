/*
  qjacobianldg.cpp

  This file assembles the element-local LDG Jacobian block

      dQ_K / dU_K = M_K^{-1} ( C_K - E_K d\hat{U}_K/dU_K )

  for each element K.

  Existing backend operators are reused as follows:

  1. qEquationElem(...) computes and stores:
       - res.Minv2 : scalar M_K^{-1}
       - res.C     : scalar M_K^{-1} C_K for each spatial direction

  2. qEquationElemFace(...) computes and stores:
       - res.E     : scalar M_K^{-1} E_K for each spatial direction

  Therefore the actual assembly implemented below uses the already
  premultiplied operators:

      dQ_K / dU_K = C^*_K - E^*_K d\hat{U}_K/dU_K

  where
      C^*_K = M_K^{-1} C_K,
      E^*_K = M_K^{-1} E_K.

  Interior-face trace derivative:

      \hat{u} = 0.5 (u^+ + u^-)

  which is the behavior implemented by GetFaceNodes(..., opts=0). Hence,
  for an interior local face f of element K,

      d\hat{U}_{K,f}/dU_K = 0.5 P_{K,f},

  where P_{K,f} extracts the face DOFs from the element DOFs.

  Boundary-face trace derivative:

  This file leaves the corresponding rows of d\hat{U}_K/dU_K equal to zero.
  Replace qJacobianLDGBoundaryTraceBlock(...) with a future boundary Jacobian
  implementation once UbouDriver provides the needed derivatives.
*/
#ifndef __QJACOBIANLDG
#define __QJACOBIANLDG

inline void qJacobianLDGSetupOperators(solstruct &sol, resstruct &res, appstruct &app,
        masterstruct &master, meshstruct &mesh, tempstruct &tmp,
        commonstruct &common, cublasHandle_t handle, Int backend)
{
    Int npe = common.npe;
    Int npf = common.npf;
    Int nfe = common.nfe;
    Int nd = common.nd;
    Int ne = common.ne;

    if ((res.Minv2 == nullptr) || (res.C == nullptr) ||
        (res.szMinv2 != npe*npe*ne) || (res.szC != npe*npe*ne*nd)) {
        qEquationElem(sol, res, app, master, mesh, tmp, common, handle, backend);
    }

    if ((res.E == nullptr) || (res.szE != npe*npf*nfe*ne*nd)) {
        qEquationElemFace(sol, res, app, master, mesh, tmp, common, handle, backend);
    }
}

inline void qJacobianLDGBoundaryTraceBlock(dstype *DUHf, solstruct &sol, resstruct &res,
        appstruct &app, ExasimDriverABI& driver_abi, masterstruct &master,
        meshstruct &mesh, tempstruct &tmp, commonstruct &common, Int e, Int lf,
        Int face, Int backend)
{
    (void)DUHf;
    (void)sol;
    (void)res;
    (void)app;
    (void)driver_abi;
    (void)master;
    (void)mesh;
    (void)tmp;
    (void)common;
    (void)e;
    (void)lf;
    (void)face;
    (void)backend;
}

inline void qJacobianLDGTrace(dstype *DUH, solstruct &sol, resstruct &res,
        appstruct &app, ExasimDriverABI& driver_abi, masterstruct &master,
        meshstruct &mesh, tempstruct &tmp, commonstruct &common, Int e1, Int e2,
        Int backend)
{
    if (DUH == nullptr)
        return;

    Int npe = common.npe;
    Int npf = common.npf;
    Int nfe = common.nfe;
    Int ncu = common.ncu;
    Int ne = e2 - e1;
    Int ndf = npf*nfe;
    Int nlocu = npe*ncu;
    Int nlocuh = ndf*ncu;

    ArraySetValue(DUH, zero, nlocuh*nlocu*ne);

    Int N = ne*nfe*npf*ncu;
    Kokkos::parallel_for("qJacobianLDGTraceInterior", N, KOKKOS_LAMBDA(const size_t idx) {
        Int c = idx % ncu;
        Int t = idx / ncu;
        Int a = t % npf;
        t = t / npf;
        Int lf = t % nfe;
        Int e = t / nfe;
        Int eg = e + e1;
        Int i = a + npf*lf;

        Int face = mesh.elemcon[i + ndf*eg] / npf;
        Int isboundary = (mesh.f2e[4*face + 2] < 0) ? 1 : 0;
        if (isboundary == 0) {
            Int row = i + ndf*c;
            Int col = mesh.perm[i] + npe*c;
            DUH[row + nlocuh*col + nlocuh*nlocu*e] = 0.5;
        }
    });

    for (Int e = e1; e < e2; e++) {
        for (Int lf = 0; lf < nfe; lf++) {
            Int i = npf*lf;
            Int face = mesh.elemcon[i + ndf*e] / npf;
            if (mesh.f2e[4*face + 2] < 0) {
                dstype *DUHf = &DUH[(npf*lf) + nlocuh*nlocu*(e - e1)];
                qJacobianLDGBoundaryTraceBlock(DUHf, sol, res, app, driver_abi,
                        master, mesh, tmp, common, e, lf, face, backend);
            }
        }
    }
}

inline void qJacobianLDGAssemble(dstype *DQDU, const dstype *DUH, resstruct &res,
        commonstruct &common, Int e1, Int e2, Int backend)
{
    if ((DQDU == nullptr) || (DUH == nullptr))
        return;

    Int npe = common.npe;
    Int npf = common.npf;
    Int nfe = common.nfe;
    Int ncu = common.ncu;
    Int ncq = common.ncq;
    Int nd = common.nd;
    Int ne = e2 - e1;
    Int ndf = npf*nfe;
    Int nlocu = npe*ncu;
    Int nlocuh = ndf*ncu;
    Int nlocq = npe*ncq;
    Int npenf = npe*ndf;
    Int net = common.ne;

    Int N = ne*nd*ncu*npe*nlocu;
    Kokkos::parallel_for("qJacobianLDGAssemble", N, KOKKOS_LAMBDA(const size_t idx) {
        Int colu = idx % nlocu;
        Int t = idx / nlocu;
        Int a = t % npe;
        t = t / npe;
        Int c = t % ncu;
        t = t / ncu;
        Int d = t % nd;
        Int e = t / nd;
        Int eg = e + e1;

        Int colnode = colu % npe;
        Int colcomp = colu / npe;
        Int row = a + npe*(c + ncu*d);

        dstype value = zero;
        if (colcomp == c) {
            value = res.C[a + npe*colnode + npe*npe*(eg + net*d)];
        }

        const dstype *Ed = &res.E[npenf*(eg + net*d)];
        for (Int b = 0; b < ndf; b++) {
            Int duhrow = b + ndf*c;
            value -= Ed[a + npe*b] *
                     DUH[duhrow + nlocuh*colu + nlocuh*nlocu*e];
        }

        DQDU[row + nlocq*colu + nlocq*nlocu*e] = value;
    });
}

inline void qJacobianLDG(dstype *DQDU, solstruct &sol, resstruct &res,
        appstruct &app, ExasimDriverABI& driver_abi, masterstruct &master,
        meshstruct &mesh, tempstruct &tmp, commonstruct &common,
        cublasHandle_t handle, Int e1, Int e2, Int backend)
{
    (void)MinvFull;
    (void)CFull;
    (void)EFull;

    qJacobianLDGSetupOperators(sol, res, app, master, mesh, tmp, common, handle, backend);

    dstype *DUHwork = DUH;
    if ((DQDU != nullptr) && (DUHwork == nullptr)) {
        Int npe = common.npe;
        Int npf = common.npf;
        Int nfe = common.nfe;
        Int ncu = common.ncu;
        Int ndf = npf*nfe;
        Int nlocu = npe*ncu;
        Int nlocuh = ndf*ncu;
        Int ne = e2 - e1;
        TemplateMalloc(&DUHwork, nlocuh*nlocu*ne, backend);
    }

    qJacobianLDGTrace(DUHwork, sol, res, app, driver_abi, master, mesh, tmp, common, e1, e2, backend);
    qJacobianLDGAssemble(DQDU, DUHwork, res, common, e1, e2, backend);

    if ((DQDU != nullptr) && (DUH == nullptr) && (DUHwork != nullptr)) {
        TemplateFree(DUHwork, backend);
    }
}

inline void qJacobianLDG(dstype *DQDU, solstruct &sol, resstruct &res,
        appstruct &app, ExasimDriverABI& driver_abi, masterstruct &master,
        meshstruct &mesh, tempstruct &tmp, commonstruct &common,
        cublasHandle_t handle, Int backend)
{
    qJacobianLDG(DQDU, nullptr, nullptr, nullptr, nullptr, sol, res, app,
            driver_abi, master, mesh, tmp, common, handle, 0, common.ne, backend);
}

#endif
