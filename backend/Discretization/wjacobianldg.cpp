/*
  wjacobianldg.cpp

  This file assembles the element-local LDG Jacobian block

      dW_K / dU_K

  where U_K contains the degrees of freedom of u_h only, with local size
      npe * ncu,
  and W_K contains the degrees of freedom of w_h, with local size
      npe * ncw.

  The implementation matches the active GetW(...) branches in residual.cpp:

  1. wave == 1
       w = (1/dtfactor) * (u + wsrc)

  2. wave == 0 and dae_alpha == dae_beta == 0
       Eos(x,u,o,w) = 0
       dw/du = - (dEos/dw)^{-1} (dEos/du)

  3. wave == 0 and not (dae_alpha == dae_beta == 0)
       explicit update used by GetW(...):

         if dae_steps == 0:
             w = scalar * (dae_alpha * wsrc + Sourcew(x,u,o,w))

         else:
             w = scalar * (dae_gamma * wdual + dae_alpha * wsrc + Sourcew(x,u,o,w))

       For the implemented update map, dW/dU is therefore:

         dW/dU = scalar * dSourcew/du

       evaluated at the current stored sol.wdg.

  This file intentionally does not modify global solver state. It assumes
  sol.wdg already contains the state W at which the Jacobian should be
  evaluated.
*/
#ifndef __WJACOBIANLDG
#define __WJACOBIANLDG

inline void wJacobianLDGZero(dstype* DWDU, commonstruct& common, Int e1, Int e2)
{
    if (DWDU == nullptr)
        return;

    Int npe = common.npe;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int ne = e2 - e1;
    ArraySetValue(DWDU, zero, npe*ncw*npe*ncu*ne);
}

inline void wJacobianLDGInsertPointwise(dstype* DWDU, const dstype* dwdudg,
        commonstruct& common, Int e1, Int e2, Int ncols)
{
    if ((DWDU == nullptr) || (dwdudg == nullptr))
        return;

    Int npe = common.npe;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int ne = e2 - e1;
    Int nn = npe*ne;
    Int nlocu = npe*ncu;
    Int nlocw = npe*ncw;

    ArraySetValue(DWDU, zero, nlocw*nlocu*ne);

    Int cols = (ncols < ncu) ? ncols : ncu;
    Int N = ne*npe*ncw*cols;
    Kokkos::parallel_for("wJacobianLDGInsertPointwise", N, KOKKOS_LAMBDA(const size_t idx) {
        Int uc = idx % cols;
        Int t = idx / cols;
        Int wc = t % ncw;
        t = t / ncw;
        Int a = t % npe;
        Int e = t / npe;

        Int row = a + npe*wc;
        Int col = a + npe*uc;
        Int i = a + npe*e;

        DWDU[row + nlocw*col + nlocw*nlocu*e] =
            dwdudg[i + nn*wc + nn*ncw*uc];
    });
}

inline void wJacobianLDGWave(dstype* DWDU, commonstruct& common, Int e1, Int e2)
{
    if (DWDU == nullptr)
        return;

    Int npe = common.npe;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int ne = e2 - e1;
    Int nlocu = npe*ncu;
    Int nlocw = npe*ncw;
    dstype scalar = one/common.dtfactor;

    ArraySetValue(DWDU, zero, nlocw*nlocu*ne);

    Int ncopy = (ncw < ncu) ? ncw : ncu;
    Int N = ne*npe*ncopy;
    Kokkos::parallel_for("wJacobianLDGWave", N, KOKKOS_LAMBDA(const size_t idx) {
        Int c = idx % ncopy;
        Int t = idx / ncopy;
        Int a = t % npe;
        Int e = t / npe;

        Int row = a + npe*c;
        Int col = a + npe*c;
        DWDU[row + nlocw*col + nlocw*nlocu*e] = scalar;
    });
}

inline void wJacobianLDGEos(dstype* DWDU, dstype* Fu, dstype* Fw,
        solstruct& sol, resstruct& res, appstruct& app, ExasimDriverABI& driver_abi,
        masterstruct& master, meshstruct& mesh, tempstruct& tmp,
        commonstruct& common, Int e1, Int e2, Int backend)
{
    (void)res;

    if (DWDU == nullptr)
        return;

    Int npe = common.npe;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nc = common.nc;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int ne = e2 - e1;
    Int nn = npe*ne;

    EosduDriver(Fu, &sol.xdg[npe*ncx*e1], &sol.udg[npe*nc*e1], &sol.odg[npe*nco*e1],
            &sol.wdg[npe*ncw*e1], driver_abi, mesh, master, app, sol, tmp, common,
            npe, e1, e2, backend);

    EosdwDriver(Fw, &sol.xdg[npe*ncx*e1], &sol.udg[npe*nc*e1], &sol.odg[npe*nco*e1],
            &sol.wdg[npe*ncw*e1], driver_abi, mesh, master, app, sol, tmp, common,
            npe, e1, e2, backend);

    if (ncw==1)
        SmallMatrixSolve11(Fu, Fw, nn, ncu);
    else if (ncw==2)
        SmallMatrixSolve22(Fu, Fw, nn, ncu);
    else if (ncw==3)
        SmallMatrixSolve33(Fu, Fw, nn, ncu);
    else if (ncw==4)
        SmallMatrixSolve44(Fu, Fw, nn, ncu);
    else if (ncw==5)
        SmallMatrixSolve55(Fu, Fw, nn, ncu);
    else
        error("wjacobianldg supports at most five w variables.");

    ArrayMultiplyScalar(Fu, minusone, nn*ncw*ncu);
    wJacobianLDGInsertPointwise(DWDU, Fu, common, e1, e2, ncu);
}

inline void wJacobianLDGSourcew(dstype* DWDU, dstype* Su, dstype* Sw,
        solstruct& sol, resstruct& res, appstruct& app, ExasimDriverABI& driver_abi,
        masterstruct& master, meshstruct& mesh, tempstruct& tmp,
        commonstruct& common, Int e1, Int e2, Int backend)
{
    (void)res;
    (void)Sw;

    if (DWDU == nullptr)
        return;

    Int npe = common.npe;
    Int nc = common.nc;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int ne = e2 - e1;
    Int nn = npe*ne;

    SourcewDriver(&tmp.tempn[0], Su, Sw, &sol.xdg[npe*ncx*e1], &sol.udg[npe*nc*e1],
            &sol.odg[npe*nco*e1], &sol.wdg[npe*ncw*e1], driver_abi, mesh, master,
            app, sol, tmp, common, npe, e1, e2, backend);

    dstype scalar = zero;
    if (common.dae_steps==0)
        scalar = one/(common.dae_alpha*common.dtfactor + common.dae_beta);
    else
        scalar = one/(common.dae_alpha*common.dtfactor + common.dae_beta + common.dae_gamma);

    ArrayMultiplyScalar(Su, scalar, nn*ncw*nc);
    wJacobianLDGInsertPointwise(DWDU, Su, common, e1, e2, common.ncu);
}

inline void wJacobianLDG(dstype* DWDU, solstruct& sol, resstruct& res,
        appstruct& app, ExasimDriverABI& driver_abi, masterstruct& master,
        meshstruct& mesh, tempstruct& tmp, commonstruct& common,
        cublasHandle_t handle, Int e1, Int e2, Int backend)
{
    (void)handle;

    if (DWDU == nullptr)
        return;

    if ((common.subproblem != 0) || (common.ncw <= 0)) {
        wJacobianLDGZero(DWDU, common, e1, e2);
        return;
    }

    if (common.wave == 1) {
        wJacobianLDGWave(DWDU, common, e1, e2);
        return;
    }

    Int npe = common.npe;
    Int ncu = common.ncu;
    Int nc = common.nc;
    Int ncw = common.ncw;
    Int ne = e2 - e1;
    Int nn = npe*ne;

    if ((fabs(common.dae_alpha) < 1e-10) && (fabs(common.dae_beta) < 1e-10)) {
        dstype* Fu = nullptr;
        dstype* Fw = nullptr;
        TemplateMalloc(&Fu, nn*ncw*ncu, backend);
        TemplateMalloc(&Fw, nn*ncw*ncw, backend);

        wJacobianLDGEos(DWDU, Fu, Fw, sol, res, app, driver_abi, master, mesh,
                tmp, common, e1, e2, backend);

        TemplateFree(Fu, backend);
        TemplateFree(Fw, backend);
    }
    else {
        dstype* Su = nullptr;
        dstype* Sw = nullptr;
        TemplateMalloc(&Su, nn*ncw*nc, backend);
        TemplateMalloc(&Sw, nn*ncw*ncw, backend);

        wJacobianLDGSourcew(DWDU, Su, Sw, sol, res, app, driver_abi, master, mesh,
                tmp, common, e1, e2, backend);

        TemplateFree(Su, backend);
        TemplateFree(Sw, backend);
    }
}

inline void wJacobianLDG(dstype* DWDU, solstruct& sol, resstruct& res,
        appstruct& app, ExasimDriverABI& driver_abi, masterstruct& master,
        meshstruct& mesh, tempstruct& tmp, commonstruct& common,
        cublasHandle_t handle, Int backend)
{
    wJacobianLDG(DWDU, sol, res, app, driver_abi, master, mesh, tmp, common,
            handle, 0, common.ne, backend);
}

#endif
