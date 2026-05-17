/*
  ujacobianldg.cpp

  This file assembles the element-local LDG Jacobian block

      dRtilde_elem_K/dU_K
    + dRtilde_elem_K/dQ_K * dQ_K/dU_K
    + dRtilde_elem_K/dW_K * dW_K/dU_K

  together with the matching face contribution

      dRtilde_face_Kf/dU_K
    + dRtilde_face_Kf/dQ_K * dQ_K/dU_K
    + dRtilde_face_Kf/dW_K * dW_K/dU_K

  for each element K.

  The implementation reuses the same Gauss-point derivative drivers and metric
  operators used by the existing backend assembly paths, but keeps the direct
  U, Q, and W contributions separate so they can be composed with
  qJacobianLDG(...) and wJacobianLDG(...).

  This file does not modify the primal solver state. The face contribution is
  assembled from direct derivative face-driver blocks and then composed with
  qJacobianLDG(...) and wJacobianLDG(...).
*/
#ifndef __UJACOBIANLDG
#define __UJACOBIANLDG

inline void qJacobianLDG(dstype *DQDU, solstruct &sol, resstruct &res,
        appstruct &app, ExasimDriverABI& driver_abi, masterstruct &master,
        meshstruct &mesh, tempstruct &tmp, commonstruct &common,
        cublasHandle_t handle, Int e1, Int e2, Int backend);

inline void wJacobianLDG(dstype* DWDU, solstruct& sol, resstruct& res,
        appstruct& app, ExasimDriverABI& driver_abi, masterstruct& master,
        meshstruct& mesh, tempstruct& tmp, commonstruct& common,
        cublasHandle_t handle, Int e1, Int e2, Int backend);

// Future LDG face-derivative drivers. These are assumed to return only the
// direct derivatives needed for
//
//   dRface/dUface + dRface/dQface * dQ/dU + dRface/dWface * dW/dU
//
// and therefore exclude the uhat chain.
inline void FhatDriver(dstype* f, dstype* f_udg1, dstype* f_udg2,
        dstype* f_wdg1, dstype* f_wdg2, const dstype* xg, const dstype* ug1,
        const dstype* ug2, const dstype* og1, const dstype* og2,
        const dstype* wg1, const dstype* wg2, const dstype* uhg,
        const dstype* nl, ExasimDriverABI& driver_abi, meshstruct& mesh,
        masterstruct& master, appstruct& app, solstruct& sol,
        tempstruct& tmp, commonstruct& common, Int ngf, Int f1, Int f2,
        Int backend);

inline void FbouDriver(dstype* f, dstype* f_udg, dstype* f_wdg,
        const dstype* xg, const dstype* udg, const dstype* odg,
        const dstype* wdg, const dstype* uhg, const dstype* nl,
        ExasimDriverABI& driver_abi, meshstruct& mesh, masterstruct& master,
        appstruct& app, solstruct& sol, tempstruct& tmp,
        commonstruct& common, Int ngf, Int f1, Int f2, Int ib, Int backend);

inline void uJacobianLDGReorderElemMatrix(dstype* A, const dstype* Atmp,
        Int npe, Int ncu, Int ncol, Int ne)
{
    if ((A == nullptr) || (Atmp == nullptr) || (ncol <= 0))
        return;

    Int nlocr = npe*ncu;
    Int nlocc = npe*ncol;
    Int M = npe*npe;
    Int L = M*ne;
    Int K = L*ncu;
    Int N = nlocr*nlocc*ne;

    Kokkos::parallel_for("uJacobianLDGReorderElemMatrix", N, KOKKOS_LAMBDA(const size_t idx) {
        Int i = idx % npe;
        Int t = idx / npe;
        Int m = t % ncu;
        t = t / ncu;
        Int j = t % npe;
        t = t / npe;
        Int n = t % ncol;
        Int e = t / ncol;

        A[idx] = Atmp[i + npe*j + M*e + L*m + K*n];
    });
}

inline void uJacobianLDGDirect(dstype* AU, dstype* AQ, dstype* AW,
        solstruct& sol, resstruct& res, appstruct& app, ExasimDriverABI& driver_abi,
        masterstruct& master, meshstruct& mesh, tempstruct& tmp,
        commonstruct& common, cublasHandle_t handle, Int e1, Int e2, Int backend)
{
    (void)res;
    (void)tmp;

    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncq = common.ncq;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int ncs = common.ncs;
    Int ncw = common.ncw;
    Int nd = common.nd;
    Int npe = common.npe;
    Int nge = common.nge;

    Int ne = e2-e1;
    Int nn = npe*ne;
    Int nga = nge*ne;

    Int nlocu = npe*ncu;
    Int nlocq = npe*ncq;
    Int nlocw = npe*ncw;
    Int nrawU = npe*npe*ne*ncu*ncu;
    Int nrawQ = npe*npe*ne*ncu*ncq;
    Int nrawW = npe*npe*ne*ncu*ncw;

    Int nm = nge*e1*(ncx+nd*nd+1);
    dstype *xg = &sol.elemg[nm];
    dstype *Xx = &sol.elemg[nm + nga*ncx];
    dstype *jac = &sol.elemg[nm + nga*(ncx + nd*nd)];
    dstype *og = &sol.odgg[nge*nco*e1];

    dstype *uqg = nullptr;
    dstype *wg = nullptr;
    dstype *sg = nullptr;
    dstype *fg = nullptr;
    dstype *sg_uq = nullptr;
    dstype *fg_uq = nullptr;
    dstype *sg_w = nullptr;
    dstype *fg_w = nullptr;
    dstype *td = nullptr;
    dstype *rg = nullptr;
    dstype *rawU = nullptr;
    dstype *rawQ = nullptr;
    dstype *rawW = nullptr;

    TemplateMalloc(&uqg, nga*nc, backend);
    if (ncw > 0) TemplateMalloc(&wg, nga*ncw, backend);
    TemplateMalloc(&sg, nga*ncu, backend);
    TemplateMalloc(&fg, nga*ncu*nd, backend);
    TemplateMalloc(&sg_uq, nga*ncu*nc, backend);
    TemplateMalloc(&fg_uq, nga*ncu*nd*nc, backend);
    if (ncw > 0) {
        TemplateMalloc(&sg_w, nga*ncu*ncw, backend);
        TemplateMalloc(&fg_w, nga*ncu*nd*ncw, backend);
    }
    if (common.tdep) TemplateMalloc(&td, nga*ncu, backend);

    Int maxcol = ncu;
    if (ncq > maxcol) maxcol = ncq;
    if (ncw > maxcol) maxcol = ncw;
    if (maxcol > 0) TemplateMalloc(&rg, nge*(nd+1)*ne*ncu*maxcol, backend);

    if (AU != nullptr) TemplateMalloc(&rawU, nrawU, backend);
    if ((AQ != nullptr) && (ncq > 0)) TemplateMalloc(&rawQ, nrawQ, backend);
    if ((AW != nullptr) && (ncw > 0)) TemplateMalloc(&rawW, nrawW, backend);

    GetArrayAtIndex(tmp.tempn, sol.udg, &mesh.eindudg1[npe*nc*e1], nn*nc);
    Node2Gauss(handle, uqg, tmp.tempn, master.shapegt, nge, npe, ne*nc, backend);

    if (ncw > 0) {
        GetElemNodes(tmp.tempn, sol.wdg, npe, ncw, 0, ncw, e1, e2);
        Node2Gauss(handle, wg, tmp.tempn, master.shapegt, nge, npe, ne*ncw, backend);
    }

    ArraySetValue(sg, zero, nga*ncu);
    ArraySetValue(sg_uq, zero, nga*ncu*nc);
    if (ncw > 0) ArraySetValue(sg_w, zero, nga*ncu*ncw);
    SourceDriver(sg, sg_uq, sg_w, xg, uqg, og, wg, driver_abi, mesh, master,
            app, sol, tmp, common, nge, e1, e2, backend);

    if (common.tdep) {
        ArrayAXPBY(td, &sol.sdgg[nge*ncs*e1], uqg, one, -common.dtfactor, nga*ncu);

        if (common.tdfunc==1)
            TdfuncDriver(fg, xg, uqg, og, wg, driver_abi, mesh, master, app, sol, tmp, common, nge, e1, e2, backend);
        else
            ArraySetValue(fg, one, nga*ncu);

        ArrayAXY(td, td, fg, one, nga*ncu);
        ArrayAXPBY(sg, sg, td, one, one, nga*ncu);
        ApplyDtcoef(sg_uq, fg, -common.dtfactor, nga, ncu);
    }

    ArraySetValue(fg, zero, nga*ncu*nd);
    ArraySetValue(fg_uq, zero, nga*ncu*nd*nc);
    if (ncw > 0) ArraySetValue(fg_w, zero, nga*ncu*nd*ncw);
    FluxDriver(fg, fg_uq, fg_w, xg, uqg, og, wg, driver_abi, mesh, master,
            app, sol, tmp, common, nge, e1, e2, backend);

    if (AU != nullptr) {
        ApplyXxJac(rg, sg_uq, fg_uq, Xx, jac, nge, nd, ncu, ncu, ne);
        Gauss2Node(handle, rawU, rg, master.shapegwdotshapeg, nge*(nd+1), npe*npe, ncu*ncu*ne, backend);
        uJacobianLDGReorderElemMatrix(AU, rawU, npe, ncu, ncu, ne);
    }

    if ((AQ != nullptr) && (ncq > 0)) {
        ApplyXxJac(rg, &sg_uq[nga*ncu*ncu], &fg_uq[nga*ncu*nd*ncu], Xx, jac, nge, nd, ncu, ncq, ne);
        Gauss2Node(handle, rawQ, rg, master.shapegwdotshapeg, nge*(nd+1), npe*npe, ncu*ncq*ne, backend);
        uJacobianLDGReorderElemMatrix(AQ, rawQ, npe, ncu, ncq, ne);
    }

    if ((AW != nullptr) && (ncw > 0)) {
        ApplyXxJac(rg, sg_w, fg_w, Xx, jac, nge, nd, ncu, ncw, ne);
        Gauss2Node(handle, rawW, rg, master.shapegwdotshapeg, nge*(nd+1), npe*npe, ncu*ncw*ne, backend);
        uJacobianLDGReorderElemMatrix(AW, rawW, npe, ncu, ncw, ne);
    }

    if (rawW != nullptr) TemplateFree(rawW, backend);
    if (rawQ != nullptr) TemplateFree(rawQ, backend);
    if (rawU != nullptr) TemplateFree(rawU, backend);
    if (rg != nullptr) TemplateFree(rg, backend);
    if (td != nullptr) TemplateFree(td, backend);
    if (fg_w != nullptr) TemplateFree(fg_w, backend);
    if (sg_w != nullptr) TemplateFree(sg_w, backend);
    TemplateFree(fg_uq, backend);
    TemplateFree(sg_uq, backend);
    TemplateFree(fg, backend);
    TemplateFree(sg, backend);
    if (wg != nullptr) TemplateFree(wg, backend);
    TemplateFree(uqg, backend);
}

inline void uJacobianLDGCompose(dstype* JU, const dstype* AU, const dstype* AQ,
        const dstype* AW, const dstype* DQDU, const dstype* DWDU,
        commonstruct& common, cublasHandle_t handle, Int e1, Int e2, Int backend)
{
    if ((JU == nullptr) || (AU == nullptr))
        return;

    Int npe = common.npe;
    Int ncu = common.ncu;
    Int ncq = common.ncq;
    Int ncw = common.ncw;
    Int ne = e2 - e1;
    Int nlocr = npe*ncu;
    Int nlocu = npe*ncu;
    Int nlocq = npe*ncq;
    Int nlocw = npe*ncw;

    ArrayCopy(JU, AU, nlocr*nlocu*ne);

    if ((ncq > 0) && (AQ != nullptr) && (DQDU != nullptr)) {
        PGEMNMStridedBached(handle, nlocr, nlocu, nlocq, one,
                AQ, nlocr, DQDU, nlocq, one, JU, nlocr, ne, backend);
    }

    if ((ncw > 0) && (AW != nullptr) && (DWDU != nullptr)) {
        PGEMNMStridedBached(handle, nlocr, nlocu, nlocw, one,
                AW, nlocr, DWDU, nlocw, one, JU, nlocr, ne, backend);
    }
}

inline void uJacobianLDGAssembleFaceSide(dstype* A, const dstype* Atmp,
        meshstruct& mesh, commonstruct& common, Int e1, Int e2, Int f1, Int f2,
        Int ncol, Int side, dstype sign)
{
    if ((A == nullptr) || (Atmp == nullptr) || (ncol <= 0))
        return;

    Int npe = common.npe;
    Int npf = common.npf;
    Int ncu = common.ncu;
    Int nf = f2 - f1;
    Int ne = e2 - e1;
    Int nlocu = npe*ncu;
    Int nlocc = npe*ncol;
    Int N = npf*npf*nf*ncu*ncol;

    Kokkos::parallel_for("uJacobianLDGAssembleFaceSide", N, KOKKOS_LAMBDA(const size_t idx) {
        Int p = idx % npf;
        Int t = idx / npf;
        Int q = t % npf;
        t = t / npf;
        Int f = t % nf;
        t = t / nf;
        Int r = t % ncu;
        Int c = t / ncu;

        Int face = f1 + f;
        Int sideoffset = side - 1;

        Int mr = npf*face + p;
        Int mc = npf*face + q;
        Int kr = mesh.facecon[2*mr + sideoffset];
        Int kc = mesh.facecon[2*mc + sideoffset];

        Int rownode = kr % npe;
        Int rowelem = (kr - rownode) / npe;
        Int colnode = kc % npe;
        Int colelem = (kc - colnode) / npe;

        if ((rowelem >= e1) && (rowelem < e2) && (rowelem == colelem)) {
            Int elocal = rowelem - e1;
            Int row = rownode + npe*r;
            Int col = colnode + npe*c;
            Int M = npf*npf;
            dstype value = sign * Atmp[p + npf*q + M*(f + nf*(r + ncu*c))];
            Kokkos::atomic_add(&A[row + nlocu*col + nlocu*nlocc*elocal], value);
        }
    });
}

inline void uJacobianLDGFace(dstype* JU, dstype* JQ, dstype* JW,
        solstruct& sol, resstruct& res, appstruct& app, ExasimDriverABI& driver_abi,
        masterstruct& master, meshstruct& mesh, tempstruct& tmp,
        commonstruct& common, cublasHandle_t handle, Int e1, Int e2, Int backend)
{
    (void)res;
    (void)tmp;

    Int npe = common.npe;
    Int npf = common.npf;
    Int ncu = common.ncu;
    Int nc = common.nc;
    Int ncq = common.ncq;
    Int nco = common.nco;
    Int ncw = common.ncw;
    Int nlocu = npe*ncu;
    Int nlocq = npe*ncq;
    Int nlocw = npe*ncw;
    Int ne = e2 - e1;

    if (JU != nullptr) ArraySetValue(JU, zero, nlocu*nlocu*ne);
    if ((JQ != nullptr) && (ncq > 0)) ArraySetValue(JQ, zero, nlocu*nlocq*ne);
    if ((JW != nullptr) && (ncw > 0)) ArraySetValue(JW, zero, nlocu*nlocw*ne);

    for (Int j = 0; j < common.nbf; ++j) {
        Int f1 = common.fblks[3*j] - 1;
        Int f2 = common.fblks[3*j + 1];
        Int ib = common.fblks[3*j + 2];
        Int nf = f2 - f1;
        Int nn = npf*nf;
        Int ngf = common.ngf;
        Int nga = ngf*nf;
        Int ncx = common.ncx;
        Int nd = common.nd;

        Int nodalCols = (ib == 0) ? (ncu + 2*nc + 2*ncw) : (ncu + nc + ncw);
        dstype* fn = nullptr;
        dstype* fg = nullptr;
        dstype* fh = nullptr;
        dstype* fh_udg1 = nullptr;
        dstype* fh_udg2 = nullptr;
        dstype* fh_wdg1 = nullptr;
        dstype* fh_wdg2 = nullptr;
        dstype* Utmp1 = nullptr;
        dstype* Utmp2 = nullptr;
        dstype* Qtmp1 = nullptr;
        dstype* Qtmp2 = nullptr;
        dstype* Wtmp1 = nullptr;
        dstype* Wtmp2 = nullptr;

        TemplateMalloc(&fn, nn*nodalCols, backend);
        TemplateMalloc(&fg, nga*nodalCols, backend);
        TemplateMalloc(&fh, nga*ncu, backend);
        TemplateMalloc(&fh_udg1, nga*ncu*nc, backend);
        if ((JU != nullptr) || ((JQ != nullptr) && (ncq > 0)))
            TemplateMalloc(&Utmp1, npf*npf*nf*ncu*ncu, backend);
        if ((JQ != nullptr) && (ncq > 0))
            TemplateMalloc(&Qtmp1, npf*npf*nf*ncu*ncq, backend);
        if ((JW != nullptr) && (ncw > 0))
            TemplateMalloc(&Wtmp1, npf*npf*nf*ncu*ncw, backend);
        if (ncw > 0)
            TemplateMalloc(&fh_wdg1, nga*ncu*ncw, backend);

        GetElemNodes(fn, sol.uh, npf, ncu, 0, ncu, f1, f2);
        GetArrayAtIndex(&fn[nn*ncu], sol.udg, &mesh.findudg1[npf*nc*f1], nn*nc);
        if (ncw > 0)
            GetFaceNodes(&fn[nn*(ncu+nc)], sol.wdg, mesh.facecon, npf, ncw, npe, ncw, f1, f2, 1);

        if (ib == 0) {
            TemplateMalloc(&fh_udg2, nga*ncu*nc, backend);
            if ((JU != nullptr) || ((JQ != nullptr) && (ncq > 0)))
                TemplateMalloc(&Utmp2, npf*npf*nf*ncu*ncu, backend);
            if ((JQ != nullptr) && (ncq > 0))
                TemplateMalloc(&Qtmp2, npf*npf*nf*ncu*ncq, backend);
            if ((JW != nullptr) && (ncw > 0))
                TemplateMalloc(&Wtmp2, npf*npf*nf*ncu*ncw, backend);
            if (ncw > 0)
                TemplateMalloc(&fh_wdg2, nga*ncu*ncw, backend);

            GetArrayAtIndex(&fn[nn*(ncu+nc+ncw)], sol.udg, &mesh.findudg2[npf*nc*f1], nn*nc);
            if (ncw > 0)
                GetFaceNodes(&fn[nn*(ncu+nc+ncw+nc)], sol.wdg, mesh.facecon, npf, ncw, npe, ncw, f1, f2, 2);

            Node2Gauss(handle, fg, fn, master.shapfgt, ngf, npf, nf*nodalCols, backend);

            dstype* uhg = &fg[0];
            dstype* ug1 = &fg[nga*ncu];
            dstype* wg1 = &fg[nga*(ncu+nc)];
            dstype* ug2 = &fg[nga*(ncu+nc+ncw)];
            dstype* wg2 = &fg[nga*(ncu+2*nc+ncw)];
            const dstype* xg = &sol.faceg[ngf*f1*(ncx+nd+1)];
            const dstype* nlg = &sol.faceg[ngf*f1*(ncx+nd+1) + nga*ncx];
            const dstype* jac = &sol.faceg[ngf*f1*(ncx+nd+1) + nga*(ncx+nd)];

            ArraySetValue(fh, zero, nga*ncu);
            ArraySetValue(fh_udg1, zero, nga*ncu*nc);
            ArraySetValue(fh_udg2, zero, nga*ncu*nc);
            if (ncw > 0) {
                ArraySetValue(fh_wdg1, zero, nga*ncu*ncw);
                ArraySetValue(fh_wdg2, zero, nga*ncu*ncw);
            }

            FhatDriver(fh, fh_udg1, fh_udg2, fh_wdg1, fh_wdg2, xg, ug1, ug2,
                    &sol.og1[ngf*nco*f1], &sol.og2[ngf*nco*f1], wg1, wg2, uhg,
                    nlg, driver_abi, mesh, master, app, sol, tmp, common,
                    ngf, f1, f2, backend);

            columnwiseMultiply(fh_udg1, fh_udg1, jac, nga, ncu*nc);
            columnwiseMultiply(fh_udg2, fh_udg2, jac, nga, ncu*nc);
            if (ncw > 0) {
                columnwiseMultiply(fh_wdg1, fh_wdg1, jac, nga, ncu*ncw);
                columnwiseMultiply(fh_wdg2, fh_wdg2, jac, nga, ncu*ncw);
            }

            if (JU != nullptr) {
                Gauss2Node(handle, Utmp1, fh_udg1, master.shapfgwdotshapfg, ngf, npf*npf, nf*ncu*ncu, backend);
                Gauss2Node(handle, Utmp2, fh_udg2, master.shapfgwdotshapfg, ngf, npf*npf, nf*ncu*ncu, backend);
                uJacobianLDGAssembleFaceSide(JU, Utmp1, mesh, common, e1, e2, f1, f2, ncu, 1, minusone);
                uJacobianLDGAssembleFaceSide(JU, Utmp2, mesh, common, e1, e2, f1, f2, ncu, 2, one);
            }

            if ((JQ != nullptr) && (ncq > 0)) {
                Gauss2Node(handle, Qtmp1, &fh_udg1[nga*ncu*ncu], master.shapfgwdotshapfg, ngf, npf*npf, nf*ncu*ncq, backend);
                Gauss2Node(handle, Qtmp2, &fh_udg2[nga*ncu*ncu], master.shapfgwdotshapfg, ngf, npf*npf, nf*ncu*ncq, backend);
                uJacobianLDGAssembleFaceSide(JQ, Qtmp1, mesh, common, e1, e2, f1, f2, ncq, 1, minusone);
                uJacobianLDGAssembleFaceSide(JQ, Qtmp2, mesh, common, e1, e2, f1, f2, ncq, 2, one);
            }

            if ((JW != nullptr) && (ncw > 0)) {
                Gauss2Node(handle, Wtmp1, fh_wdg1, master.shapfgwdotshapfg, ngf, npf*npf, nf*ncu*ncw, backend);
                Gauss2Node(handle, Wtmp2, fh_wdg2, master.shapfgwdotshapfg, ngf, npf*npf, nf*ncu*ncw, backend);
                uJacobianLDGAssembleFaceSide(JW, Wtmp1, mesh, common, e1, e2, f1, f2, ncw, 1, minusone);
                uJacobianLDGAssembleFaceSide(JW, Wtmp2, mesh, common, e1, e2, f1, f2, ncw, 2, one);
            }
        }
        else {
            Node2Gauss(handle, fg, fn, master.shapfgt, ngf, npf, nf*nodalCols, backend);

            dstype* uhg = &fg[0];
            dstype* ug1 = &fg[nga*ncu];
            dstype* wg1 = &fg[nga*(ncu+nc)];
            const dstype* xg = &sol.faceg[ngf*f1*(ncx+nd+1)];
            const dstype* nlg = &sol.faceg[ngf*f1*(ncx+nd+1) + nga*ncx];
            const dstype* jac = &sol.faceg[ngf*f1*(ncx+nd+1) + nga*(ncx+nd)];

            ArraySetValue(fh, zero, nga*ncu);
            ArraySetValue(fh_udg1, zero, nga*ncu*nc);
            if (ncw > 0)
                ArraySetValue(fh_wdg1, zero, nga*ncu*ncw);

            FbouDriver(fh, fh_udg1, fh_wdg1, xg, ug1, &sol.og1[ngf*nco*f1], wg1,
                    uhg, nlg, driver_abi, mesh, master, app, sol, tmp, common,
                    ngf, f1, f2, ib, backend);

            columnwiseMultiply(fh_udg1, fh_udg1, jac, nga, ncu*nc);
            if (ncw > 0)
                columnwiseMultiply(fh_wdg1, fh_wdg1, jac, nga, ncu*ncw);

            if (JU != nullptr) {
                Gauss2Node(handle, Utmp1, fh_udg1, master.shapfgwdotshapfg, ngf, npf*npf, nf*ncu*ncu, backend);
                uJacobianLDGAssembleFaceSide(JU, Utmp1, mesh, common, e1, e2, f1, f2, ncu, 1, minusone);
            }

            if ((JQ != nullptr) && (ncq > 0)) {
                Gauss2Node(handle, Qtmp1, &fh_udg1[nga*ncu*ncu], master.shapfgwdotshapfg, ngf, npf*npf, nf*ncu*ncq, backend);
                uJacobianLDGAssembleFaceSide(JQ, Qtmp1, mesh, common, e1, e2, f1, f2, ncq, 1, minusone);
            }

            if ((JW != nullptr) && (ncw > 0)) {
                Gauss2Node(handle, Wtmp1, fh_wdg1, master.shapfgwdotshapfg, ngf, npf*npf, nf*ncu*ncw, backend);
                uJacobianLDGAssembleFaceSide(JW, Wtmp1, mesh, common, e1, e2, f1, f2, ncw, 1, minusone);
            }
        }

        if (Wtmp2 != nullptr) TemplateFree(Wtmp2, backend);
        if (Wtmp1 != nullptr) TemplateFree(Wtmp1, backend);
        if (Qtmp2 != nullptr) TemplateFree(Qtmp2, backend);
        if (Qtmp1 != nullptr) TemplateFree(Qtmp1, backend);
        if (Utmp2 != nullptr) TemplateFree(Utmp2, backend);
        if (Utmp1 != nullptr) TemplateFree(Utmp1, backend);
        if (fh_wdg2 != nullptr) TemplateFree(fh_wdg2, backend);
        if (fh_wdg1 != nullptr) TemplateFree(fh_wdg1, backend);
        if (fh_udg2 != nullptr) TemplateFree(fh_udg2, backend);
        TemplateFree(fh_udg1, backend);
        TemplateFree(fh, backend);
        TemplateFree(fg, backend);
        TemplateFree(fn, backend);
    }
}

inline void uJacobianLDG(dstype* JU, dstype* AU, dstype* AQ, dstype* AW,
        dstype* DQDU, dstype* DWDU, solstruct& sol, resstruct& res,
        appstruct& app, ExasimDriverABI& driver_abi, masterstruct& master,
        meshstruct& mesh, tempstruct& tmp, commonstruct& common,
        cublasHandle_t handle, Int e1, Int e2, Int backend)
{
    Int npe = common.npe;
    Int ncu = common.ncu;
    Int ncq = common.ncq;
    Int ncw = common.ncw;
    Int ne = e2 - e1;
    Int nlocu = npe*ncu;
    Int nlocq = npe*ncq;
    Int nlocw = npe*ncw;

    dstype* AUwork = AU;
    dstype* AQwork = AQ;
    dstype* AWwork = AW;
    dstype* DQDUwork = DQDU;
    dstype* DWDUwork = DWDU;
    dstype* JUface = nullptr;
    dstype* JQface = nullptr;
    dstype* JWface = nullptr;
    dstype* JFwork = nullptr;

    if ((JU != nullptr) && (AUwork == nullptr))
        TemplateMalloc(&AUwork, nlocu*nlocu*ne, backend);
    if ((JU != nullptr) && (ncq > 0) && (AQwork == nullptr))
        TemplateMalloc(&AQwork, nlocu*nlocq*ne, backend);
    if ((JU != nullptr) && (ncw > 0) && (AWwork == nullptr))
        TemplateMalloc(&AWwork, nlocu*nlocw*ne, backend);
    if ((JU != nullptr) && (ncq > 0) && (DQDUwork == nullptr))
        TemplateMalloc(&DQDUwork, nlocq*nlocu*ne, backend);
    if ((JU != nullptr) && (ncw > 0) && (DWDUwork == nullptr))
        TemplateMalloc(&DWDUwork, nlocw*nlocu*ne, backend);

    uJacobianLDGDirect(AUwork, AQwork, AWwork, sol, res, app, driver_abi, master,
            mesh, tmp, common, handle, e1, e2, backend);

    if ((ncq > 0) && (DQDUwork != nullptr))
        qJacobianLDG(DQDUwork, sol, res, app, driver_abi, master, mesh, tmp,
                common, handle, e1, e2, backend);

    if ((ncw > 0) && (DWDUwork != nullptr))
        wJacobianLDG(DWDUwork, sol, res, app, driver_abi, master, mesh, tmp,
                common, handle, e1, e2, backend);

    uJacobianLDGCompose(JU, AUwork, AQwork, AWwork, DQDUwork, DWDUwork,
            common, handle, e1, e2, backend);

    if (JU != nullptr) {
        TemplateMalloc(&JFwork, nlocu*nlocu*ne, backend);

        TemplateMalloc(&JUface, nlocu*nlocu*ne, backend);
        if (ncq > 0)
            TemplateMalloc(&JQface, nlocu*nlocq*ne, backend);
        if (ncw > 0)
            TemplateMalloc(&JWface, nlocu*nlocw*ne, backend);

        uJacobianLDGFace(JUface, JQface, JWface, sol, res, app, driver_abi,
                master, mesh, tmp, common, handle, e1, e2, backend);

        uJacobianLDGCompose(JFwork, JUface, JQface, JWface, DQDUwork, DWDUwork,
                common, handle, e1, e2, backend);

        ArrayAXPBY(JU, JU, JFwork, one, one, nlocu*nlocu*ne);

        if (JWface != nullptr) TemplateFree(JWface, backend);
        if (JQface != nullptr) TemplateFree(JQface, backend);
        if (JUface != nullptr) TemplateFree(JUface, backend);
        TemplateFree(JFwork, backend);
    }

    if ((JU != nullptr) && (DWDU == nullptr) && (DWDUwork != nullptr))
        TemplateFree(DWDUwork, backend);
    if ((JU != nullptr) && (DQDU == nullptr) && (DQDUwork != nullptr))
        TemplateFree(DQDUwork, backend);
    if ((JU != nullptr) && (AW == nullptr) && (AWwork != nullptr))
        TemplateFree(AWwork, backend);
    if ((JU != nullptr) && (AQ == nullptr) && (AQwork != nullptr))
        TemplateFree(AQwork, backend);
    if ((JU != nullptr) && (AU == nullptr) && (AUwork != nullptr))
        TemplateFree(AUwork, backend);
}

inline void uJacobianLDG(dstype* JU, solstruct& sol, resstruct& res,
        appstruct& app, ExasimDriverABI& driver_abi, masterstruct& master,
        meshstruct& mesh, tempstruct& tmp, commonstruct& common,
        cublasHandle_t handle, Int e1, Int e2, Int backend)
{
    uJacobianLDG(JU, nullptr, nullptr, nullptr, nullptr, nullptr, sol, res, app,
            driver_abi, master, mesh, tmp, common, handle, e1, e2, backend);
}

inline void uJacobianLDG(dstype* JU, solstruct& sol, resstruct& res,
        appstruct& app, ExasimDriverABI& driver_abi, masterstruct& master,
        meshstruct& mesh, tempstruct& tmp, commonstruct& common,
        cublasHandle_t handle, Int backend)
{
    uJacobianLDG(JU, nullptr, nullptr, nullptr, nullptr, nullptr, sol, res, app,
            driver_abi, master, mesh, tmp, common, handle, 0, common.ne, backend);
}

#endif
