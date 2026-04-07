#ifndef __POINTINTERPOLATION_H__
#define __POINTINTERPOLATION_H__

#include "../Common/common.h"

#include <vector>

// The nodal shape-evaluation routine is implemented in backend/Preprocessing/makemaster*.cpp.
// This declaration allows the interpolation utilities to reuse the same basis definition.
void mkshape(std::vector<double> &shap, std::vector<double> &plocal, std::vector<double> &pts,
             int npoints, int elemtype, int porder, int nd, int numNodes);

static inline void EvaluateVolumeShapeAndGradientsAtReferencePoint(
    dstype* N,
    dstype* dN,
    const dstype* xi,
    const dstype* xpe,
    Int nd,
    Int npe,
    Int elemtype,
    Int porder)
{
    std::vector<double> plocal(static_cast<size_t>(npe * nd));
    std::vector<double> pts(static_cast<size_t>(nd));
    std::vector<double> shap(static_cast<size_t>(npe * (nd + 1)));

    for (Int i = 0; i < npe * nd; ++i) plocal[i] = static_cast<double>(xpe[i]);
    for (Int d = 0; d < nd; ++d) pts[d] = static_cast<double>(xi[d]);

    mkshape(shap, plocal, pts, 1, static_cast<int>(elemtype), static_cast<int>(porder),
            static_cast<int>(nd), static_cast<int>(npe));

    for (Int a = 0; a < npe; ++a) N[a] = static_cast<dstype>(shap[a]);
    for (Int k = 0; k < nd; ++k)
        for (Int a = 0; a < npe; ++a)
            dN[k * npe + a] = static_cast<dstype>(shap[a + npe * (k + 1)]);
}

static inline void EvaluateVolumeShapeAtReferencePoint(
    dstype* N,
    const dstype* xi,
    const dstype* xpe,
    Int nd,
    Int npe,
    Int elemtype,
    Int porder,
    Int npoints)
{
    if (N == nullptr || xi == nullptr || xpe == nullptr || npoints <= 0) return;

    std::vector<double> plocal(static_cast<size_t>(npe * nd));
    std::vector<double> pts(static_cast<size_t>(nd * npoints));
    std::vector<double> shap(static_cast<size_t>(npe * npoints * (nd + 1)));

    for (Int i = 0; i < npe * nd; ++i) plocal[i] = static_cast<double>(xpe[i]);
    for (Int i = 0; i < nd * npoints; ++i) pts[i] = static_cast<double>(xi[i]);

    mkshape(shap, plocal, pts, static_cast<int>(npoints), static_cast<int>(elemtype),
            static_cast<int>(porder), static_cast<int>(nd), static_cast<int>(npe));

    for (Int ip = 0; ip < npoints; ++ip)
        for (Int a = 0; a < npe; ++a)
            N[a + npe * ip] = static_cast<dstype>(shap[a + npe * ip]);
}

static inline void EvaluateVolumeShapeAtReferencePoint(
    dstype* N,
    const dstype* xi,
    const dstype* xpe,
    Int nd,
    Int npe,
    Int elemtype,
    Int porder)
{
    EvaluateVolumeShapeAtReferencePoint(
        N, xi, xpe, nd, npe, elemtype, porder, 1);
}

static inline void InterpolateFieldAtReferencePoint(
    dstype* u,
    const dstype* Uelem,
    const dstype* N,
    Int npe,
    Int nc)
{
    for (Int c = 0; c < nc; ++c) {
        dstype val = 0.0;
        for (Int a = 0; a < npe; ++a)
            val += N[a] * Uelem[a + npe * c];
        u[c] = val;
    }
}

static inline void InterpolateFieldBatch(
    dstype* u,
    const dstype* U,
    const Int* elems,
    const dstype* shap,
    Int npe,
    Int nc,
    Int nelems,
    Int npoints)
{
    for (Int ip = 0; ip < npoints; ++ip) {
        const Int e = elems[ip];
        const dstype* Uelem = &U[npe * nc * e];
        const dstype* N = &shap[npe * ip];
        InterpolateFieldAtReferencePoint(&u[nc * ip], Uelem, N, npe, nc);
    }
}

#endif
