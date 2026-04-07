#include "wallmodelsampling.h"

#include <algorithm>
#include <cmath>
#include <limits>

// Defined in makemasterexasim.cpp and already part of backend preprocessing.
void mkshape(std::vector<double>& shap, std::vector<double>& plocal, std::vector<double>& pts,
             int npoints, int elemtype, int porder, int nd, int numNodes);

namespace exasim {
namespace wm {

namespace {

static inline dstype VecNorm(const dstype* a, Int n)
{
    dstype s = 0.0;
    for (Int i = 0; i < n; ++i) s += a[i] * a[i];
    return std::sqrt(s);
}

static inline void DefaultInitialGuess(dstype* xi, Int nd, Int elemtype)
{
    if (elemtype == 0) {
        const dstype val = 1.0 / static_cast<dstype>(nd + 1);
        for (Int d = 0; d < nd; ++d) xi[d] = val;
    } else {
        for (Int d = 0; d < nd; ++d) xi[d] = 0.5;
    }
}

static bool SolveSmallLinearSystem(dstype* x, const dstype* A, const dstype* b, Int nd)
{
    if (nd == 1) {
        if (std::abs(A[0]) <= std::numeric_limits<dstype>::epsilon()) return false;
        x[0] = b[0] / A[0];
        return true;
    }

    if (nd == 2) {
        const dstype det = A[0] * A[3] - A[1] * A[2];
        if (std::abs(det) <= std::numeric_limits<dstype>::epsilon()) return false;
        x[0] = ( A[3] * b[0] - A[1] * b[1]) / det;
        x[1] = (-A[2] * b[0] + A[0] * b[1]) / det;
        return true;
    }

    if (nd == 3) {
        const dstype det =
            A[0] * (A[4] * A[8] - A[5] * A[7]) -
            A[1] * (A[3] * A[8] - A[5] * A[6]) +
            A[2] * (A[3] * A[7] - A[4] * A[6]);
        if (std::abs(det) <= std::numeric_limits<dstype>::epsilon()) return false;

        const dstype inv[9] = {
             (A[4] * A[8] - A[5] * A[7]) / det,
            -(A[1] * A[8] - A[2] * A[7]) / det,
             (A[1] * A[5] - A[2] * A[4]) / det,
            -(A[3] * A[8] - A[5] * A[6]) / det,
             (A[0] * A[8] - A[2] * A[6]) / det,
            -(A[0] * A[5] - A[2] * A[3]) / det,
             (A[3] * A[7] - A[4] * A[6]) / det,
            -(A[0] * A[7] - A[1] * A[6]) / det,
             (A[0] * A[4] - A[1] * A[3]) / det
        };

        for (Int i = 0; i < 3; ++i) {
            x[i] = 0.0;
            for (Int j = 0; j < 3; ++j)
                x[i] += inv[i * 3 + j] * b[j];
        }
        return true;
    }

    return false;
}

static void CopyMasterNodes(std::vector<double>& plocal, const masterstruct& master, Int npe, Int nd)
{
    plocal.resize(static_cast<size_t>(npe * nd));
    for (Int i = 0; i < npe * nd; ++i) plocal[i] = static_cast<double>(master.xpe[i]);
}

} // namespace

void EvaluateVolumeShapeAndGradientsAtXi(
    dstype* N,
    dstype* dN,
    const dstype* xi,
    const masterstruct& master,
    const Int nd,
    const Int elemtype,
    const Int porder,
    const Int npe)
{
    std::vector<double> plocal;
    std::vector<double> pts(static_cast<size_t>(nd));
    std::vector<double> shap;

    CopyMasterNodes(plocal, master, npe, nd);
    for (Int d = 0; d < nd; ++d) pts[d] = static_cast<double>(xi[d]);

    mkshape(shap, plocal, pts, 1, static_cast<int>(elemtype), static_cast<int>(porder),
            static_cast<int>(nd), static_cast<int>(npe));

    for (Int a = 0; a < npe; ++a) N[a] = static_cast<dstype>(shap[a]);
    for (Int k = 0; k < nd; ++k)
        for (Int a = 0; a < npe; ++a)
            dN[k * npe + a] = static_cast<dstype>(shap[a + npe * (k + 1)]);
}

void MapReferenceToPhysical(
    dstype* x,
    dstype* J,
    const dstype* xi,
    const dstype* xdg_elem,
    const masterstruct& master,
    const Int nd,
    const Int elemtype,
    const Int porder,
    const Int npe)
{
    std::vector<dstype> N(static_cast<size_t>(npe));
    std::vector<dstype> dN(static_cast<size_t>(nd * npe));

    EvaluateVolumeShapeAndGradientsAtXi(N.data(), dN.data(), xi, master, nd, elemtype, porder, npe);

    for (Int d = 0; d < nd; ++d) x[d] = 0.0;
    for (Int i = 0; i < nd * nd; ++i) J[i] = 0.0;

    for (Int a = 0; a < npe; ++a) {
        for (Int d = 0; d < nd; ++d) {
            const dstype xa = xdg_elem[a + npe * d];
            x[d] += N[a] * xa;
            for (Int k = 0; k < nd; ++k)
                J[d * nd + k] += dN[k * npe + a] * xa;
        }
    }
}

bool ComputeReferenceCoordinatesInElement(
    dstype* xi,
    const dstype* xphys,
    const dstype* xdg_elem,
    const masterstruct& master,
    const Int nd,
    const Int elemtype,
    const Int porder,
    const Int npe,
    const Int maxNewtonIter,
    const dstype tol)
{
    DefaultInitialGuess(xi, nd, elemtype);

    std::vector<dstype> xcur(static_cast<size_t>(nd));
    std::vector<dstype> rhs(static_cast<size_t>(nd));
    std::vector<dstype> delta(static_cast<size_t>(nd));
    std::vector<dstype> J(static_cast<size_t>(nd * nd));

    for (Int it = 0; it < maxNewtonIter; ++it) {
        MapReferenceToPhysical(xcur.data(), J.data(), xi, xdg_elem, master, nd, elemtype, porder, npe);

        for (Int d = 0; d < nd; ++d) rhs[d] = xphys[d] - xcur[d];
        if (VecNorm(rhs.data(), nd) <= tol) return true;

        if (!SolveSmallLinearSystem(delta.data(), J.data(), rhs.data(), nd))
            return false;

        for (Int d = 0; d < nd; ++d) xi[d] += delta[d];
        if (VecNorm(delta.data(), nd) <= tol) return true;
    }

    std::fill(xi, xi + nd, std::numeric_limits<dstype>::quiet_NaN());
    return false;
}

void EvaluateVolumeBasisAtXi(
    dstype* N,
    const dstype* xi,
    const masterstruct& master,
    const Int nd,
    const Int elemtype,
    const Int porder,
    const Int npe)
{
    std::vector<dstype> dN(static_cast<size_t>(nd * npe));
    EvaluateVolumeShapeAndGradientsAtXi(N, dN.data(), xi, master, nd, elemtype, porder, npe);
}

void EvaluateVolumeBasisBatch(
    std::vector<dstype>& shap1,
    const std::vector<dstype>& xi1,
    const masterstruct& master,
    const Int nd,
    const Int elemtype,
    const Int porder,
    const Int npe,
    const Int npoints)
{
    shap1.resize(static_cast<size_t>(npe * npoints));
    std::vector<dstype> N(static_cast<size_t>(npe));
    for (Int ip = 0; ip < npoints; ++ip) {
        EvaluateVolumeBasisAtXi(N.data(), &xi1[nd * ip], master, nd, elemtype, porder, npe);
        std::copy(N.begin(), N.end(), shap1.begin() + static_cast<size_t>(npe * ip));
    }
}

bool PointInReferenceElement(
    const dstype* xi,
    const Int nd,
    const Int elemtype,
    const dstype tol)
{
    if (elemtype == 0) {
        dstype sum = 0.0;
        for (Int d = 0; d < nd; ++d) {
            if (xi[d] < -tol) return false;
            sum += xi[d];
        }
        return (sum <= 1.0 + tol);
    }

    for (Int d = 0; d < nd; ++d) {
        if (xi[d] < -tol || xi[d] > 1.0 + tol) return false;
    }
    return true;
}

} // namespace wm
} // namespace exasim
