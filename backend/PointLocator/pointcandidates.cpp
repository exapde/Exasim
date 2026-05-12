#ifndef __POINTCANDIDATES_H__
#define __POINTCANDIDATES_H__

#include "../Common/common.h"

#include <cmath>
#include <limits>
#include <vector>

enum PointCandidateMethod {
    POINTCANDIDATES_FROM_FACE_F2E = 0,
    POINTCANDIDATES_FROM_ELEMENT_FACES = 1,
    POINTCANDIDATES_FROM_BOUNDING_BOXES = 2,
    POINTCANDIDATES_FROM_ELLIPSOIDS = 3
};

static inline bool CandidateListContains(const Int* candidateElems, Int ncandidates, Int elem)
{
    for (Int i = 0; i < ncandidates; ++i)
        if (candidateElems[i] == elem)
            return true;
    return false;
}

static inline Int AppendUniqueCandidate(Int* candidateElems, Int ncandidates, Int maxCandidates, Int elem)
{
    if (candidateElems == nullptr || maxCandidates <= 0) return ncandidates;
    if (elem < 0) return ncandidates;
    if (CandidateListContains(candidateElems, ncandidates, elem)) return ncandidates;
    if (ncandidates >= maxCandidates) return ncandidates;
    candidateElems[ncandidates] = elem;
    return ncandidates + 1;
}

static inline Int BuildAllCandidateElements(Int* candidateElems, Int maxCandidates, Int ne)
{
    if (candidateElems == nullptr || maxCandidates <= 0 || ne <= 0) return 0;

    const Int n = (ne < maxCandidates) ? ne : maxCandidates;
    for (Int e = 0; e < n; ++e)
        candidateElems[e] = e;
    return n;
}

static inline Int BuildCandidateElementsFromFaceF2E(
    Int* candidateElems,
    Int maxCandidates,
    Int face,
    const Int* f2e)
{
    if (candidateElems == nullptr || f2e == nullptr || maxCandidates <= 0 || face < 0)
        return 0;

    Int ncandidates = 0;
    ncandidates = AppendUniqueCandidate(candidateElems, ncandidates, maxCandidates, f2e[4 * face + 0]);
    ncandidates = AppendUniqueCandidate(candidateElems, ncandidates, maxCandidates, f2e[4 * face + 2]);
    return ncandidates;
}

static inline Int BuildCandidateElementsFromElementFaces(
    Int* candidateElems,
    Int maxCandidates,
    Int elem,
    const Int* e2f,
    const Int* f2e,
    Int nfe)
{
    if (candidateElems == nullptr || e2f == nullptr || f2e == nullptr || maxCandidates <= 0)
        return 0;
    if (elem < 0 || nfe <= 0) return 0;

    for (Int i = 0; i < maxCandidates; ++i)
        candidateElems[i] = -1;

    Int ncandidates = 0;
    ncandidates = AppendUniqueCandidate(candidateElems, ncandidates, maxCandidates, elem);

    Int ifront = 0;
    while (ifront < ncandidates && ncandidates < maxCandidates) {
        const Int ecur = candidateElems[ifront++];

        for (Int lf = 0; lf < nfe && ncandidates < maxCandidates; ++lf) {
            const Int face = e2f[ecur * nfe + lf];
            if (face < 0) continue;

            const Int e1 = f2e[4 * face + 0];
            const Int e2 = f2e[4 * face + 2];
            ncandidates = AppendUniqueCandidate(candidateElems, ncandidates, maxCandidates, e1);
            ncandidates = AppendUniqueCandidate(candidateElems, ncandidates, maxCandidates, e2);
        }
    }

    return ncandidates;
}

static inline void ComputeElementBoundingBox(
    dstype* xmin,
    dstype* xmax,
    const dstype* xdgElem,
    Int nd,
    Int npe)
{
    for (Int d = 0; d < nd; ++d) {
        xmin[d] = xdgElem[npe * d];
        xmax[d] = xdgElem[npe * d];
        for (Int a = 1; a < npe; ++a) {
            const dstype x = xdgElem[a + npe * d];
            if (x < xmin[d]) xmin[d] = x;
            if (x > xmax[d]) xmax[d] = x;
        }
    }
}

static inline bool InvertSmallMatrix(
    dstype* Ainv,
    dstype* A,
    Int nd)
{
    if (Ainv == nullptr || A == nullptr || nd <= 0) return false;

    for (Int i = 0; i < nd; ++i)
        for (Int j = 0; j < nd; ++j)
            Ainv[i * nd + j] = (i == j) ? static_cast<dstype>(1.0) : static_cast<dstype>(0.0);

    for (Int k = 0; k < nd; ++k) {
        Int ipiv = k;
        dstype pivmax = std::abs(A[k * nd + k]);
        for (Int i = k + 1; i < nd; ++i) {
            const dstype val = std::abs(A[i * nd + k]);
            if (val > pivmax) {
                pivmax = val;
                ipiv = i;
            }
        }
        if (pivmax <= std::numeric_limits<dstype>::epsilon()) return false;

        if (ipiv != k) {
            for (Int j = 0; j < nd; ++j) {
                std::swap(A[k * nd + j], A[ipiv * nd + j]);
                std::swap(Ainv[k * nd + j], Ainv[ipiv * nd + j]);
            }
        }

        const dstype piv = A[k * nd + k];
        for (Int j = 0; j < nd; ++j) {
            A[k * nd + j] /= piv;
            Ainv[k * nd + j] /= piv;
        }

        for (Int i = 0; i < nd; ++i) {
            if (i == k) continue;
            const dstype factor = A[i * nd + k];
            if (factor == static_cast<dstype>(0.0)) continue;
            for (Int j = 0; j < nd; ++j) {
                A[i * nd + j] -= factor * A[k * nd + j];
                Ainv[i * nd + j] -= factor * Ainv[k * nd + j];
            }
        }
    }
    return true;
}

static inline Int GetEnclosingEllipsoidScratchSize(Int nd, Int npe)
{
    const Int m = nd + 1;
    return 2 * npe + 2 * m * m + m + 2 * nd * nd;
}

static inline void ComputeEnclosingEllipsoid(
    dstype* center,
    dstype* Aell,
    const dstype* xdgElem,
    Int nd,
    Int npe,
    dstype* tm,
    Int maxIter = 200,
    dstype tol = static_cast<dstype>(1.0e-8))
{
    const dstype eps = static_cast<dstype>(1.0e-12);
    const Int m = nd + 1;

    if (tm == nullptr) return;

    dstype* u    = &tm[0];
    dstype* unew = &tm[npe];
    dstype* X    = &tm[2 * npe];
    dstype* Xinv = &tm[2 * npe + m * m];
    dstype* q    = &tm[2 * npe + 2 * m * m];
    dstype* S    = &tm[2 * npe + 2 * m * m + m];
    dstype* Sinv = &tm[2 * npe + 2 * m * m + m + nd * nd];

    for (Int a = 0; a < npe; ++a) {
        u[a] = static_cast<dstype>(1.0) / static_cast<dstype>(npe);
        unew[a] = 0.0;
    }

    for (Int it = 0; it < maxIter; ++it) {
        for (Int i = 0; i < m * m; ++i) X[i] = 0.0;

        for (Int a = 0; a < npe; ++a) {
            for (Int d = 0; d < nd; ++d) q[d] = xdgElem[a + npe * d];
            q[nd] = static_cast<dstype>(1.0);

            for (Int i = 0; i < m; ++i)
                for (Int j = 0; j < m; ++j)
                    X[i * m + j] += u[a] * q[i] * q[j];
        }

        if (!InvertSmallMatrix(Xinv, X, m)) break;

        Int imax = 0;
        dstype mmax = -1.0;
        for (Int a = 0; a < npe; ++a) {
            for (Int d = 0; d < nd; ++d) q[d] = xdgElem[a + npe * d];
            q[nd] = static_cast<dstype>(1.0);

            dstype Mi = 0.0;
            for (Int i = 0; i < m; ++i)
                for (Int j = 0; j < m; ++j)
                    Mi += q[i] * Xinv[i * m + j] * q[j];

            if (Mi > mmax) {
                mmax = Mi;
                imax = a;
            }
        }

        if (mmax <= static_cast<dstype>(m) * (static_cast<dstype>(1.0) + tol))
            break;

        dstype step = (mmax - static_cast<dstype>(m)) /
                       (static_cast<dstype>(m) * (mmax - static_cast<dstype>(1.0)));
        if (step < 0.0) step = 0.0;
        if (step > 1.0) step = 1.0;

        for (Int a = 0; a < npe; ++a) unew[a] = (static_cast<dstype>(1.0) - step) * u[a];
        unew[imax] += step;
        for (Int a = 0; a < npe; ++a) {
            const dstype tmp = u[a];
            u[a] = unew[a];
            unew[a] = tmp;
        }
    }

    for (Int d = 0; d < nd; ++d) {
        center[d] = 0.0;
        for (Int a = 0; a < npe; ++a)
            center[d] += u[a] * xdgElem[a + npe * d];
    }

    for (Int i = 0; i < nd * nd; ++i) S[i] = 0.0;
    for (Int a = 0; a < npe; ++a) {
        for (Int i = 0; i < nd; ++i) {
            const dstype di = xdgElem[a + npe * i] - center[i];
            for (Int j = 0; j < nd; ++j) {
                const dstype dj = xdgElem[a + npe * j] - center[j];
                S[i * nd + j] += u[a] * di * dj;
            }
        }
    }

    for (Int d = 0; d < nd; ++d)
        S[d * nd + d] += eps;

    if (!InvertSmallMatrix(Sinv, S, nd)) {
        for (Int i = 0; i < nd * nd; ++i) Aell[i] = 0.0;
        for (Int d = 0; d < nd; ++d) Aell[d * nd + d] = static_cast<dstype>(1.0) / eps;
        return;
    }

    const dstype scale = static_cast<dstype>(1.0) / static_cast<dstype>(nd);
    for (Int i = 0; i < nd * nd; ++i)
        Aell[i] = scale * Sinv[i];
}

static inline void BuildElementEllipsoids(
    dstype* elemEllipsoids,
    const dstype* xdg,
    Int nd,
    Int npe,
    Int ncx,
    Int ne)
{
    if (elemEllipsoids == nullptr || xdg == nullptr) return;

    std::vector<dstype> xdgElem(static_cast<size_t>(npe * nd));
    std::vector<dstype> tm(static_cast<size_t>(GetEnclosingEllipsoidScratchSize(nd, npe)), 0.0);
    for (Int e = 0; e < ne; ++e) {
        const Int offset = npe * ncx * e;
        for (Int d = 0; d < nd; ++d)
            for (Int a = 0; a < npe; ++a)
                xdgElem[a + npe * d] = xdg[offset + a + npe * d];

        dstype* center = &elemEllipsoids[(nd + nd * nd) * e];
        dstype* Aell   = &elemEllipsoids[(nd + nd * nd) * e + nd];
        ComputeEnclosingEllipsoid(center, Aell, xdgElem.data(), nd, npe, tm.data());
    }
}

static inline bool PointInBoundingBox(
    const dstype* xphys,
    const dstype* xmin,
    const dstype* xmax,
    Int nd,
    dstype tol)
{
    for (Int d = 0; d < nd; ++d) {
        if (xphys[d] < xmin[d] - tol || xphys[d] > xmax[d] + tol)
            return false;
    }
    return true;
}

static inline bool PointInEllipsoid(
    const dstype* xphys,
    const dstype* center,
    const dstype* Aell,
    Int nd,
    dstype tol)
{
    dstype s = 0.0;
    for (Int i = 0; i < nd; ++i) {
        const dstype dxi = xphys[i] - center[i];
        for (Int j = 0; j < nd; ++j)
            s += dxi * Aell[i * nd + j] * (xphys[j] - center[j]);
    }

    return (s <= static_cast<dstype>(1.0) + tol);
}

static inline Int BuildCandidateElementsFromBoundingBoxes(
    Int* candidateElems,
    Int maxCandidates,
    const dstype* xphys,
    const dstype* elemBoxes,
    Int nd,
    Int ne,
    dstype tol)
{
    if (candidateElems == nullptr || xphys == nullptr || elemBoxes == nullptr || maxCandidates <= 0)
        return 0;

    Int ncandidates = 0;
    for (Int e = 0; e < ne; ++e) {
        const dstype* xmin = &elemBoxes[(2 * nd) * e];
        const dstype* xmax = &elemBoxes[(2 * nd) * e + nd];
        if (PointInBoundingBox(xphys, xmin, xmax, nd, tol)) {
            ncandidates = AppendUniqueCandidate(candidateElems, ncandidates, maxCandidates, e);
            if (ncandidates >= maxCandidates) break;
        }
    }
    return ncandidates;
}

static inline Int BuildCandidateElementsFromEllipsoids(
    Int* candidateElems,
    Int maxCandidates,
    const dstype* xphys,
    const dstype* elemEllipsoids,
    Int nd,
    Int ne,
    dstype tol)
{
    if (candidateElems == nullptr || xphys == nullptr || elemEllipsoids == nullptr || maxCandidates <= 0)
        return 0;

    Int ncandidates = 0;
    for (Int e = 0; e < ne; ++e) {
        const dstype* center = &elemEllipsoids[(nd + nd * nd) * e];
        const dstype* Aell   = &elemEllipsoids[(nd + nd * nd) * e + nd];
        if (PointInEllipsoid(xphys, center, Aell, nd, tol)) {
            ncandidates = AppendUniqueCandidate(candidateElems, ncandidates, maxCandidates, e);
            if (ncandidates >= maxCandidates) break;
        }
    }
    return ncandidates;
}

static inline void BuildElementBoundingBoxes(
    dstype* elemBoxes,
    const dstype* xdg,
    Int nd,
    Int npe,
    Int ncx,
    Int ne)
{
    if (elemBoxes == nullptr || xdg == nullptr) return;

    for (Int e = 0; e < ne; ++e) {
        const Int offset = npe * ncx * e;
        dstype* xmin = &elemBoxes[(2 * nd) * e];
        dstype* xmax = &elemBoxes[(2 * nd) * e + nd];

        for (Int d = 0; d < nd; ++d) {
            xmin[d] = xdg[offset + npe * d];
            xmax[d] = xdg[offset + npe * d];
            for (Int a = 1; a < npe; ++a) {
                const dstype x = xdg[offset + a + npe * d];
                if (x < xmin[d]) xmin[d] = x;
                if (x > xmax[d]) xmax[d] = x;
            }
        }
    }
}

static inline Int BuildCandidateElements(
    Int* candidateElems,
    Int maxCandidates,
    Int method,
    Int face,
    Int elem,
    const dstype* xphys,
    const Int* f2e,
    const Int* e2f,
    Int nfe,
    const dstype* elemBoxes,
    const dstype* elemEllipsoids,
    Int nd,
    Int ne,
    dstype tol)
{
    switch (method) {
        case POINTCANDIDATES_FROM_FACE_F2E:
            return BuildCandidateElementsFromFaceF2E(candidateElems, maxCandidates, face, f2e);

        case POINTCANDIDATES_FROM_ELEMENT_FACES:
            return BuildCandidateElementsFromElementFaces(
                candidateElems, maxCandidates, elem, e2f, f2e, nfe);

        case POINTCANDIDATES_FROM_BOUNDING_BOXES:
            return BuildCandidateElementsFromBoundingBoxes(
                candidateElems, maxCandidates, xphys, elemBoxes, nd, ne, tol);

        case POINTCANDIDATES_FROM_ELLIPSOIDS:
            return BuildCandidateElementsFromEllipsoids(
                candidateElems, maxCandidates, xphys, elemEllipsoids, nd, ne, tol);

        default:
            return 0;
    }
}

#endif
