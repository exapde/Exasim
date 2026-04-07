#ifndef __POINTLOCATION_H__
#define __POINTLOCATION_H__

#include "../Common/common.h"
#include "../Preprocessing/makemaster.cpp"
#include "pointellipsoidgrid.cpp"
#include "pointcandidates.cpp"
#include "pointinterpolation.cpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

static inline bool PointInReferenceElement(const dstype* xi, Int nd, Int elemtype, dstype tol)
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

static inline void ProjectPointToReferenceElement(
    dstype* xproj,
    const dstype* xi,
    dstype* u,
    Int nd,
    Int elemtype)
{
    if (elemtype == 0) {
        if (xproj == nullptr || xi == nullptr || u == nullptr || nd <= 0)
            return;

        dstype sumpos = 0.0;
        for (Int d = 0; d < nd; ++d) {
            xproj[d] = (xi[d] > 0.0) ? xi[d] : 0.0;
            sumpos += xproj[d];
        }

        if (sumpos <= 1.0) {
            return;
        }

        for (Int d = 0; d < nd; ++d)
            u[d] = xproj[d];
        std::sort(u, u + nd, std::greater<dstype>());

        dstype csum = 0.0;
        dstype theta = 0.0;
        for (Int i = 0; i < nd; ++i) {
            csum += u[i];
            const dstype t = (csum - 1.0) / static_cast<dstype>(i + 1);
            if (u[i] > t) theta = t;
        }

        for (Int d = 0; d < nd; ++d)
            xproj[d] = std::max(xproj[d] - theta, static_cast<dstype>(0.0));
        return;
    }

    if (xproj == nullptr || xi == nullptr || nd <= 0)
        return;

    for (Int d = 0; d < nd; ++d) {
        if (xi[d] < 0.0) xproj[d] = 0.0;
        else if (xi[d] > 1.0) xproj[d] = 1.0;
        else xproj[d] = xi[d];
    }
}

static inline dstype ReferenceElementDistance(
    const dstype* xi,
    dstype* xproj,
    dstype* u,
    Int nd,
    Int elemtype)
{
    if (xi == nullptr || xproj == nullptr || nd <= 0)
        return std::numeric_limits<dstype>::max();

    ProjectPointToReferenceElement(xproj, xi, u, nd, elemtype);

    dstype dist2 = 0.0;
    for (Int d = 0; d < nd; ++d) {
        const dstype diff = xi[d] - xproj[d];
        dist2 += diff * diff;
    }
    return std::sqrt(dist2);
}

static inline void ExtractElementGeometry(dstype* xdgElem, const dstype* xdg, Int elem, Int npe, Int ncx, Int nd)
{
    const Int offset = npe * ncx * elem;
    for (Int d = 0; d < nd; ++d)
        for (Int a = 0; a < npe; ++a)
            xdgElem[a + npe * d] = xdg[offset + a + npe * d];
}

static inline dstype PointLocationNorm(const dstype* a, Int n)
{
    dstype s = 0.0;
    for (Int i = 0; i < n; ++i) s += a[i] * a[i];
    return std::sqrt(s);
}

static inline bool SolveSmallLinearSystem(dstype* x, const dstype* A, const dstype* b, Int nd)
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

static inline void DefaultReferenceGuess(dstype* xi, Int nd, Int elemtype)
{
    if (elemtype == 0) {
        const dstype val = 1.0 / static_cast<dstype>(nd + 1);
        for (Int d = 0; d < nd; ++d) xi[d] = val;
    }
    else {
        for (Int d = 0; d < nd; ++d) xi[d] = 0.5;
    }
}

static inline void MapReferenceToPhysical(
    dstype* x,
    dstype* J,
    dstype* N,
    dstype* dN,
    const dstype* xi,
    const dstype* xdgElem,
    const dstype* xpe,
    Int nd,
    Int npe,
    Int elemtype,
    Int porder)
{
    if (x == nullptr || J == nullptr || N == nullptr || dN == nullptr ||
        xi == nullptr || xdgElem == nullptr || xpe == nullptr || nd <= 0 || npe <= 0)
        return;

    EvaluateVolumeShapeAndGradientsAtReferencePoint(
        N, dN, xi, xpe, nd, npe, elemtype, porder);

    for (Int d = 0; d < nd; ++d) x[d] = 0.0;
    for (Int i = 0; i < nd * nd; ++i) J[i] = 0.0;

    for (Int a = 0; a < npe; ++a) {
        for (Int d = 0; d < nd; ++d) {
            const dstype xa = xdgElem[a + npe * d];
            x[d] += N[a] * xa;
            for (Int k = 0; k < nd; ++k)
                J[d * nd + k] += dN[k * npe + a] * xa;
        }
    }
}

static inline Int GetReferenceCoordinateScratchSize(Int nd, Int npe)
{
    return 3 * nd + nd * nd + npe + nd * npe;
}

static inline bool ComputeReferenceCoordinates(
    dstype* xi,
    const dstype* xphys,
    const dstype* xdgElem,
    const dstype* xpe,
    Int nd,
    Int npe,
    dstype* tm,
    Int elemtype,
    Int porder,
    Int maxNewtonIter,
    dstype newtonTol)
{
    if (xi == nullptr || xphys == nullptr || xdgElem == nullptr || xpe == nullptr ||
        tm == nullptr || nd <= 0 || npe <= 0)
        return false;

    DefaultReferenceGuess(xi, nd, elemtype);

    dstype* xcur = &tm[0];
    dstype* rhs = &tm[nd];
    dstype* delta = &tm[2 * nd];
    dstype* J = &tm[3 * nd];
    dstype* N = &tm[3 * nd + nd * nd];
    dstype* dN = &tm[3 * nd + nd * nd + npe];

    for (Int it = 0; it < maxNewtonIter; ++it) {
        MapReferenceToPhysical(
            xcur, J, N, dN, xi, xdgElem, xpe,
            nd, npe, elemtype, porder);

        for (Int d = 0; d < nd; ++d) rhs[d] = xphys[d] - xcur[d];
        if (PointLocationNorm(rhs, nd) <= newtonTol) return true;

        if (!SolveSmallLinearSystem(delta, J, rhs, nd))
            return false;

        for (Int d = 0; d < nd; ++d) xi[d] += delta[d];
        if (PointLocationNorm(delta, nd) <= newtonTol) return true;
    }
    return false;
}

static inline bool FindPointInElement(
    Int* elem,
    dstype* xi,
    const dstype* xphys,
    const dstype* xdg,
    const Int* candidateElems,
    Int ncandidates,
    const dstype* xpe,
    Int nd,
    Int npe,
    Int ncx,
    Int elemtype,
    Int porder,
    Int maxNewtonIter,
    dstype newtonTol,
    dstype insideTol)
{
    if (elem == nullptr || xi == nullptr || xphys == nullptr || xdg == nullptr ||
        candidateElems == nullptr || xpe == nullptr || ncandidates <= 0)
        return false;

    *elem = -1;
    std::vector<dstype> xdgElem(static_cast<size_t>(npe * nd));
    std::vector<dstype> bestXi(static_cast<size_t>(nd), 0.0);
    std::vector<dstype> xproj(static_cast<size_t>(nd), 0.0);
    std::vector<dstype> uproj(static_cast<size_t>(nd), 0.0);
    std::vector<dstype> tm(static_cast<size_t>(GetReferenceCoordinateScratchSize(nd, npe)), 0.0);
    dstype bestDist = std::numeric_limits<dstype>::max();
    bool haveBest = false;

    for (Int ic = 0; ic < ncandidates; ++ic) {
        const Int e = candidateElems[ic];
        if (e < 0) continue;

        ExtractElementGeometry(xdgElem.data(), xdg, e, npe, ncx, nd);
        if (!ComputeReferenceCoordinates(
                xi, xphys, xdgElem.data(), xpe, nd, npe, tm.data(), elemtype, porder,
                maxNewtonIter, newtonTol))
            continue;

        if (PointInReferenceElement(xi, nd, elemtype, insideTol)) {
            *elem = e;
            return true;
        }

        const dstype dist = ReferenceElementDistance(xi, xproj.data(), uproj.data(), nd, elemtype);
        if (!haveBest || dist < bestDist) {
            bestDist = dist;
            *elem = e;
            for (Int d = 0; d < nd; ++d) bestXi[d] = xi[d];
            haveBest = true;
        }
    }

    if (haveBest) {
        for (Int d = 0; d < nd; ++d) xi[d] = bestXi[d];
    }
    return false;
}

static inline bool FindPointInElementIterative(
    Int* elem,
    dstype* xi,
    const dstype* xphys,
    const dstype* xdg,
    Int* candidateElems,
    Int ncandidates,
    Int maxCandidates,
    const Int* e2f,
    const Int* f2e,
    Int nfe,
    const dstype* xpe,
    Int nd,
    Int npe,
    Int ncx,
    Int ne,
    Int elemtype,
    Int porder,
    Int maxNewtonIter,
    dstype newtonTol,
    dstype insideTol,
    bool allowGlobalFallback)
{
    if (elem == nullptr || xi == nullptr || xphys == nullptr || xdg == nullptr ||
        candidateElems == nullptr || xpe == nullptr || maxCandidates <= 0)
        return false;

    *elem = -1;
    std::vector<Int> current(static_cast<size_t>(maxCandidates), -1);
    std::vector<Int> next(static_cast<size_t>(maxCandidates), -1);
    std::vector<Int> visited(static_cast<size_t>(maxCandidates), -1);
    std::vector<dstype> bestXi(static_cast<size_t>(nd), 0.0);
    std::vector<dstype> xproj(static_cast<size_t>(nd), 0.0);
    std::vector<dstype> uproj(static_cast<size_t>(nd), 0.0);
    Int bestElem = -1;
    dstype bestDist = std::numeric_limits<dstype>::max();
    bool haveBest = false;

    Int ncurrent = 0;
    Int nvisited = 0;
    for (Int i = 0; i < ncandidates && i < maxCandidates; ++i) {
        ncurrent = AppendUniqueCandidate(current.data(), ncurrent, maxCandidates, candidateElems[i]);
        nvisited = AppendUniqueCandidate(visited.data(), nvisited, maxCandidates, candidateElems[i]);
    }

    while (ncurrent > 0) {
        if (FindPointInElement(
                elem, xi, xphys, xdg, current.data(), ncurrent, xpe, nd, npe, ncx,
                elemtype, porder, maxNewtonIter, newtonTol, insideTol))
            return true;

        if (*elem >= 0) {
            const dstype dist = ReferenceElementDistance(xi, xproj.data(), uproj.data(), nd, elemtype);
            if (!haveBest || dist < bestDist) {
                bestDist = dist;
                bestElem = *elem;
                for (Int d = 0; d < nd; ++d) bestXi[d] = xi[d];
                haveBest = true;
            }
        }

        Int nnext = 0;
        std::fill(next.begin(), next.end(), -1);
        for (Int i = 0; i < ncurrent; ++i) {
            const Int e = current[i];
            std::vector<Int> local(static_cast<size_t>(maxCandidates), -1);
            const Int nlocal = BuildCandidateElementsFromElementFaces(
                local.data(), maxCandidates, e, e2f, f2e, nfe > 0 ? nfe : 0);
            for (Int j = 0; j < nlocal; ++j)
                nnext = AppendUniqueCandidate(next.data(), nnext, maxCandidates, local[j]);
        }

        // Remove elements already tested in the current frontier.
        Int nfiltered = 0;
        std::vector<Int> filtered(static_cast<size_t>(maxCandidates), -1);
        for (Int i = 0; i < nnext; ++i) {
            const Int e = next[i];
            if (!CandidateListContains(visited.data(), nvisited, e))
                nfiltered = AppendUniqueCandidate(filtered.data(), nfiltered, maxCandidates, e);
        }

        if (nfiltered == 0) break;
        for (Int i = 0; i < nfiltered; ++i)
            nvisited = AppendUniqueCandidate(visited.data(), nvisited, maxCandidates, filtered[i]);
        current.swap(filtered);
        ncurrent = nfiltered;
    }

    if (allowGlobalFallback) {
        std::vector<Int> all(static_cast<size_t>(maxCandidates), -1);
        const Int nall = BuildAllCandidateElements(all.data(), maxCandidates, ne);
        if (FindPointInElement(
                elem, xi, xphys, xdg, all.data(), nall, xpe, nd, npe, ncx,
                elemtype, porder, maxNewtonIter, newtonTol, insideTol))
            return true;

        if (*elem >= 0) {
            const dstype dist = ReferenceElementDistance(xi, xproj.data(), uproj.data(), nd, elemtype);
            if (!haveBest || dist < bestDist) {
                bestDist = dist;
                bestElem = *elem;
                for (Int d = 0; d < nd; ++d) bestXi[d] = xi[d];
                haveBest = true;
            }
        }
    }

    if (haveBest) {
        *elem = bestElem;
        for (Int d = 0; d < nd; ++d) xi[d] = bestXi[d];
    }
    else {
        *elem = -1;
    }
    return false;
}

static inline bool FindPointAndShapeFunctionsFromElementFaces(
    Int* elem,
    dstype* xi,
    dstype* N,
    const dstype* xphys,
    const dstype* xdg,
    Int seedElem,
    Int* candidateElems,
    Int maxCandidates,
    const Int* e2f,
    const Int* f2e,
    Int nfe,
    const dstype* xpe,
    Int nd,
    Int npe,
    Int ncx,
    Int elemtype,
    Int porder,
    Int maxNewtonIter,
    dstype newtonTol,
    dstype insideTol)
{
    if (elem == nullptr || xi == nullptr || N == nullptr || xphys == nullptr || xdg == nullptr ||
        candidateElems == nullptr || e2f == nullptr || f2e == nullptr || xpe == nullptr ||
        maxCandidates <= 0 || nd <= 0 || npe <= 0 || ncx <= 0)
        return false;

    const Int ncandidates = BuildCandidateElementsFromElementFaces(
        candidateElems, maxCandidates, seedElem, e2f, f2e, nfe);

    const bool found = FindPointInElement(
        elem, xi, xphys, xdg, candidateElems, ncandidates, xpe, nd, npe, ncx,
        elemtype, porder, maxNewtonIter, newtonTol, insideTol);

    if (*elem >= 0)
        EvaluateVolumeShapeAtReferencePoint(N, xi, xpe, nd, npe, elemtype, porder);

    return found;
}

static inline bool FindPointsAndShapeFunctionsFromElementFaces(
    Int* elem,
    dstype* xi,
    dstype* N,
    const dstype* xphys,
    Int npoints,
    const dstype* xdg,
    const Int* seedElems,
    Int* candidateElems,
    Int maxCandidates,
    const Int* e2f,
    const Int* f2e,
    Int nfe,
    const dstype* xpe,
    Int nd,
    Int npe,
    Int ncx,
    Int elemtype,
    Int porder,
    Int maxNewtonIter,
    dstype newtonTol,
    dstype insideTol)
{
    if (elem == nullptr || xi == nullptr || N == nullptr || xphys == nullptr ||
        xdg == nullptr || seedElems == nullptr || candidateElems == nullptr ||
        e2f == nullptr || f2e == nullptr || xpe == nullptr ||
        npoints <= 0 || maxCandidates <= 0 || nd <= 0 || npe <= 0 || ncx <= 0)
        return false;

    bool allFound = true;
    for (Int ip = 0; ip < npoints; ++ip) {
        const bool found = FindPointAndShapeFunctionsFromElementFaces(
            &elem[ip], &xi[nd * ip], &N[npe * ip], &xphys[nd * ip],
            xdg, seedElems[ip], candidateElems, maxCandidates,
            e2f, f2e, nfe, xpe, nd, npe, ncx, elemtype, porder,
            maxNewtonIter, newtonTol, insideTol);
        allFound = allFound && found;
    }

    return allFound;
}

static inline bool FindPointShapeAndFieldFromEllipsoidGrid(
    Int* elem,
    dstype* xi,
    dstype* N,
    dstype* u,
    const dstype* xphys,
    const dstype* U,
    const dstype* xdg,
    Int* candidateElems,
    Int maxCandidates,
    const Int* binOffsets,
    const Int* binElems,
    const dstype* gridMin,
    const dstype* gridH,
    const Int* nbins,
    const dstype* elemEllipsoids,
    const dstype* xpe,
    Int nd,
    Int npe,
    Int nc,
    Int ncx,
    Int elemtype,
    Int porder,
    Int maxNewtonIter,
    dstype newtonTol,
    dstype insideTol,
    dstype ellipTol)
{
    if (elem == nullptr || xi == nullptr || N == nullptr || u == nullptr ||
        xphys == nullptr || U == nullptr || xdg == nullptr || candidateElems == nullptr ||
        binOffsets == nullptr || binElems == nullptr || gridMin == nullptr ||
        gridH == nullptr || nbins == nullptr || elemEllipsoids == nullptr ||
        xpe == nullptr || maxCandidates <= 0 || nd <= 0 || npe <= 0 || nc <= 0 || ncx <= 0)
        return false;

    const Int ncandidates = BuildCandidateElementsFromEllipsoidGridExact(
        candidateElems, maxCandidates, xphys, binOffsets, binElems, gridMin, gridH,
        nbins, elemEllipsoids, nd, ellipTol);

    const bool found = FindPointInElement(
        elem, xi, xphys, xdg, candidateElems, ncandidates, xpe, nd, npe, ncx,
        elemtype, porder, maxNewtonIter, newtonTol, insideTol);

    if (*elem >= 0) {
        EvaluateVolumeShapeAtReferencePoint(N, xi, xpe, nd, npe, elemtype, porder);
        const dstype* Uelem = &U[npe * nc * (*elem)];
        InterpolateFieldAtReferencePoint(u, Uelem, N, npe, nc);
    }

    return found;
}

static inline bool FindPointsShapesAndFieldsFromEllipsoidGrid(
    Int* elem,
    dstype* xi,
    dstype* N,
    dstype* u,
    const dstype* xphys,
    Int npoints,
    const dstype* U,
    const dstype* xdg,
    Int* candidateElems,
    Int maxCandidates,
    const Int* binOffsets,
    const Int* binElems,
    const dstype* gridMin,
    const dstype* gridH,
    const Int* nbins,
    const dstype* elemEllipsoids,
    const dstype* xpe,
    Int nd,
    Int npe,
    Int nc,
    Int ncx,
    Int elemtype,
    Int porder,
    Int maxNewtonIter,
    dstype newtonTol,
    dstype insideTol,
    dstype ellipTol)
{
    if (elem == nullptr || xi == nullptr || N == nullptr || u == nullptr ||
        xphys == nullptr || U == nullptr || xdg == nullptr || candidateElems == nullptr ||
        binOffsets == nullptr || binElems == nullptr || gridMin == nullptr ||
        gridH == nullptr || nbins == nullptr || elemEllipsoids == nullptr ||
        xpe == nullptr || npoints <= 0 || maxCandidates <= 0 ||
        nd <= 0 || npe <= 0 || nc <= 0 || ncx <= 0)
        return false;

    bool allFound = true;
    for (Int ip = 0; ip < npoints; ++ip) {
        const bool found = FindPointShapeAndFieldFromEllipsoidGrid(
            &elem[ip], &xi[nd * ip], &N[npe * ip], &u[nc * ip], &xphys[nd * ip],
            U, xdg, candidateElems, maxCandidates, binOffsets, binElems,
            gridMin, gridH, nbins, elemEllipsoids, xpe, nd, npe, nc, ncx,
            elemtype, porder, maxNewtonIter, newtonTol, insideTol, ellipTol);
        allFound = allFound && found;
    }

    return allFound;
}

#endif
