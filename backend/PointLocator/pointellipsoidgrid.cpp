#ifndef __ELLIPSOIDGRID_H__
#define __ELLIPSOIDGRID_H__

#include "../Common/common.h"
#include "pointcandidates.cpp"

/*
 * Uniform-grid broad-phase index for element ellipsoids.
 *
 * Data layout
 * -----------
 *
 * 1. elemEllipsoids
 *    Per element e, the storage is
 *
 *      [c_0, ..., c_{nd-1}, A_00, A_01, ..., A_{nd-1,nd-1}]
 *
 *    with stride
 *
 *      ellStride = nd + nd*nd
 *
 *    so
 *
 *      center = &elemEllipsoids[ellStride * e]
 *      Aell   = &elemEllipsoids[ellStride * e + nd]
 *
 * 2. ellipBoxes
 *    Per element e, the enclosing axis-aligned bounding box is
 *
 *      [xmin_0, ..., xmin_{nd-1}, xmax_0, ..., xmax_{nd-1}]
 *
 *    with stride
 *
 *      boxStride = 2*nd
 *
 *    so
 *
 *      xmin = &ellipBoxes[boxStride * e]
 *      xmax = &ellipBoxes[boxStride * e + nd]
 *
 * 3. Uniform grid geometry
 *
 *      gridMin[nd]   lower corner of the grid
 *      gridH[nd]     bin spacing in each coordinate direction
 *      nbins[nd]     number of bins in each direction
 *
 *    The total number of bins is
 *
 *      nbinsTotal = prod_{d=0}^{nd-1} nbins[d]
 *
 * 4. Bin-to-ellipsoid adjacency in CSR form
 *
 *      binOffsets[nbinsTotal + 1]
 *      binElems[nbinElems]
 *
 *    For a flat bin index b, the ellipsoid ids in that bin are
 *
 *      binElems[ binOffsets[b] : binOffsets[b+1]-1 ]
 */

static inline Int EllipsoidStride(Int nd)
{
    return nd + nd * nd;
}

static inline Int EllipsoidBoxStride(Int nd)
{
    return 2 * nd;
}

static inline Int EllipsoidGridTotalBins(const Int* nbins, Int nd)
{
    if (nbins == nullptr || nd <= 0) return 0;

    Int ntotal = 1;
    for (Int d = 0; d < nd; ++d)
        ntotal *= nbins[d];
    return ntotal;
}

/*
 * Flatten a multi-index ibin[0:nd-1] into a single bin id using row-major
 * ordering:
 *
 *   flat = ibin[0]
 *        + nbins[0] * (ibin[1]
 *        + nbins[1] * (ibin[2] + ... ))
 */
static inline Int EllipsoidGridFlattenIndex(const Int* ibin, const Int* nbins, Int nd)
{
    if (ibin == nullptr || nbins == nullptr || nd <= 0) return -1;

    Int flat = 0;
    Int stride = 1;
    for (Int d = 0; d < nd; ++d) {
        flat += ibin[d] * stride;
        stride *= nbins[d];
    }
    return flat;
}

/*
 * Inverse of EllipsoidGridFlattenIndex.
 */
static inline void EllipsoidGridUnflattenIndex(Int* ibin, Int flat, const Int* nbins, Int nd)
{
    if (ibin == nullptr || nbins == nullptr || nd <= 0) return;

    for (Int d = 0; d < nd; ++d) {
        ibin[d] = flat % nbins[d];
        flat /= nbins[d];
    }
}

/*
 * Scratch size for BuildEllipsoidBoundingBoxes.
 * A minimal implementation needs space for a working copy of A and its inverse.
 */
static inline Int GetEllipsoidBoundingBoxScratchSize(Int nd)
{
    return 2 * nd * nd;
}

/*
 * Integer scratch size for building the CSR grid connectivity.
 * A simple two-pass implementation needs
 *
 *   binCounts[nbinsTotal] and binCursor[nbinsTotal].
 */
static inline Int GetEllipsoidUniformGridIntScratchSize(const Int* nbins, Int nd)
{
    const Int nbinsTotal = EllipsoidGridTotalBins(nbins, nd);
    return 2 * nbinsTotal;
}

/*
 * Build axis-aligned bounding boxes for the ellipsoids.
 *
 * Inputs:
 *   elemEllipsoids : [ne * (nd + nd*nd)]
 *   nd, ne
 *   tm             : scratch of size GetEllipsoidBoundingBoxScratchSize(nd)
 *
 * Output:
 *   ellipBoxes     : [ne * 2*nd]
 */
static inline void BuildEllipsoidBoundingBoxes(
    dstype* ellipBoxes,
    const dstype* elemEllipsoids,
    Int nd,
    Int ne,
    dstype* tm)
{
    if (ellipBoxes == nullptr || elemEllipsoids == nullptr || tm == nullptr || nd <= 0 || ne <= 0)
        return;

    const Int ellStride = EllipsoidStride(nd);
    const Int boxStride = EllipsoidBoxStride(nd);
    dstype* Awork = &tm[0];
    dstype* Ainv = &tm[nd * nd];

    for (Int e = 0; e < ne; ++e) {
        const dstype* center = &elemEllipsoids[ellStride * e];
        const dstype* Aell = &elemEllipsoids[ellStride * e + nd];
        dstype* xmin = &ellipBoxes[boxStride * e];
        dstype* xmax = &ellipBoxes[boxStride * e + nd];

        for (Int i = 0; i < nd * nd; ++i)
            Awork[i] = Aell[i];

        const bool ok = InvertSmallMatrix(Ainv, Awork, nd);
        for (Int d = 0; d < nd; ++d) {
            dstype rad = static_cast<dstype>(0.0);
            if (ok) {
                const dstype val = Ainv[d * nd + d];
                rad = (val > static_cast<dstype>(0.0)) ? std::sqrt(val) : static_cast<dstype>(0.0);
            }
            xmin[d] = center[d] - rad;
            xmax[d] = center[d] + rad;
        }
    }
}

/*
 * Compute the global bounds of the uniform grid from the element geometry xdg.
 *
 * Inputs:
 *   xdg : [ne * npe * ncx]
 *
 * Outputs:
 *   gridMin[nd], gridMax[nd]
 */
static inline void ComputeDomainGridBounds(
    dstype* gridMin,
    dstype* gridMax,
    const dstype* xdg,
    Int nd,
    Int npe,
    Int ncx,
    Int ne)
{
    if (gridMin == nullptr || gridMax == nullptr || xdg == nullptr ||
        nd <= 0 || npe <= 0 || ncx <= 0 || ne <= 0)
        return;

    const Int offset0 = 0;
    for (Int d = 0; d < nd; ++d)
        gridMin[d] = gridMax[d] = xdg[offset0 + npe * d];

    for (Int e = 0; e < ne; ++e) {
        const Int offset = npe * ncx * e;
        for (Int d = 0; d < nd; ++d) {
            for (Int a = 0; a < npe; ++a) {
                const dstype x = xdg[offset + a + npe * d];
                if (x < gridMin[d]) gridMin[d] = x;
                if (x > gridMax[d]) gridMax[d] = x;
            }
        }
    }
}

/*
 * Compute grid spacing from user-selected nbins.
 *
 * Inputs:
 *   gridMin[nd], gridMax[nd], nbins[nd]
 *
 * Output:
 *   gridH[nd]
 */
static inline void ComputeEllipsoidGridSpacing(
    dstype* gridH,
    const dstype* gridMin,
    const dstype* gridMax,
    const Int* nbins,
    Int nd)
{
    if (gridH == nullptr || gridMin == nullptr || gridMax == nullptr || nbins == nullptr || nd <= 0)
        return;

    const dstype eps = static_cast<dstype>(1.0e-12);
    for (Int d = 0; d < nd; ++d) {
        const Int nb = (nbins[d] > 0) ? nbins[d] : 1;
        const dstype span = gridMax[d] - gridMin[d];
        gridH[d] = (span > eps) ? span / static_cast<dstype>(nb) : static_cast<dstype>(1.0);
    }
}

/*
 * Compute nbins and grid spacing from a target isotropic spacing h0.
 *
 * Inputs:
 *   gridMin[nd], gridMax[nd], h0
 *
 * Outputs:
 *   nbins[nd], gridH[nd]
 */
static inline void ComputeEllipsoidGridResolution(
    Int* nbins,
    dstype* gridH,
    const dstype* gridMin,
    const dstype* gridMax,
    dstype h0,
    Int nd)
{
    if (nbins == nullptr || gridH == nullptr || gridMin == nullptr || gridMax == nullptr || nd <= 0)
        return;

    const dstype eps = static_cast<dstype>(1.0e-12);
    const dstype htarget = (h0 > eps) ? h0 : static_cast<dstype>(1.0);
    for (Int d = 0; d < nd; ++d) {
        const dstype span = gridMax[d] - gridMin[d];
        Int nb = 1;
        if (span > eps)
            nb = static_cast<Int>(std::ceil(span / htarget));
        if (nb < 1) nb = 1;
        nbins[d] = nb;
        gridH[d] = (span > eps) ? span / static_cast<dstype>(nb) : static_cast<dstype>(1.0);
    }
}

/*
 * Compute the inclusive integer bin range overlapped by one bounding box.
 *
 * Inputs:
 *   xmin[nd], xmax[nd], gridMin[nd], gridH[nd], nbins[nd]
 *
 * Outputs:
 *   ibeg[nd], iend[nd]
 */
static inline void ComputeOverlappedBinsForBox(
    Int* ibeg,
    Int* iend,
    const dstype* xmin,
    const dstype* xmax,
    const dstype* gridMin,
    const dstype* gridH,
    const Int* nbins,
    Int nd)
{
    if (ibeg == nullptr || iend == nullptr || xmin == nullptr || xmax == nullptr ||
        gridMin == nullptr || gridH == nullptr || nbins == nullptr || nd <= 0)
        return;

    for (Int d = 0; d < nd; ++d) {
        const Int nb = (nbins[d] > 0) ? nbins[d] : 1;
        const dstype h = (gridH[d] > static_cast<dstype>(0.0)) ? gridH[d] : static_cast<dstype>(1.0);
        Int i0 = static_cast<Int>(std::floor((xmin[d] - gridMin[d]) / h));
        Int i1 = static_cast<Int>(std::floor((xmax[d] - gridMin[d]) / h));
        if (i0 < 0) i0 = 0;
        if (i1 < 0) i1 = 0;
        if (i0 > nb - 1) i0 = nb - 1;
        if (i1 > nb - 1) i1 = nb - 1;
        if (i1 < i0) i1 = i0;
        ibeg[d] = i0;
        iend[d] = i1;
    }
}

/*
 * First pass of the CSR build: count how many ellipsoids overlap each bin.
 *
 * Inputs:
 *   ellipBoxes : [ne * 2*nd]
 *   gridMin[nd], gridH[nd], nbins[nd]
 *
 * Input/output:
 *   binCounts[nbinsTotal] : initialized to zero before the call
 */
static inline void CountEllipsoidGridOccupancy(
    Int* binCounts,
    const dstype* ellipBoxes,
    const dstype* gridMin,
    const dstype* gridH,
    const Int* nbins,
    Int nd,
    Int ne)
{
    if (binCounts == nullptr || ellipBoxes == nullptr || gridMin == nullptr ||
        gridH == nullptr || nbins == nullptr || nd <= 0 || ne <= 0)
        return;

    const Int boxStride = EllipsoidBoxStride(nd);
    std::vector<Int> ibeg(static_cast<size_t>(nd), 0);
    std::vector<Int> iend(static_cast<size_t>(nd), 0);
    std::vector<Int> idx(static_cast<size_t>(nd), 0);

    for (Int e = 0; e < ne; ++e) {
        const dstype* xmin = &ellipBoxes[boxStride * e];
        const dstype* xmax = &ellipBoxes[boxStride * e + nd];
        ComputeOverlappedBinsForBox(ibeg.data(), iend.data(), xmin, xmax, gridMin, gridH, nbins, nd);
        for (Int d = 0; d < nd; ++d) idx[d] = ibeg[d];

        while (true) {
            const Int flat = EllipsoidGridFlattenIndex(idx.data(), nbins, nd);
            if (flat >= 0) binCounts[flat] += 1;

            Int d = 0;
            for (; d < nd; ++d) {
                idx[d] += 1;
                if (idx[d] <= iend[d]) break;
                idx[d] = ibeg[d];
            }
            if (d == nd) break;
        }
    }
}

/*
 * Build CSR offsets from bin counts.
 *
 * Inputs:
 *   binCounts[nbinsTotal]
 *
 * Output:
 *   binOffsets[nbinsTotal + 1]
 */
static inline void BuildEllipsoidGridOffsets(
    Int* binOffsets,
    const Int* binCounts,
    Int nbinsTotal)
{
    if (binOffsets == nullptr || binCounts == nullptr || nbinsTotal < 0)
        return;

    binOffsets[0] = 0;
    for (Int b = 0; b < nbinsTotal; ++b)
        binOffsets[b + 1] = binOffsets[b] + binCounts[b];
}

/*
 * Second pass of the CSR build: fill the adjacency.
 *
 * Inputs:
 *   binOffsets[nbinsTotal + 1]
 *   ellipBoxes [ne * 2*nd]
 *   gridMin[nd], gridH[nd], nbins[nd]
 *
 * Input/output:
 *   binCursor[nbinsTotal] : initialize from binOffsets[0:nbinsTotal-1]
 *
 * Output:
 *   binElems[binOffsets[nbinsTotal]]
 */
static inline void FillEllipsoidGridAdjacency(
    Int* binElems,
    Int* binCursor,
    const Int* binOffsets,
    const dstype* ellipBoxes,
    const dstype* gridMin,
    const dstype* gridH,
    const Int* nbins,
    Int nd,
    Int ne)
{
    if (binElems == nullptr || binCursor == nullptr || binOffsets == nullptr ||
        ellipBoxes == nullptr || gridMin == nullptr || gridH == nullptr ||
        nbins == nullptr || nd <= 0 || ne <= 0)
        return;

    const Int nbinsTotal = EllipsoidGridTotalBins(nbins, nd);
    for (Int b = 0; b < nbinsTotal; ++b)
        binCursor[b] = binOffsets[b];

    const Int boxStride = EllipsoidBoxStride(nd);
    std::vector<Int> ibeg(static_cast<size_t>(nd), 0);
    std::vector<Int> iend(static_cast<size_t>(nd), 0);
    std::vector<Int> idx(static_cast<size_t>(nd), 0);

    for (Int e = 0; e < ne; ++e) {
        const dstype* xmin = &ellipBoxes[boxStride * e];
        const dstype* xmax = &ellipBoxes[boxStride * e + nd];
        ComputeOverlappedBinsForBox(ibeg.data(), iend.data(), xmin, xmax, gridMin, gridH, nbins, nd);
        for (Int d = 0; d < nd; ++d) idx[d] = ibeg[d];

        while (true) {
            const Int flat = EllipsoidGridFlattenIndex(idx.data(), nbins, nd);
            if (flat >= 0) {
                const Int pos = binCursor[flat];
                binElems[pos] = e;
                binCursor[flat] = pos + 1;
            }

            Int d = 0;
            for (; d < nd; ++d) {
                idx[d] += 1;
                if (idx[d] <= iend[d]) break;
                idx[d] = ibeg[d];
            }
            if (d == nd) break;
        }
    }
}

/*
 * Convenience builder for the full uniform-grid broad-phase index.
 *
 * Inputs:
 *   elemEllipsoids : [ne * (nd + nd*nd)]
 *   xdg            : [ne * npe * ncx]
 *   nbinsInput[nd]
 *   tm             : scratch of size GetEllipsoidBoundingBoxScratchSize(nd)
 *   itmp           : scratch of size GetEllipsoidUniformGridIntScratchSize(nbinsInput, nd)
 *
 * Outputs:
 *   ellipBoxes : [ne * 2*nd]
 *   gridMin[nd], gridH[nd], nbins[nd]
 *   binOffsets[nbinsTotal + 1]
 *   binElems[binOffsets[nbinsTotal]]
 */
static inline void BuildEllipsoidUniformGrid(
    dstype* ellipBoxes,
    dstype* gridMin,
    dstype* gridH,
    Int* nbins,
    Int* binOffsets,
    Int* binElems,
    const dstype* elemEllipsoids,
    const dstype* xdg,
    Int nd,
    Int npe,
    Int ncx,
    Int ne,
    const Int* nbinsInput,
    dstype* tm,
    Int* itmp)
{
    if (ellipBoxes == nullptr || gridMin == nullptr || gridH == nullptr || nbins == nullptr ||
        binOffsets == nullptr || elemEllipsoids == nullptr || xdg == nullptr ||
        nbinsInput == nullptr || tm == nullptr || itmp == nullptr ||
        nd <= 0 || npe <= 0 || ncx <= 0 || ne <= 0)
        return;

    for (Int d = 0; d < nd; ++d)
        nbins[d] = (nbinsInput[d] > 0) ? nbinsInput[d] : 1;

    std::vector<dstype> gridMax(static_cast<size_t>(nd), 0.0);
    BuildEllipsoidBoundingBoxes(ellipBoxes, elemEllipsoids, nd, ne, tm);
    ComputeDomainGridBounds(gridMin, gridMax.data(), xdg, nd, npe, ncx, ne);
    ComputeEllipsoidGridSpacing(gridH, gridMin, gridMax.data(), nbins, nd);

    const Int nbinsTotal = EllipsoidGridTotalBins(nbins, nd);
    Int* binCounts = &itmp[0];
    Int* binCursor = &itmp[nbinsTotal];
    for (Int b = 0; b < nbinsTotal; ++b)
        binCounts[b] = 0;

    CountEllipsoidGridOccupancy(binCounts, ellipBoxes, gridMin, gridH, nbins, nd, ne);
    BuildEllipsoidGridOffsets(binOffsets, binCounts, nbinsTotal);

    if (binElems != nullptr)
        FillEllipsoidGridAdjacency(binElems, binCursor, binOffsets, ellipBoxes, gridMin, gridH, nbins, nd, ne);
}

/*
 * Locate the grid bin containing xphys.
 *
 * Inputs:
 *   xphys[nd], gridMin[nd], gridH[nd], nbins[nd]
 *
 * Outputs:
 *   ibin[nd]   : per-dimension bin index
 *   flatBin    : flattened bin index
 *
 * Returns false if xphys lies outside the grid domain.
 */
static inline bool FindEllipsoidGridBin(
    Int* ibin,
    Int* flatBin,
    const dstype* xphys,
    const dstype* gridMin,
    const dstype* gridH,
    const Int* nbins,
    Int nd)
{
    if (xphys == nullptr || gridMin == nullptr || gridH == nullptr || nbins == nullptr || nd <= 0)
        return false;
    
    Int flat = 0;
    Int stride = 1;

    for (Int d = 0; d < nd; ++d) {
        const Int nb = (nbins[d] > 0) ? nbins[d] : 1;
        const dstype h = (gridH[d] > static_cast<dstype>(0.0)) ? gridH[d] : static_cast<dstype>(1.0);
        const Int i = static_cast<Int>(std::floor((xphys[d] - gridMin[d]) / h));
        if (i < 0 || i >= nb)
            return false;

        if (ibin != nullptr)
            ibin[d] = i;
        flat += i * stride;
        stride *= nb;
    }

    if (flatBin != nullptr)
        *flatBin = flat;
    return true;
}

/*
 * Broad-phase candidate query using only the bin adjacency.
 *
 * Inputs:
 *   xphys[nd]
 *   binOffsets[nbinsTotal + 1], binElems[nbinElems]
 *   gridMin[nd], gridH[nd], nbins[nd]
 *
 * Output:
 *   candidateElems[maxCandidates]
 *
 * Returns the number of broad-phase candidates copied into candidateElems.
 */
static inline Int BuildCandidateElementsFromEllipsoidGrid(
    Int* candidateElems,
    Int maxCandidates,
    const dstype* xphys,
    const Int* binOffsets,
    const Int* binElems,
    const dstype* gridMin,
    const dstype* gridH,
    const Int* nbins,
    Int nd)
{
    if (candidateElems == nullptr || xphys == nullptr || binOffsets == nullptr ||
        binElems == nullptr || gridMin == nullptr || gridH == nullptr ||
        nbins == nullptr || maxCandidates <= 0)
        return 0;

    Int flatBin = -1;
    if (!FindEllipsoidGridBin(nullptr, &flatBin, xphys, gridMin, gridH, nbins, nd))
        return 0;

    Int ncandidates = 0;
    const Int start = binOffsets[flatBin];
    const Int end = binOffsets[flatBin + 1];
    for (Int i = start; i < end; ++i) {
        ncandidates = AppendUniqueCandidate(candidateElems, ncandidates, maxCandidates, binElems[i]);
        if (ncandidates >= maxCandidates) break;
    }
    return ncandidates;
}

/*
 * Broad-phase + exact ellipsoid query.
 *
 * This is the intended fast replacement for the full O(ne) scan:
 *   1. find the bin containing xphys
 *   2. fetch ellipsoid ids in that bin
 *   3. run PointInEllipsoid(...) only on those ids
 *
 * Inputs:
 *   xphys[nd]
 *   binOffsets[nbinsTotal + 1], binElems[nbinElems]
 *   gridMin[nd], gridH[nd], nbins[nd]
 *   elemEllipsoids[ne * (nd + nd*nd)]
 *
 * Output:
 *   candidateElems[maxCandidates]
 *
 * Returns the number of exact candidates copied into candidateElems.
 */
static inline Int BuildCandidateElementsFromEllipsoidGridExact(
    Int* candidateElems,
    Int maxCandidates,
    const dstype* xphys,
    const Int* binOffsets,
    const Int* binElems,
    const dstype* gridMin,
    const dstype* gridH,
    const Int* nbins,
    const dstype* elemEllipsoids,
    Int nd,
    dstype tol)
{
    if (candidateElems == nullptr || xphys == nullptr || binOffsets == nullptr ||
        binElems == nullptr || gridMin == nullptr || gridH == nullptr || nbins == nullptr ||
        elemEllipsoids == nullptr || maxCandidates <= 0)
        return 0;

    Int flatBin = -1;
    if (!FindEllipsoidGridBin(nullptr, &flatBin, xphys, gridMin, gridH, nbins, nd))
        return 0;

    const Int ellStride = EllipsoidStride(nd);
    Int ncandidates = 0;
    const Int start = binOffsets[flatBin];
    const Int end = binOffsets[flatBin + 1];
    for (Int i = start; i < end; ++i) {
        const Int e = binElems[i];
        const dstype* center = &elemEllipsoids[ellStride * e];
        const dstype* Aell = &elemEllipsoids[ellStride * e + nd];
        if (PointInEllipsoid(xphys, center, Aell, nd, tol)) {
            ncandidates = AppendUniqueCandidate(candidateElems, ncandidates, maxCandidates, e);
            if (ncandidates >= maxCandidates) break;
        }
    }
    return ncandidates;
}

#endif
