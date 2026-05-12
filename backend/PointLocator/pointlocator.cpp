#include "pointlocator.h"
#include "pointlocation.cpp"
#include "pointwallmodel.cpp"

CPointLocator::CPointLocator()
    : m_xdg(nullptr),
      m_xpe(nullptr),
      m_U(nullptr),
      m_e2f(nullptr),
      m_f2e(nullptr),
      m_binOffsets(nullptr),
      m_binElems(nullptr),
      m_nbins(nullptr),
      m_gridMin(nullptr),
      m_gridH(nullptr),
      m_elemEllipsoids(nullptr),
      m_nd(0),
      m_npe(0),
      m_ncx(0),
      m_nc(0),
      m_nfe(0),
      m_elemtype(0),
      m_porder(0),
      m_maxNewtonIter(20),
      m_newtonTol(static_cast<dstype>(1.0e-8)),
      m_insideTol(static_cast<dstype>(1.0e-6)),
      m_ellipTol(static_cast<dstype>(1.0e-6))
{
}

CPointLocator::CPointLocator(
    const dstype* xdg,
    const dstype* xpe,
    Int nd,
    Int npe,
    Int ncx,
    Int elemtype,
    Int porder)
    : CPointLocator()
{
    SetGeometry(xdg, xpe, nd, npe, ncx, elemtype, porder);
}

void CPointLocator::SetGeometry(
    const dstype* xdg,
    const dstype* xpe,
    Int nd,
    Int npe,
    Int ncx,
    Int elemtype,
    Int porder)
{
    m_xdg = xdg;
    m_xpe = xpe;
    m_nd = nd;
    m_npe = npe;
    m_ncx = ncx;
    m_elemtype = elemtype;
    m_porder = porder;
}

void CPointLocator::SetElementFaceConnectivity(
    const Int* e2f,
    const Int* f2e,
    Int nfe)
{
    m_e2f = e2f;
    m_f2e = f2e;
    m_nfe = nfe;
}

void CPointLocator::SetEllipsoidGrid(
    const Int* binOffsets,
    const Int* binElems,
    const dstype* gridMin,
    const dstype* gridH,
    const Int* nbins,
    const dstype* elemEllipsoids,
    dstype ellipTol)
{
    m_binOffsets = binOffsets;
    m_binElems = binElems;
    m_gridMin = gridMin;
    m_gridH = gridH;
    m_nbins = nbins;
    m_elemEllipsoids = elemEllipsoids;
    m_ellipTol = ellipTol;
}

void CPointLocator::SetField(
    const dstype* U,
    Int nc)
{
    m_U = U;
    m_nc = nc;
}

void CPointLocator::SetSearchParameters(
    Int maxNewtonIter,
    dstype newtonTol,
    dstype insideTol)
{
    m_maxNewtonIter = maxNewtonIter;
    m_newtonTol = newtonTol;
    m_insideTol = insideTol;
}

bool CPointLocator::BuildWallModelSamplingData(
    CDiscretization& disc,
    Int ibc,
    dstype y1)
{
    SetGeometry(
        disc.sol.xdg, disc.master.xpe, disc.common.nd, disc.common.npe,
        disc.common.ncx, disc.common.elemtype, disc.common.porder);
    SetElementFaceConnectivity(disc.mesh.e2f, disc.mesh.f2e, disc.common.nfe);

    return ::BuildWallModelSamplingData(
        wm, disc, ibc, y1, m_maxNewtonIter, m_newtonTol, m_insideTol);
}

bool CPointLocator::FindPointAndShapeFunctionsFromElementFaces(
    Int* elem,
    dstype* xi,
    dstype* N,
    const dstype* xphys,
    Int seedElem,
    Int* candidateElems,
    Int maxCandidates) const
{
    if (m_xdg == nullptr || m_xpe == nullptr || m_e2f == nullptr || m_f2e == nullptr)
        return false;

    return ::FindPointAndShapeFunctionsFromElementFaces(
        elem, xi, N, xphys, m_xdg, seedElem, candidateElems, maxCandidates,
        m_e2f, m_f2e, m_nfe, m_xpe, m_nd, m_npe, m_ncx, m_elemtype, m_porder,
        m_maxNewtonIter, m_newtonTol, m_insideTol);
}

bool CPointLocator::FindPointShapeAndFieldFromEllipsoidGrid(
    Int* elem,
    dstype* xi,
    dstype* N,
    dstype* u,
    const dstype* xphys,
    Int* candidateElems,
    Int maxCandidates) const
{
    if (m_xdg == nullptr || m_xpe == nullptr || m_U == nullptr ||
        m_binOffsets == nullptr || m_binElems == nullptr ||
        m_gridMin == nullptr || m_gridH == nullptr || m_nbins == nullptr ||
        m_elemEllipsoids == nullptr)
        return false;

    return ::FindPointShapeAndFieldFromEllipsoidGrid(
        elem, xi, N, u, xphys, m_U, m_xdg, candidateElems, maxCandidates,
        m_binOffsets, m_binElems, m_gridMin, m_gridH, m_nbins, m_elemEllipsoids,
        m_xpe, m_nd, m_npe, m_nc, m_ncx, m_elemtype, m_porder,
        m_maxNewtonIter, m_newtonTol, m_insideTol, m_ellipTol);
}
