#ifndef __CPOINTLOCATOR_H__
#define __CPOINTLOCATOR_H__

#include "../Common/common.h"
#include <vector>

//class CDiscretization;

struct WallModelSamplingData {
    Int ibc = -1;
    Int nd = 0;
    Int ncx = 0;
    Int npe = 0;
    Int npf = 0;
    Int ngf = 0;
    Int nfe = 0;
    Int nfaces = 0;
    Int npoints = 0;
    Int nbe1 = 0;
    dstype y1 = 0.0;

    std::vector<Int> faces;        // size = nfaces
    std::vector<Int> nextfaces;    // size = nbe1+1
    std::vector<Int> elems;        // size = npoints, owner element of each point in xw
    std::vector<dstype> xw;        // size = nd * npoints
    std::vector<dstype> nw;        // size = nd * npoints
    std::vector<dstype> x1;        // size = nd * npoints
    std::vector<dstype> xi1;       // size = nd * npoints
    std::vector<dstype> shap1;     // size = npe * npoints
};

class CPointLocator {
public:
    WallModelSamplingData wm;
  
    CPointLocator();

    CPointLocator(
        const dstype* xdg,
        const dstype* xpe,
        Int nd,
        Int npe,
        Int ncx,
        Int elemtype,
        Int porder);

    void SetGeometry(
        const dstype* xdg,
        const dstype* xpe,
        Int nd,
        Int npe,
        Int ncx,
        Int elemtype,
        Int porder);

    void SetElementFaceConnectivity(
        const Int* e2f,
        const Int* f2e,
        Int nfe);

    void SetEllipsoidGrid(
        const Int* binOffsets,
        const Int* binElems,
        const dstype* gridMin,
        const dstype* gridH,
        const Int* nbins,
        const dstype* elemEllipsoids,
        dstype ellipTol);

    void SetField(
        const dstype* U,
        Int nc);

    void SetSearchParameters(
        Int maxNewtonIter,
        dstype newtonTol,
        dstype insideTol);

    bool BuildWallModelSamplingData(
        CDiscretization& disc,
        Int ibc,
        dstype y1);

    bool FindPointAndShapeFunctionsFromElementFaces(
        Int* elem,
        dstype* xi,
        dstype* N,
        const dstype* xphys,
        Int seedElem,
        Int* candidateElems,
        Int maxCandidates) const;

    bool FindPointShapeAndFieldFromEllipsoidGrid(
        Int* elem,
        dstype* xi,
        dstype* N,
        dstype* u,
        const dstype* xphys,
        Int* candidateElems,
        Int maxCandidates) const;

private:
    const dstype* m_xdg;
    const dstype* m_xpe;
    const dstype* m_U;

    const Int* m_e2f;
    const Int* m_f2e;

    const Int* m_binOffsets;
    const Int* m_binElems;
    const Int* m_nbins;
    const dstype* m_gridMin;
    const dstype* m_gridH;
    const dstype* m_elemEllipsoids;

    Int m_nd;
    Int m_npe;
    Int m_ncx;
    Int m_nc;
    Int m_nfe;
    Int m_elemtype;
    Int m_porder;
    Int m_maxNewtonIter;

    dstype m_newtonTol;
    dstype m_insideTol;
    dstype m_ellipTol;
};

#endif
