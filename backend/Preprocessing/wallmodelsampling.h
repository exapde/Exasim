#pragma once

#include <vector>

#include "../Common/common.h"

class CDiscretization;

namespace exasim {
namespace wm {

struct WallModelSamplingOptions {
    Int ibc = -1;
    dstype y1 = 0.0;
    Int maxNewtonIter = 20;
    dstype newtonTol = 1.0e-12;
    dstype insideTol = 1.0e-10;
    bool useAdjacentFaceElemFirst = true;
    bool allowNeighborSearch = true;
    bool allowGlobalSearch = false;
};

struct WallModelSamplingData {
    Int ibc = -1;
    Int nd = 0;
    Int ncx = 0;
    Int npe = 0;
    Int npf = 0;
    Int ngf = 0;
    Int nfaces = 0;
    Int npoints = 0;
    dstype y1 = 0.0;

    std::vector<Int> faces;        // size = nfaces
    std::vector<dstype> xw;        // size = nd * npoints
    std::vector<dstype> nw;        // size = nd * npoints
    std::vector<dstype> x1;        // size = nd * npoints
    std::vector<Int> e1;           // size = npoints
    std::vector<dstype> xi1;       // size = nd * npoints
    std::vector<dstype> shap1;     // size = npe * npoints
    std::vector<dstype> x1_check;  // size = nd * npoints
};

inline Int PointIndex(const Int ig, const Int iface, const Int ngf)
{
    return ig + ngf * iface;
}

inline Int VecIndex(const Int d, const Int ip, const Int nd)
{
    return d + nd * ip;
}

inline Int BasisIndex(const Int a, const Int ip, const Int npe)
{
    return a + npe * ip;
}

void ComputeOffWallPoints(
    std::vector<dstype>& x1,
    const std::vector<dstype>& xw,
    const std::vector<dstype>& nw,
    dstype y1,
    Int nd,
    Int npoints);

void ComputeOffWallPoints(
    dstype* x1,
    const dstype* xw,
    const dstype* nw,
    dstype y1,
    Int nd,
    Int npoints);

Int GetInteriorElementAdjacentToBoundaryFace(
    const meshstruct& mesh,
    Int face);

void ExtractElementNodes(
    std::vector<dstype>& xdg_elem,
    const dstype* xdg,
    Int elem,
    Int npe,
    Int ncx,
    Int nd);

bool PointInReferenceElement(
    const dstype* xi,
    Int nd,
    Int elemtype,
    dstype tol);

void EvaluateVolumeShapeAndGradientsAtXi(
    dstype* N,
    dstype* dN,
    const dstype* xi,
    const masterstruct& master,
    Int nd,
    Int elemtype,
    Int porder,
    Int npe);

void MapReferenceToPhysical(
    dstype* x,
    dstype* J,
    const dstype* xi,
    const dstype* xdg_elem,
    const masterstruct& master,
    Int nd,
    Int elemtype,
    Int porder,
    Int npe);

bool ComputeReferenceCoordinatesInElement(
    dstype* xi,
    const dstype* xphys,
    const dstype* xdg_elem,
    const masterstruct& master,
    Int nd,
    Int elemtype,
    Int porder,
    Int npe,
    Int maxNewtonIter,
    dstype tol);

void EvaluateVolumeBasisAtXi(
    dstype* N,
    const dstype* xi,
    const masterstruct& master,
    Int nd,
    Int elemtype,
    Int porder,
    Int npe);

void EvaluateVolumeBasisBatch(
    std::vector<dstype>& shap1,
    const std::vector<dstype>& xi1,
    const masterstruct& master,
    Int nd,
    Int elemtype,
    Int porder,
    Int npe,
    Int npoints);

bool FindContainingElementForPoint(
    Int& elemOut,
    dstype* xiOut,
    const dstype* xphys,
    Int face,
    const dstype* xdg,
    const meshstruct& mesh,
    const masterstruct& master,
    const WallModelSamplingOptions& opts,
    Int nd,
    Int ncx,
    Int npe,
    Int ne,
    Int elemtype,
    Int porder);

void FindContainingElements(
    std::vector<Int>& e1,
    std::vector<dstype>& xi1,
    const std::vector<dstype>& x1,
    const std::vector<Int>& faces,
    const dstype* xdg,
    const meshstruct& mesh,
    const masterstruct& master,
    const WallModelSamplingOptions& opts,
    Int nd,
    Int ncx,
    Int npe,
    Int ne,
    Int ngf,
    Int elemtype,
    Int porder);

void GatherWallGaussPointsAndNormals(
    WallModelSamplingData& wm,
    CDiscretization& disc,
    Int ibc);

void BuildWallModelSamplingData(
    WallModelSamplingData& wm,
    CDiscretization& disc,
    const WallModelSamplingOptions& opts);

void ValidateWallModelSamplingData(
    WallModelSamplingData& wm,
    const dstype* xdg,
    const masterstruct& master,
    Int ncx,
    Int elemtype,
    Int porder,
    dstype tol = 1.0e-10);

void EvaluateOffWallStateAtPoint(
    dstype* U1p,
    const dstype* udg,
    Int elem,
    const dstype* N,
    Int npe,
    Int nc);

void EvaluateOffWallState(
    dstype* U1,
    const dstype* udg,
    const Int* e1,
    const dstype* shap1,
    Int npe,
    Int nc,
    Int npoints);

void ConservativeToPrimitiveAtOffWallPoint(
    dstype& rho1,
    dstype* vel1,
    dstype& p1,
    dstype& T1,
    const dstype* U1,
    const dstype* param,
    Int nd,
    Int nc);

void EvaluateOffWallStateForBoundaryBlock(
    dstype* U1,
    const dstype* udg,
    const WallModelSamplingData& wm,
    Int npe,
    Int nc);

dstype ComputeTangentialSpeed(
    dstype* that,
    const dstype* vel1,
    const dstype* n,
    Int nd);

dstype ComputeBFWMInput(
    dstype y1,
    dstype utmag1,
    dstype rho1,
    dstype mu1,
    dstype k1,
    dstype cp);

} // namespace wm
} // namespace exasim
