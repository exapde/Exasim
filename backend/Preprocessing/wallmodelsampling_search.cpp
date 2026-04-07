#include "wallmodelsampling.h"

#include <algorithm>
#include <limits>
#include <sstream>

namespace exasim {
namespace wm {

Int GetInteriorElementAdjacentToBoundaryFace(
    const meshstruct& mesh,
    const Int face)
{
    if (mesh.f2e == nullptr) return -1;
    return mesh.f2e[4 * face + 0];
}

void ExtractElementNodes(
    std::vector<dstype>& xdg_elem,
    const dstype* xdg,
    const Int elem,
    const Int npe,
    const Int ncx,
    const Int nd)
{
    xdg_elem.resize(static_cast<size_t>(npe * nd));
    const Int offset = npe * ncx * elem;
    for (Int d = 0; d < nd; ++d)
        for (Int a = 0; a < npe; ++a)
            xdg_elem[a + npe * d] = xdg[offset + a + npe * d];
}

namespace {

bool TryCandidateElement(
    Int elem,
    dstype* xiOut,
    const dstype* xphys,
    const dstype* xdg,
    const masterstruct& master,
    const WallModelSamplingOptions& opts,
    Int nd,
    Int ncx,
    Int npe,
    Int elemtype,
    Int porder)
{
    if (elem < 0) return false;

    std::vector<dstype> xdg_elem;
    ExtractElementNodes(xdg_elem, xdg, elem, npe, ncx, nd);

    if (!ComputeReferenceCoordinatesInElement(
            xiOut, xphys, xdg_elem.data(), master, nd, elemtype, porder, npe,
            opts.maxNewtonIter, opts.newtonTol))
        return false;

    return PointInReferenceElement(xiOut, nd, elemtype, opts.insideTol);
}

} // namespace

bool FindContainingElementForPoint(
    Int& elemOut,
    dstype* xiOut,
    const dstype* xphys,
    const Int face,
    const dstype* xdg,
    const meshstruct& mesh,
    const masterstruct& master,
    const WallModelSamplingOptions& opts,
    const Int nd,
    const Int ncx,
    const Int npe,
    const Int ne,
    const Int elemtype,
    const Int porder)
{
    elemOut = -1;
    std::fill(xiOut, xiOut + nd, std::numeric_limits<dstype>::quiet_NaN());

    if (opts.useAdjacentFaceElemFirst) {
        const Int elem = GetInteriorElementAdjacentToBoundaryFace(mesh, face);
        if (TryCandidateElement(elem, xiOut, xphys, xdg, master, opts, nd, ncx, npe, elemtype, porder)) {
            elemOut = elem;
            return true;
        }
    }

    if (!(opts.allowNeighborSearch || opts.allowGlobalSearch))
        return false;

    for (Int elem = 0; elem < ne; ++elem) {
        if (TryCandidateElement(elem, xiOut, xphys, xdg, master, opts, nd, ncx, npe, elemtype, porder)) {
            elemOut = elem;
            return true;
        }
    }

    return false;
}

void FindContainingElements(
    std::vector<Int>& e1,
    std::vector<dstype>& xi1,
    const std::vector<dstype>& x1,
    const std::vector<Int>& faces,
    const dstype* xdg,
    const meshstruct& mesh,
    const masterstruct& master,
    const WallModelSamplingOptions& opts,
    const Int nd,
    const Int ncx,
    const Int npe,
    const Int ne,
    const Int ngf,
    const Int elemtype,
    const Int porder)
{
    const Int nfaces = static_cast<Int>(faces.size());
    const Int npoints = ngf * nfaces;

    e1.assign(static_cast<size_t>(npoints), -1);
    xi1.resize(static_cast<size_t>(nd * npoints));

    for (Int iface = 0; iface < nfaces; ++iface) {
        const Int face = faces[iface];
        for (Int ig = 0; ig < ngf; ++ig) {
            const Int ip = PointIndex(ig, iface, ngf);
            dstype* xi = &xi1[nd * ip];
            Int elem = -1;
            const bool found = FindContainingElementForPoint(
                elem, xi, &x1[nd * ip], face, xdg, mesh, master, opts,
                nd, ncx, npe, ne, elemtype, porder);

            if (!found || elem < 0) {
                std::ostringstream oss;
                oss << "FindContainingElements failed for wall sample point ip=" << ip
                    << " on face=" << face;
                error(oss.str());
            }
            e1[ip] = elem;
        }
    }
}

} // namespace wm
} // namespace exasim
