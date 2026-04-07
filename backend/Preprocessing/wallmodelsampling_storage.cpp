#include "wallmodelsampling.h"

#include "../Discretization/discretization.h"

#include <cmath>
#include <sstream>

namespace exasim {
namespace wm {

namespace {

void CopyFieldToHost(std::vector<dstype>& host, dstype* data, Int n, Int backend)
{
    host.resize(static_cast<size_t>(n));
    TemplateCopytoHost(host.data(), data, n, backend);
}

void CopyVolumeCoordinatesToHost(std::vector<dstype>& xdg_host, CDiscretization& disc)
{
    const Int n = disc.common.npe * disc.common.ncx * disc.common.ne;
    xdg_host.resize(static_cast<size_t>(n));
    TemplateCopytoHost(xdg_host.data(), disc.sol.xdg, n, disc.common.backend);
}

} // namespace

void GatherWallGaussPointsAndNormals(
    WallModelSamplingData& wm,
    CDiscretization& disc,
    const Int ibc)
{
    wm.ibc = ibc;
    wm.nd = disc.common.nd;
    wm.ncx = disc.common.ncx;
    wm.npe = disc.common.npe;
    wm.npf = disc.common.npf;
    wm.ngf = disc.common.ngf;

    Int* faces_raw = nullptr;
    wm.nfaces = disc.getFacesOnInterface(&faces_raw, ibc);
    wm.npoints = wm.ngf * wm.nfaces;

    if (wm.nfaces <= 0) {
        std::ostringstream oss;
        oss << "GatherWallGaussPointsAndNormals found no faces for boundary condition " << ibc;
        error(oss.str());
    }

    wm.faces.assign(faces_raw, faces_raw + wm.nfaces);
    CPUFREE(faces_raw);

    dstype* xdgint = nullptr;
    dstype* nlint = nullptr;
    dstype* xdggint = nullptr;
    dstype* nlgint = nullptr;

    const Int backend = disc.common.backend;
    TemplateMalloc(&xdgint, wm.npf * wm.nfaces * wm.ncx, backend);
    TemplateMalloc(&nlint, wm.npf * wm.nfaces * wm.nd, backend);
    TemplateMalloc(&xdggint, wm.ngf * wm.nfaces * wm.ncx, backend);
    TemplateMalloc(&nlgint, wm.ngf * wm.nfaces * wm.nd, backend);

    disc.getDGNodesOnInterface(xdgint, wm.faces.data(), wm.nfaces);
    disc.getNormalVectorOnInterface(nlint, xdgint, wm.nfaces);
    disc.getFieldsAtGaussPointsOnInterface(xdggint, xdgint, wm.nfaces, wm.ncx);
    disc.getFieldsAtGaussPointsOnInterface(nlgint, nlint, wm.nfaces, wm.nd);

    std::vector<dstype> xg_full;
    std::vector<dstype> ng_full;
    CopyFieldToHost(xg_full, xdggint, wm.npoints * wm.ncx, backend);
    CopyFieldToHost(ng_full, nlgint, wm.npoints * wm.nd, backend);

    wm.xw.resize(static_cast<size_t>(wm.npoints * wm.nd));
    wm.nw.resize(static_cast<size_t>(wm.npoints * wm.nd));
    for (Int ip = 0; ip < wm.npoints; ++ip) {
        for (Int d = 0; d < wm.nd; ++d) {
            wm.xw[VecIndex(d, ip, wm.nd)] = xg_full[d + wm.ncx * ip];
            wm.nw[VecIndex(d, ip, wm.nd)] = ng_full[d + wm.nd * ip];
        }
    }

    TemplateFree(xdgint, backend);
    TemplateFree(nlint, backend);
    TemplateFree(xdggint, backend);
    TemplateFree(nlgint, backend);
}

void BuildWallModelSamplingData(
    WallModelSamplingData& wm,
    CDiscretization& disc,
    const WallModelSamplingOptions& opts)
{
    GatherWallGaussPointsAndNormals(wm, disc, opts.ibc);
    wm.y1 = opts.y1;

    ComputeOffWallPoints(wm.x1, wm.xw, wm.nw, opts.y1, wm.nd, wm.npoints);

    std::vector<dstype> xdg_host;
    CopyVolumeCoordinatesToHost(xdg_host, disc);

    std::vector<dstype> xpe_host(static_cast<size_t>(disc.master.szxpe));
    TemplateCopytoHost(xpe_host.data(), disc.master.xpe, disc.master.szxpe, disc.common.backend);
    masterstruct master_host = disc.master;
    master_host.xpe = xpe_host.data();

    FindContainingElements(
        wm.e1, wm.xi1, wm.x1, wm.faces, xdg_host.data(),
        disc.mesh, master_host, opts, wm.nd, wm.ncx, wm.npe,
        disc.common.ne, wm.ngf, disc.common.elemtype, disc.common.porder);

    EvaluateVolumeBasisBatch(
        wm.shap1, wm.xi1, master_host, wm.nd, disc.common.elemtype,
        disc.common.porder, wm.npe, wm.npoints);

    ValidateWallModelSamplingData(
        wm, xdg_host.data(), master_host, wm.ncx, disc.common.elemtype,
        disc.common.porder, opts.insideTol);
}

void ValidateWallModelSamplingData(
    WallModelSamplingData& wm,
    const dstype* xdg,
    const masterstruct& master,
    const Int ncx,
    const Int elemtype,
    const Int porder,
    const dstype tol)
{
    wm.x1_check.resize(static_cast<size_t>(wm.nd * wm.npoints));

    for (Int ip = 0; ip < wm.npoints; ++ip) {
        std::vector<dstype> xdg_elem;
        ExtractElementNodes(xdg_elem, xdg, wm.e1[ip], wm.npe, ncx, wm.nd);

        dstype xchk[3] = {0.0, 0.0, 0.0};
        dstype J[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        MapReferenceToPhysical(
            xchk, J, &wm.xi1[wm.nd * ip], xdg_elem.data(), master,
            wm.nd, elemtype, porder, wm.npe);

        dstype err2 = 0.0;
        dstype shape_sum = 0.0;
        for (Int d = 0; d < wm.nd; ++d) {
            wm.x1_check[VecIndex(d, ip, wm.nd)] = xchk[d];
            const dstype diff = xchk[d] - wm.x1[VecIndex(d, ip, wm.nd)];
            err2 += diff * diff;
        }
        for (Int a = 0; a < wm.npe; ++a)
            shape_sum += wm.shap1[BasisIndex(a, ip, wm.npe)];

        if (std::sqrt(err2) > 10.0 * tol || std::abs(shape_sum - 1.0) > 100.0 * tol) {
            std::ostringstream oss;
            oss << "ValidateWallModelSamplingData failed at ip=" << ip
                << ", reconstruction error=" << std::sqrt(err2)
                << ", basis sum error=" << std::abs(shape_sum - 1.0);
            error(oss.str());
        }
    }
}

} // namespace wm
} // namespace exasim
