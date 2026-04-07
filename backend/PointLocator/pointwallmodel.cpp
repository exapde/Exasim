#ifndef __POINTWALLMODEL_H__
#define __POINTWALLMODEL_H__

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

static inline void cpuFaceGeom1D(dstype* jacg, dstype* nlg, const dstype* Jg, const Int na)
{
    (void) Jg;
    for (Int i = 0; i < na; ++i) {
        jacg[i] = 1.0;
        nlg[i] = 1.0;
    }
}

static inline void cpuFaceGeom2D(dstype* jacg, dstype* nlg, const dstype* Jg, const Int na)
{
    for (Int i = 0; i < na; ++i) {
        const Int j = i + na;
        jacg[i] = std::sqrt(Jg[i] * Jg[i] + Jg[j] * Jg[j]);
        nlg[i] = Jg[j] / jacg[i];
        nlg[j] = -Jg[i] / jacg[i];
    }
}

static inline void cpuFaceGeom3D(dstype* jacg, dstype* nlg, const dstype* Jg, const Int na)
{
    const Int n11 = 0;
    const Int n21 = na;
    const Int n31 = 2 * na;
    const Int n12 = 3 * na;
    const Int n22 = 4 * na;
    const Int n32 = 5 * na;

    for (Int i = 0; i < na; ++i) {
        const Int j = i + na;
        const Int k = i + 2 * na;
        nlg[i] = Jg[i + n21] * Jg[i + n32] - Jg[i + n31] * Jg[i + n22];
        nlg[j] = Jg[i + n31] * Jg[i + n12] - Jg[i + n11] * Jg[i + n32];
        nlg[k] = Jg[i + n11] * Jg[i + n22] - Jg[i + n21] * Jg[i + n12];
        jacg[i] = std::sqrt(nlg[i] * nlg[i] + nlg[j] * nlg[j] + nlg[k] * nlg[k]);
        nlg[i] /= jacg[i];
        nlg[j] /= jacg[i];
        nlg[k] /= jacg[i];
    }
}

static inline void ComputeOffWallPoints(
    dstype* x1,
    const dstype* xw,
    const dstype* nw,
    const dstype y1,
    const Int nd,
    const Int npoints)
{
    if (x1 == nullptr || xw == nullptr || nw == nullptr || nd <= 0 || npoints <= 0)
        return;

    for (Int ip = 0; ip < npoints; ++ip) {
        for (Int d = 0; d < nd; ++d) {
            const Int k = VecIndex(d, ip, nd);
            x1[k] = xw[k] - y1 * nw[k];
        }
    }
}

static inline Int getFacesOnInterface(
    std::vector<Int>& faces,
    std::vector<Int>& nextfaces,
    CDiscretization& disc,
    const Int boundarycondition)
{
    const Int nintfaces = getinterfacefaces(
        disc.mesh.bf, disc.common.nfe, disc.common.ne1, boundarycondition);

    faces.resize(static_cast<size_t>(nintfaces));
    getinterfacefaces(
        faces.data(), disc.mesh.bf, disc.common.nfe, disc.common.ne1,
        boundarycondition, nintfaces);

    nextfaces.resize(static_cast<size_t>(disc.common.nbe1 + 1));
    Int nfacestotal = 0;
    nextfaces[0] = nfacestotal;
    for (Int j = 0; j < disc.common.nbe1; ++j) {
        const Int e1 = disc.common.eblks[3 * j] - 1;
        const Int e2 = disc.common.eblks[3 * j + 1];
        const Int nfe = disc.common.nfe;
        const Int nfaces = getinterfacefaces(
            &disc.mesh.bf[nfe * e1], nfe, e2 - e1, boundarycondition);
        nfacestotal += nfaces;
        nextfaces[j + 1] = nfacestotal;
    }

    if (nfacestotal != nintfaces)
        error("getFacesOnInterface is wrong because the numbers of faces do not match.");

    return nintfaces;
}

static inline void getNodesOnInterface(
    dstype* ub,
    const dstype* udg,
    const Int* boufaces,
    const Int* perm,
    const Int nfe,
    const Int npf,
    const Int npe,
    const Int ncu,
    const Int nc,
    const Int nf)
{
    const Int K = npf * nf;
    const Int N = K * ncu;
    const Int M = npe * nc;
    for (Int idx = 0; idx < N; ++idx) {
        const Int i = idx / K;
        const Int j = idx % K;
        const Int n = j % npf;
        const Int k = j / npf;
        const Int e1 = boufaces[k] / nfe;
        const Int l1 = boufaces[k] % nfe;
        const Int m = perm[n + npf * l1];
        ub[idx] = udg[m + i * npe + e1 * M];
    }
}

static inline void getFieldsAtGaussPointsOnInterface(
    dstype* ugint,
    dstype* uint,
    CDiscretization& disc,
    const Int nfaces,
    const Int ncomp)
{
    Node2Gauss(
        disc.common.cublasHandle, ugint, uint, disc.master.shapfgt,
        disc.common.ngf, disc.common.npf, nfaces * ncomp, 0);
}

static inline void getNormalVectorOnInterface(
    dstype* nlint,
    dstype* xdgint,
    CDiscretization& disc,
    const Int nfaces)
{  
    const Int nd = disc.common.nd;
    const Int npf = disc.common.npf;
    const Int nn = npf * nfaces;
    const Int n2 = 0;
    const Int n3 = nn;

    if (nd == 1) {
        std::vector<dstype> jacg(static_cast<size_t>(nn), 0.0);
        cpuFaceGeom1D(jacg.data(), nlint, xdgint, nn);
        return;
    }

    const Int derivStorage = (nd == 2) ? nn * nd : 2 * nn * nd;
    std::vector<dstype> scratch(static_cast<size_t>(n3 + derivStorage), 0.0);
    dstype* jacg = scratch.data() + n2;
    dstype* Jg = scratch.data() + n3;

    Node2Gauss(
        disc.common.cublasHandle, Jg, xdgint, &disc.master.shapfnt[npf * npf],
        npf, npf, nfaces * nd, 0);

    if (nd == 2) {
        cpuFaceGeom2D(jacg, nlint, Jg, nn);
        return;
    }

    Node2Gauss(
        disc.common.cublasHandle, &Jg[nn * nd], xdgint,
        &disc.master.shapfnt[2 * npf * npf], npf, npf, nfaces * nd,
        0);
    cpuFaceGeom3D(jacg, nlint, Jg, nn);
}

static inline void GatherWallGaussPointsAndNormals(
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
    wm.nfe = disc.common.nfe;
    wm.nbe1 = disc.common.nbe1;

    wm.nfaces = getFacesOnInterface(wm.faces, wm.nextfaces, disc, ibc);
    wm.npoints = wm.ngf * wm.nfaces;

    if (wm.nfaces <= 0) {
        std::ostringstream oss;
        oss << "GatherWallGaussPointsAndNormals found no faces for boundary condition " << ibc;
        error(oss.str());
    }

    wm.elems.resize(static_cast<size_t>(wm.npoints));
    for (Int iface = 0; iface < wm.nfaces; ++iface) {
        const Int elem = wm.faces[iface] / wm.nfe;
        for (Int ig = 0; ig < wm.ngf; ++ig)
            wm.elems[PointIndex(ig, iface, wm.ngf)] = elem;
    }

    dstype* xdgint = nullptr;
    dstype* nlint = nullptr;
    TemplateMalloc(&xdgint, wm.npf * wm.nfaces * wm.ncx, 0);
    TemplateMalloc(&nlint, wm.npf * wm.nfaces * wm.nd, 0);

    getNodesOnInterface(
        xdgint, disc.sol.xdg, wm.faces.data(), disc.mesh.perm, disc.common.nfe,
        disc.common.npf, disc.common.npe, disc.common.ncx, disc.common.ncx, wm.nfaces);
    getNormalVectorOnInterface(nlint, xdgint, disc, wm.nfaces);

    wm.xw.resize(static_cast<size_t>(wm.npoints * wm.nd));
    wm.nw.resize(static_cast<size_t>(wm.npoints * wm.nd));
    getFieldsAtGaussPointsOnInterface(wm.xw.data(), xdgint, disc, wm.nfaces, wm.nd);
    getFieldsAtGaussPointsOnInterface(wm.nw.data(), nlint, disc, wm.nfaces, wm.nd);

    wm.x1.resize(static_cast<size_t>(wm.npoints * wm.nd));
    ComputeOffWallPoints(wm.x1.data(), wm.xw.data(), wm.nw.data(), wm.y1, wm.nd, wm.npoints);

    TemplateFree(xdgint, 0);
    TemplateFree(nlint, 0);
}

static inline bool ComputeXi1AndShapeFunctionsFromElementFaces(
    WallModelSamplingData& wm,
    CDiscretization& disc,
    Int maxNewtonIter,
    dstype newtonTol,
    dstype insideTol)
{
    if (disc.sol.xdg == nullptr || disc.master.xpe == nullptr ||
        disc.mesh.e2f == nullptr || disc.mesh.f2e == nullptr)
        return false;
    if (wm.npoints <= 0 || wm.nd <= 0 || wm.npe <= 0 || wm.ncx <= 0)
        return false;
    if (static_cast<Int>(wm.x1.size()) != wm.nd * wm.npoints)
        return false;
    if (static_cast<Int>(wm.elems.size()) != wm.npoints)
        return false;

    wm.xi1.resize(static_cast<size_t>(wm.nd * wm.npoints));
    wm.shap1.resize(static_cast<size_t>(wm.npe * wm.npoints));

    const Int maxCandidates = disc.common.ne;
    if (maxCandidates <= 0)
        return false;

    std::vector<Int> pointElems(static_cast<size_t>(wm.npoints), -1);
    std::vector<Int> candidateElems(static_cast<size_t>(maxCandidates), -1);

    return FindPointsAndShapeFunctionsFromElementFaces(
        pointElems.data(),
        wm.xi1.data(),
        wm.shap1.data(),
        wm.x1.data(),
        wm.npoints,
        disc.sol.xdg,
        wm.elems.data(),
        candidateElems.data(),
        maxCandidates,
        disc.mesh.e2f,
        disc.mesh.f2e,
        wm.nfe,
        disc.master.xpe,
        wm.nd,
        wm.npe,
        wm.ncx,
        disc.common.elemtype,
        disc.common.porder,
        maxNewtonIter,
        newtonTol,
        insideTol);
}

static inline bool BuildWallModelSamplingData(
    WallModelSamplingData& wm,
    CDiscretization& disc,
    Int ibc,
    dstype y1,
    Int maxNewtonIter,
    dstype newtonTol,
    dstype insideTol)
{
    wm.y1 = y1;
    GatherWallGaussPointsAndNormals(wm, disc, ibc);
    return ComputeXi1AndShapeFunctionsFromElementFaces(
        wm, disc, maxNewtonIter, newtonTol, insideTol);
}

#endif
