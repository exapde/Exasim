#include "wallmodelsampling.h"

#include <cmath>

namespace exasim {
namespace wm {

void EvaluateOffWallStateAtPoint(
    dstype* U1p,
    const dstype* udg,
    const Int elem,
    const dstype* N,
    const Int npe,
    const Int nc)
{
    const Int eoff = npe * nc * elem;
    for (Int c = 0; c < nc; ++c) {
        dstype val = 0.0;
        for (Int a = 0; a < npe; ++a)
            val += N[a] * udg[eoff + a + npe * c];
        U1p[c] = val;
    }
}

void EvaluateOffWallState(
    dstype* U1,
    const dstype* udg,
    const Int* e1,
    const dstype* shap1,
    const Int npe,
    const Int nc,
    const Int npoints)
{
    for (Int ip = 0; ip < npoints; ++ip) {
        EvaluateOffWallStateAtPoint(
            &U1[nc * ip], udg, e1[ip], &shap1[npe * ip], npe, nc);
    }
}

void ConservativeToPrimitiveAtOffWallPoint(
    dstype& rho1,
    dstype* vel1,
    dstype& p1,
    dstype& T1,
    const dstype* U1,
    const dstype* param,
    const Int nd,
    const Int nc)
{
    (void) nc;

    const dstype gam = param[0];
    const dstype Minf = param[3];
    const dstype gam1 = gam - 1.0;

    rho1 = U1[0];
    for (Int d = 0; d < nd; ++d) vel1[d] = U1[1 + d] / rho1;

    dstype ke = 0.0;
    for (Int d = 0; d < nd; ++d) ke += vel1[d] * vel1[d];
    ke *= 0.5;

    const dstype rE = U1[nd + 1];
    p1 = gam1 * (rE - rho1 * ke);
    T1 = gam * Minf * Minf * p1 / rho1;
}

dstype ComputeTangentialSpeed(
    dstype* that,
    const dstype* vel1,
    const dstype* n,
    const Int nd)
{
    dstype un = 0.0;
    for (Int d = 0; d < nd; ++d) un += vel1[d] * n[d];

    dstype utmag2 = 0.0;
    for (Int d = 0; d < nd; ++d) {
        that[d] = vel1[d] - un * n[d];
        utmag2 += that[d] * that[d];
    }

    const dstype utmag = std::sqrt(utmag2);
    if (utmag > 0.0) {
        for (Int d = 0; d < nd; ++d) that[d] /= utmag;
    } else {
        for (Int d = 0; d < nd; ++d) that[d] = 0.0;
    }

    return utmag;
}

dstype ComputeBFWMInput(
    const dstype y1,
    const dstype utmag1,
    const dstype rho1,
    const dstype mu1,
    const dstype k1,
    const dstype cp)
{
    return y1 * utmag1 * rho1 * mu1 * cp * cp / (k1 * k1);
}

} // namespace wm
} // namespace exasim
