#include "wallmodelsampling.h"

namespace exasim {
namespace wm {

void ComputeOffWallPoints(
    std::vector<dstype>& x1,
    const std::vector<dstype>& xw,
    const std::vector<dstype>& nw,
    const dstype y1,
    const Int nd,
    const Int npoints)
{
    x1.resize(nd * npoints);
    ComputeOffWallPoints(x1.data(), xw.data(), nw.data(), y1, nd, npoints);
}

void ComputeOffWallPoints(
    dstype* x1,
    const dstype* xw,
    const dstype* nw,
    const dstype y1,
    const Int nd,
    const Int npoints)
{
    for (Int ip = 0; ip < npoints; ++ip) {
        for (Int d = 0; d < nd; ++d) {
            const Int k = VecIndex(d, ip, nd);
            x1[k] = xw[k] - y1 * nw[k];
        }
    }
}

} // namespace wm
} // namespace exasim
