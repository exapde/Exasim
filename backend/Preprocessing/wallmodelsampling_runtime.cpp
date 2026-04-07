#include "wallmodelsampling.h"

namespace exasim {
namespace wm {

// This translation unit is intentionally small. It provides a stable file
// for runtime-facing wrappers if the wall-model interpolation utilities need
// to be called from model or coupling code without pulling in the full
// preprocessing orchestration stack.

void EvaluateOffWallStateForBoundaryBlock(
    dstype* U1,
    const dstype* udg,
    const WallModelSamplingData& wm,
    const Int npe,
    const Int nc)
{
    EvaluateOffWallState(U1, udg, wm.e1.data(), wm.shap1.data(), npe, nc, wm.npoints);
}

} // namespace wm
} // namespace exasim
