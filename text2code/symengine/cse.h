#pragma once

#include <symengine/basic.h>
#include <symengine/symengine_rcp.h>
#include <symengine/visitor.h>

namespace SymEngine {

// using vec_pair = std::vector<std::pair<RCP<const Basic>, RCP<const Basic>>>;
// using vec_basic = std::vector<RCP<const Basic>>;

// Common Subexpression Elimination (CSE)
void cse(vec_pair &replacements, vec_basic &reduced_exprs, const vec_basic &exprs);

}
