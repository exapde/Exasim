function [f, tprd, t2t, f2t] = setboundaryperiodicfaces(p, t, elemtype, bndexpr, prdexpr)

[f, t2fl] = setboundary(p, t, elemtype, bndexpr);

[tprd, f] = setperiodic(p, t, f, t2fl, prdexpr);

dim = size(p,1);
[f2t, t2t] = mkf2e(tprd, elemtype, dim);

