#ifndef __QUADRATURE
#define __QUADRATURE

#include "gaussQuad.cpp"
#include "gaussLobattoQuad.cpp"

// Written by: C. Nguyen & P. Fernandez

void quadrature(vector<double> *x_p, vector<double> *w_p, Int *nq, Int pgauss, Int quadType, Int dim, Int elemtype)
{
    switch (quadType) {
        case 0:
            gaussQuad(x_p, w_p, nq, pgauss, dim, elemtype);
            break;
        case 1:
            gaussLobattoQuad(x_p, w_p, nq, pgauss, dim, elemtype);
            break;
        default:
            error("Quadrature rule not implemented.\n");
    }
}

#endif
