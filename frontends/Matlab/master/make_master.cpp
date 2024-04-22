#ifndef __MKMASTER
#define __MKMASTER

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdlib.h>

using namespace std;

#include "mkmaster.h"
#include "errormsg.cpp"
#include "mkshape.cpp"
#include "quadrature.cpp"
#include "masternodes.cpp"

// Written by: C. Nguyen & P. Fernandez

void mkmaster(masterstruct &master, Int porder, Int *pgauss, Int *quadType, Int nd, Int elemtype, Int nodetype)
{
    Int i, j, k, l;
    
    Int pgaussR = pgauss[0];
    Int pgaussJ = pgauss[1];
    Int pgaussQ = pgauss[2];
    master.pgaussR = pgaussR;
    master.pgaussJ = pgaussJ;
    master.pgaussQ = pgaussQ;
    master.quadTypeR = quadType[0];
    master.quadTypeJ = quadType[1];
    master.quadTypeQ = quadType[2];
    master.nd = nd;                     // Dimension of master element
    master.porder = porder;             // Polynomial degree
    master.elemtype = elemtype;
    master.nodetype = nodetype;
    
    // High-order nodes on the master element and face:
    masternodes(&master.plocvl, &master.plocfc, &master.npv, &master.npf, master.porder, master.nd, master.elemtype, master.nodetype);
    Int npv = master.npv;
    Int npf = master.npf;
    
    // Quadrature rule on the master element:
    quadrature(&master.gpvlR, &master.gwvlR, &master.nqvR, master.pgaussR, master.quadTypeR, nd, elemtype);
    quadrature(&master.gpvlJ, &master.gwvlJ, &master.nqvJ, master.pgaussJ, master.quadTypeJ, nd, elemtype);
    quadrature(&master.gpvlQ, &master.gwvlQ, &master.nqvQ, master.pgaussQ, master.quadTypeQ, nd, elemtype);
    Int nqvR = master.nqvR;
    Int nqvJ = master.nqvJ;
    Int nqvQ = master.nqvQ;
    
    // Quadrature rule on the master face:
    quadrature(&master.gpfcR, &master.gwfcR, &master.nqfR, pgaussR, master.quadTypeR, nd-1, elemtype);
    quadrature(&master.gpfcJ, &master.gwfcJ, &master.nqfJ, pgaussJ, master.quadTypeJ, nd-1, elemtype);
    quadrature(&master.gpfcQ, &master.gwfcQ, &master.nqfQ, pgaussQ, master.quadTypeQ, nd-1, elemtype);
    Int nqfR = master.nqfR;
    Int nqfJ = master.nqfJ;
    Int nqfQ = master.nqfQ;
    
    // Shape functions and derivatives on the master element:
    mkshape(&master.shapvlR, &master.plocvl[0], &master.gpvlR[0], nqvR, elemtype, porder, nd, npv);
    master.shapvtR.resize(npv*nqvR*(nd+1));
    master.shapvgR.resize(npv*(nd+1)*nqvR);
    master.shapvgdotshapvlR.resize(npv*npv*(nd+1)*nqvR);
    for (i = 0; i < (nd+1); i++)
        for (j = 0; j < nqvR; j++)
            for (k = 0; k < npv; k++) {
                master.shapvtR[i*npv*nqvR+k*nqvR+j] = master.shapvlR[i*nqvR*npv+j*npv+k];
                master.shapvgR[i*nqvR*npv+j*npv+k] = master.gwvlR[j] * master.shapvlR[i*nqvR*npv+j*npv+k];
            }
    for (i = 0; i < (nd+1); i++)
        for (j = 0; j < nqvR; j++)
            for (k = 0; k < npv; k++)
                for (l = 0; l < npv; l++)
                    master.shapvgdotshapvlR[i*nqvR*npv*npv + j*npv*npv + k*npv + l] = master.shapvgR[i*nqvR*npv + j*npv + l] * master.shapvlR[0*nqvR*npv + j*npv + k];
    
    mkshape(&master.shapnv, &master.plocvl[0], &master.plocvl[0], npv, elemtype, porder, nd, npv);
    master.shapnvt.resize(npv*npv*(nd+1));
    for (i = 0; i < (nd+1); i++)
        for (j = 0; j < npv; j++)
            for (k = 0; k < npv; k++) {
                master.shapnvt[i*npv*npv+k*npv+j] = master.shapnv[i*npv*npv+j*npv+k];
            }
    
    mkshape(&master.shapvlJ, &master.plocvl[0], &master.gpvlJ[0], nqvJ, elemtype, porder, nd, npv);
    master.shapvtJ.resize(npv*nqvJ*(nd+1));
    master.shapvgJ.resize(npv*(nd+1)*nqvJ);
    master.shapvgdotshapvlJ.resize(npv*npv*(nd+1)*nqvJ);
    for (i = 0; i < (nd+1); i++)
        for (j = 0; j < nqvJ; j++)
            for (k = 0; k < npv; k++) {
                master.shapvtJ[i*npv*nqvJ+k*nqvJ+j] = master.shapvlJ[i*nqvJ*npv+j*npv+k];
                master.shapvgJ[i*nqvJ*npv+j*npv+k] = master.gwvlJ[j] * master.shapvlJ[i*nqvJ*npv+j*npv+k];
            }
    for (i = 0; i < (nd+1); i++)
        for (j = 0; j < nqvJ; j++)
            for (k = 0; k < npv; k++)
                for (l = 0; l < npv; l++)
                    master.shapvgdotshapvlJ[i*nqvJ*npv*npv + j*npv*npv + k*npv + l] = master.shapvgJ[i*nqvJ*npv + j*npv + l] * master.shapvlJ[0*nqvJ*npv + j*npv + k];
    
    mkshape(&master.shapvlQ, &master.plocvl[0], &master.gpvlQ[0], nqvQ, elemtype, porder, nd, npv);
    master.shapvtQ.resize(npv*nqvQ*(nd+1));
    master.shapvgQ.resize(npv*(nd+1)*nqvQ);
    master.shapvgdotshapvlQ.resize(npv*npv*(nd+1)*nqvQ);
    for (i = 0; i < (nd+1); i++)
        for (j = 0; j < nqvQ; j++)
            for (k = 0; k < npv; k++) {
                master.shapvtQ[i*npv*nqvQ+k*nqvQ+j] = master.shapvlQ[i*nqvQ*npv+j*npv+k];
                master.shapvgQ[i*nqvQ*npv+j*npv+k] = master.gwvlQ[j] * master.shapvlQ[i*nqvQ*npv+j*npv+k];
            }
    for (i = 0; i < (nd+1); i++)
        for (j = 0; j < nqvQ; j++)
            for (k = 0; k < npv; k++)
                for (l = 0; l < npv; l++)
                    master.shapvgdotshapvlQ[i*nqvQ*npv*npv + j*npv*npv + k*npv + l] = master.shapvgQ[i*nqvQ*npv + j*npv + l] * master.shapvlQ[0*nqvQ*npv + j*npv + k];
    
    // Shape functions and derivatives on the master face:
    if (nd > 1) {
        mkshape(&master.shapfcR, &master.plocfc[0], &master.gpfcR[0], nqfR, elemtype, porder, nd-1, npf);      
        master.shapftR.resize(npf*nqfR*nd);
        master.shapfgR.resize(npf*nd*nqfR);
        master.shapfgdotshapfcR.resize(npf*npf*nd*nqfR);
        for (i = 0; i < nd; i++)
            for (j = 0; j < nqfR; j++)
                for (k = 0; k < npf; k++) {
                    master.shapftR[i*npf*nqfR+k*nqfR+j] = master.shapfcR[i*nqfR*npf+j*npf+k];
                    master.shapfgR[i*nqfR*npf+j*npf+k] = master.gwfcR[j] * master.shapfcR[i*nqfR*npf+j*npf+k];
                }
        for (i = 0; i < nd; i++)
            for (j = 0; j < nqfR; j++)
                for (k = 0; k < npf; k++)
                    for (l = 0; l < npf; l++)
                        master.shapfgdotshapfcR[i*nqfR*npf*npf + j*npf*npf + k*npf + l] = master.shapfgR[i*nqfR*npf + j*npf + l] * master.shapfcR[0*nqfR*npf + j*npf + k];

        mkshape(&master.shapfcJ, &master.plocfc[0], &master.gpfcJ[0], nqfJ, elemtype, porder, nd-1, npf);      
        master.shapftJ.resize(npf*nqfJ*nd);
        master.shapfgJ.resize(npf*nd*nqfJ);
        master.shapfgdotshapfcJ.resize(npf*npf*nd*nqfJ);
        for (i = 0; i < nd; i++)
            for (j = 0; j < nqfJ; j++)
                for (k = 0; k < npf; k++) {
                    master.shapftJ[i*npf*nqfJ+k*nqfJ+j] = master.shapfcJ[i*nqfJ*npf+j*npf+k];
                    master.shapfgJ[i*nqfJ*npf+j*npf+k] = master.gwfcJ[j] * master.shapfcJ[i*nqfJ*npf+j*npf+k];
                }
        for (i = 0; i < nd; i++)
            for (j = 0; j < nqfJ; j++)
                for (k = 0; k < npf; k++)
                    for (l = 0; l < npf; l++)
                        master.shapfgdotshapfcJ[i*nqfJ*npf*npf + j*npf*npf + k*npf + l] = master.shapfgJ[i*nqfJ*npf + j*npf + l] * master.shapfcJ[0*nqfJ*npf + j*npf + k];

        mkshape(&master.shapfcQ, &master.plocfc[0], &master.gpfcQ[0], nqfQ, elemtype, porder, nd-1, npf);      
        master.shapftQ.resize(npf*nqfQ*nd);
        master.shapfgQ.resize(npf*nd*nqfQ);
        master.shapfgdotshapfcQ.resize(npf*npf*nd*nqfQ);
        for (i = 0; i < nd; i++)
            for (j = 0; j < nqfQ; j++)
                for (k = 0; k < npf; k++) {
                    master.shapftQ[i*npf*nqfQ+k*nqfQ+j] = master.shapfcQ[i*nqfQ*npf+j*npf+k];
                    master.shapfgQ[i*nqfQ*npf+j*npf+k] = master.gwfcQ[j] * master.shapfcQ[i*nqfQ*npf+j*npf+k];
                }
        for (i = 0; i < nd; i++)
            for (j = 0; j < nqfQ; j++)
                for (k = 0; k < npf; k++)
                    for (l = 0; l < npf; l++)
                        master.shapfgdotshapfcQ[i*nqfQ*npf*npf + j*npf*npf + k*npf + l] = master.shapfgQ[i*nqfQ*npf + j*npf + l] * master.shapfcQ[0*nqfQ*npf + j*npf + k];
    }
    else {
        master.shapfcR.resize(1);
        master.shapfcR[0] = 1;
        master.shapftR.resize(1);
        master.shapftR[0] = 1;
        master.shapfgR.resize(1);
        master.shapfgR[0] = 1;
        master.shapfgdotshapfcR.resize(1);
        master.shapfgdotshapfcR[0] = 1;
        master.shapfcJ.resize(1);
        master.shapfcJ[0] = 1;
        master.shapftJ.resize(1);
        master.shapftJ[0] = 1;
        master.shapfgJ.resize(1);
        master.shapfgJ[0] = 1;
        master.shapfgdotshapfcJ.resize(1);
        master.shapfgdotshapfcJ[0] = 1;
        master.shapfcQ.resize(1);
        master.shapfcQ[0] = 1;
        master.shapftQ.resize(1);
        master.shapftQ[0] = 1;
        master.shapfgQ.resize(1);
        master.shapfgQ[0] = 1;
        master.shapfgdotshapfcQ.resize(1);
        master.shapfgdotshapfcQ[0] = 1;
    }
}

#endif
