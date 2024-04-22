#ifndef __MASTERNODES
#define __MASTERNODES

// Written by: C. Nguyen & P. Fernandez

void uniformNodes(double *xi, Int porder)
{
    // Uniform nodes on the interval [0, 1]:
    
    Int numNodes1d = porder + 1;
    
    if (porder == 0)
        xi[0] = 0.5;
    else {
        for (int k = 0; k < numNodes1d; k++)
            xi[k] = ((double) k) / ((double) numNodes1d - 1.0);
    }
}

void GaussLobattoNodes(double *xi, Int porder)
{
    // Gauss-Lobatto nodes on the interval [0, 1]. Ref: http://mathworld.wolfram.com/LobattoQuadrature.html
    
    printf("Gauss-Lobatto nodes not validated yet.\n");
    
    switch (porder) {
        case 0:
            xi[0] = 0.5;
            break;
        case 1:
            xi[0] = 0.0;
            xi[1] = 1.0;
            break;
        case 2:
            xi[0] = 0.0;
            xi[1] = 0.5;
            xi[2] = 1.0;
            break;
        case 3:
            xi[0] =  0.0;
            xi[1] = -0.5 * sqrt(5.0) / 5.0 + 0.5;
            xi[2] =  0.5 * sqrt(5.0) / 5.0 + 0.5;
            xi[3] =  1.0;
            break;
        case 4:
            xi[0] =  0.0;
            xi[1] = -0.5*sqrt(21.0) / 7.0 + 0.5;
            xi[2] =  0.5;
            xi[3] =  0.5*sqrt(21.0) / 7.0 + 0.5;
            xi[4] =  1.0;
            break;
        case 5:
            xi[0] =  0.0;
            xi[1] = -0.5 * sqrt((7.0+2.0*sqrt(7.0)) / 21.0) + 0.5;
            xi[2] = -0.5 * sqrt((7.0-2.0*sqrt(7.0)) / 21.0) + 0.5;
            xi[3] =  0.5 * sqrt((7.0-2.0*sqrt(7.0)) / 21.0) + 0.5;
            xi[4] =  0.5 * sqrt((7.0+2.0*sqrt(7.0)) / 21.0) + 0.5;
            xi[5] =  1.0;
            break;
        default:
            error("Gauss-Lobatto nodes not implemented for porder > 5.\n");
    }
}

void ChebyshevNodes(double *xi, Int porder)
{
    // Chebyshev nodes on the interval [0, 1]:
    
    printf("Chebyshev nodes not validated yet.\n");
    
    Int numNodes = porder + 1;
    double PI = 3.1415926535897932;
    
    if (porder == 0)
        xi[0] = 0.5;
    else {
        for (int k = 0; k < numNodes; k++) {
            xi[k] = - cos(PI * (2.0 * (double) k + 1.0) / (2.0 * (double) numNodes));
            xi[k] = 0.5 + 0.5*xi[k];
        }
    }
}

void extendedChebyshevNodes(double *xi, Int porder)
{
    // Extended Chebyshev nodes on the interval [0, 1]:
    
    printf("Extended Chebyshev nodes not validated yet.\n");
    
    Int numNodes = porder + 1;
    double PI = 3.1415926535897932;
    
    if (porder == 0)
        xi[0] = 0.5;
    else {
        for (int k = 0; k < numNodes; k++) {
            xi[k] = - cos(PI * (2.0 * (double) k + 1.0) / (2.0 * (double) numNodes)) / cos(PI / (2.0 * (double) numNodes));
            xi[k] = 0.5 + 0.5*xi[k];
        }
    }
}

void nodes1d(double *xi, Int porder, Int nodetype)
{
    switch (nodetype) {
        case 0:     // Uniform
            uniformNodes(xi, porder);
            break;
        case 1:     // Gauss-Lobatto
            GaussLobattoNodes(xi, porder);
            break;
        case 2:     // Chebyshev
            ChebyshevNodes(xi, porder);
            break;
        case 3:     // Extended Chebyshev
            extendedChebyshevNodes(xi, porder);
            break;
        default:
            error("nodetype not implemented.\n");
    }
}

void lineNodes1d(vector<double> *plocvl, Int *npv, Int porder, Int nodetype)
{
    plocvl[0].resize(porder+1);
    double *xi = &plocvl[0][0];
    
    *npv = porder + 1;
    nodes1d(xi, porder, nodetype);
}

void nodesTensor(vector<double> *plocvl, vector<double> *plocfc, Int *npv, Int *npf, Int porder, Int dim, Int nodetype)
{
    // x: numNodes / nd
    
    Int i, j, k;
    
    Int numNodes1d = porder+1;
    Int numNodes2d = numNodes1d * numNodes1d;
    Int numNodes3d = numNodes1d * numNodes1d * numNodes1d;
    
    vector<double> ploc1d;
    ploc1d.resize(numNodes1d*1); 
    double *xi = &ploc1d[0];
    nodes1d(xi, porder, nodetype);
    
    switch (dim) {
        case 1:
            *npv = numNodes1d;
            *npf = 1;
            plocvl[0].resize(numNodes1d*dim);
            for (i = 0; i < numNodes1d; i++) {
                plocvl[0][i] = ploc1d[i];
            }
            break;
        case 2:
            *npv = numNodes2d;
            *npf = numNodes1d;
            plocvl[0].resize(numNodes2d*dim);
            plocfc[0].resize(numNodes1d*(dim-1));
            for (i = 0; i < numNodes1d; i++) {
                for (j = 0; j < numNodes1d; j++) {
                    plocvl[0][0*numNodes2d+i*numNodes1d+j] = ploc1d[j];
                    plocvl[0][1*numNodes2d+i*numNodes1d+j] = ploc1d[i];
                }
            }   
            for (i = 0; i < numNodes1d; i++) {
                plocfc[0][i] = ploc1d[i];
            }
            break;
        case 3:
            *npv = numNodes3d;
            *npf = numNodes2d;
            plocvl[0].resize(numNodes3d*dim);
            plocfc[0].resize(numNodes2d*(dim-1));
            for (i = 0; i < numNodes1d; i++) {
                for (j = 0; j < numNodes1d; j++) {
                    for (k = 0; k < numNodes1d; k++) {
                        plocvl[0][0*numNodes3d+i*numNodes2d+j*numNodes1d+k] = ploc1d[k];
                        plocvl[0][1*numNodes3d+i*numNodes2d+j*numNodes1d+k] = ploc1d[j];
                        plocvl[0][2*numNodes3d+i*numNodes2d+j*numNodes1d+k] = ploc1d[i];
                    }
                }
            }
            for (i = 0; i < numNodes1d; i++) {
                for (j = 0; j < numNodes1d; j++) {
                    plocfc[0][0*numNodes2d+i*numNodes1d+j] = ploc1d[j];
                    plocfc[0][1*numNodes2d+i*numNodes1d+j] = ploc1d[i];
                }
            }
            break;
        default:
            error("Dimension not implemented.\n");
    }
}

void quadnodes2d(vector<double> *plocvl, vector<double> *plocfc, Int *npv, Int *npf, Int porder, Int nodetype)
{
    Int dim = 2;
    nodesTensor(plocvl, plocfc, npv, npf, porder, dim, nodetype);
}

void hexnodes3d(vector<double> *plocvl, vector<double> *plocfc, Int *npv, Int *npf, Int porder, Int nodetype)
{
    Int dim = 3;
    nodesTensor(plocvl, plocfc, npv, npf, porder, dim, nodetype);
}

void trinodes2d(vector<double> *plocvl, vector<double> *plocfc, Int *npv, Int *npf, Int porder, Int nodetype)
{
    Int in, dim = 2;
    Int numNodes1d = porder + 1;
    Int numNodesTri =  ((porder + 1) * (porder + 2)) / 2;
    
    *npv = numNodesTri;
    *npf = numNodes1d;
    
    plocvl[0].resize(numNodesTri*dim);
    plocfc[0].resize(numNodes1d*(dim-1));
    
    if (porder == 0) {
        plocvl[0][0] = 1.0 / 3.0;
        plocvl[0][1] = 1.0 / 3.0;
        plocfc[0][0] = 1.0 / 2.0;
    }
    else {
        switch (nodetype) {
            case 0:     // Uniform
                in = 0;
                for (int j = 0; j < numNodes1d; j++)
                    for (int i = 0; i < numNodes1d-j; i++) {
                        plocvl[0][0*numNodesTri+in] = ((double) i) / ((double) numNodes1d - 1.0);
                        plocvl[0][1*numNodesTri+in] = ((double) j) / ((double) numNodes1d - 1.0);
                        in++;
                    }
                for (int i = 0; i < numNodes1d; i++)
                    plocfc[0][i] = ((double) i) / ((double) numNodes1d - 1.0);
                break;
            default:
                error("nodetype not implemented for tris.\n");
        }
    }
}

void tetnodes3d(vector<double> *plocvl, vector<double> *plocfc, Int *npv, Int *npf, Int porder, Int nodetype)
{
    Int in, dim = 3;
    Int numNodes1d = porder + 1;
    Int numNodesTri =  ((porder + 1) * (porder + 2)) / 2;
    Int numNodesTet =  ((porder + 1) * (porder + 2) * (porder + 3)) / 6;
    
    *npv = numNodesTet;
    *npf = numNodesTri;
    
    plocvl[0].resize(numNodesTet*dim);
    plocfc[0].resize(numNodesTri*(dim-1));
    
    if (porder == 0) {
        plocvl[0][0] = 1.0 / 4.0;
        plocvl[0][1] = 1.0 / 4.0;
        plocvl[0][2] = 1.0 / 4.0;
        plocfc[0][0] = 1.0 / 3.0;
        plocfc[0][1] = 1.0 / 3.0;
    }
    else {
        switch (nodetype) {
            case 0:     // Uniform
                in = 0;
                for (int k = 0; k < numNodes1d; k++)
                    for (int j = 0; j < numNodes1d-k; j++)
                        for (int i = 0; i < numNodes1d-(j+k); i++) {
                            plocvl[0][0*numNodesTet+in] = ((double) i) / ((double) numNodes1d - 1.0);
                            plocvl[0][1*numNodesTet+in] = ((double) j) / ((double) numNodes1d - 1.0);
                            plocvl[0][2*numNodesTet+in] = ((double) k) / ((double) numNodes1d - 1.0);
                            in++;
                        }
                in = 0;
                for (int j = 0; j < numNodes1d; j++)
                    for (int i = 0; i < numNodes1d-j; i++) {
                        plocfc[0][0*numNodesTri+in] = ((double) i) / ((double) numNodes1d - 1.0);
                        plocfc[0][1*numNodesTri+in] = ((double) j) / ((double) numNodes1d - 1.0);
                        in++;
                    }
                break;
            default:
                error("nodetype not implemented for tets.\n");
        }
    }
}

void masternodes(vector<double> *plocvl, vector<double> *plocfc, Int *npv, Int *npf, Int porder, Int dim, Int elemtype, Int nodetype)
{
    switch (dim) {
        case 1: // 1D
            *npf = 1;
            plocfc[0].resize(1);
            plocfc[0][0] = 0.0;
            lineNodes1d(plocvl, npv, porder, nodetype);
            break;
        case 2: // 2D
            if (elemtype==0)     // tri 
                trinodes2d(plocvl, plocfc, npv, npf, porder, nodetype);
            else if (elemtype==1) // quad
                quadnodes2d(plocvl, plocfc, npv, npf, porder, nodetype);
            break;
        case 3: // 3D
            if (elemtype==0)     // tet
                tetnodes3d(plocvl, plocfc, npv, npf, porder, nodetype);
            else if (elemtype==1) // hex
                hexnodes3d(plocvl, plocfc, npv, npf, porder, nodetype);
            break;
        default:
            error("Dimension not implemented.\n");
    }
}

#endif
