#include "mex.h"
#include <math.h>
//#include <stdint.h>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

// Written by: C. Nguyen & P. Fernandez

using namespace std;
typedef int Int;

Int IsSameFace(vector<Int> face1, vector<Int> face2)
{
    Int same = 0;
    if (face1.size() == face2.size()) {
        same = 1;
        Int i, n = face1.size();
        for (i=0; i<n; i++)
           if (face1[i] != face2[i]) {
               same = 0;
               break;
           }
    }         
    return same;
}

struct elementinfo {    
    Int dim;      // spatial dimension
    Int elemtype; // type     
    
    Int nfe; // number of faces    
    Int nle;  // number of edges 
    Int nve; // number of vertices 
    Int nqf; // number of quad faces
    Int ntf; // number of tet faces
    vector< Int > facetype; // 0 triangular, 1 quadrilateral 
    vector< vector< Int > > face; // face connectivity     
    vector< vector< Int > > edge; // edge connectivity     
};

elementinfo mkeleminfo(Int dim, Int elemtype)
{    
    Int nfe, nle, nve, nqf, ntf;
    vector< Int > facetype;
    vector< vector< Int > > face;
    
    switch (dim) {
        case 1: // line element
            nfe = 2;
            nle = 0;
            nve = 2;
            nqf = 0;
            ntf = 0;
            facetype.resize(nfe,0);
            face.resize(nfe);
            face[0].resize(1);
            face[1].resize(1);
            face[0][0] = 0;
            face[1][0] = 1;
            break;
        case 2:            
            if (elemtype==0) { // triangular
                nfe = dim+1;  
                face.resize(nfe);
                
                // [[1,2];[2 0]; [0,1]]
                face[0].resize(2);
                face[0][0] = 1;
                face[0][1] = 2;
                
                face[1].resize(2);
                face[1][0] = 2;
                face[1][1] = 0;
                
                face[2].resize(2);
                face[2][0] = 0;
                face[2][1] = 1;     
            }
            else if (elemtype==1) { // quadrilateral
                nfe = 2*dim;   
                face.resize(nfe);
                
                // [[1,2];[2,3];[3,4];[4,1]] - 1;
                face[0].resize(2);
                face[0][0] = 0;
                face[0][1] = 1;
                
                face[1].resize(2);
                face[1][0] = 1;
                face[1][1] = 2;
                
                face[2].resize(2);
                face[2][0] = 2;
                face[2][1] = 3;                            
                
                face[3].resize(2);
                face[3][0] = 3;
                face[3][1] = 0;    
            }
            nle = nfe;
            nve = nfe;
            nqf = 0;
            ntf = 0;
            facetype.resize(nfe,0);
            break;
        case 3:
            if (elemtype==0) { // tet
                nfe = dim+1;
                nve = 4;
                nle = 6;
                ntf = 4;
                nqf = 0;
                
                facetype.resize(nfe,0);                
                face.resize(nfe);
                        
                // [[2,3,4];[1,4,3];[1,2,4];[1,3,2]] - 1
                face[0].resize(3);
                face[0][0] = 1;
                face[0][1] = 2;
                face[0][2] = 3;
                
                face[1].resize(3);
                face[1][0] = 0;
                face[1][1] = 3;
                face[1][2] = 2;
                
                face[2].resize(3);
                face[2][0] = 0;
                face[2][1] = 1;
                face[2][2] = 3;
                
                face[3].resize(3);
                face[3][0] = 0;
                face[3][1] = 2;
                face[3][2] = 1;                
            }
            else if (elemtype==1) { // hex
                nfe = 2*dim;       
                nve = 8;
                nle = 12;
                nqf = nfe;
                ntf = 0;
                
                facetype.resize(nfe,1);                
                face.resize(nfe);
                
                //face=[[1,4,3,2];[5,6,7,8];[1,2,6,5];[3,4,8,7];[2,3,7,6];[4,1,5,8]] - 1;                
                face[0].resize(4);
                face[0][0] = 0;
                face[0][1] = 3;
                face[0][2] = 2;
                face[0][3] = 1;
                
                face[1].resize(4);
                face[1][0] = 4;
                face[1][1] = 5;
                face[1][2] = 6;
                face[1][3] = 7;
                
                face[2].resize(4);
                face[2][0] = 0;
                face[2][1] = 1;
                face[2][2] = 5;
                face[2][3] = 4;
                
                face[3].resize(4);
                face[3][0] = 2;
                face[3][1] = 3;
                face[3][2] = 7;                
                face[3][3] = 6;                
                
                face[4].resize(4);
                face[4][0] = 1;
                face[4][1] = 2;
                face[4][2] = 6;                
                face[4][3] = 5;         
                
                face[5].resize(4);
                face[5][0] = 3;
                face[5][1] = 0;
                face[5][2] = 4;                
                face[5][3] = 7;         
                //face=[[1,4,3,2];[5,6,7,8];[1,2,6,5];[3,4,8,7];[2,3,7,6];[4,1,5,8]] - 1;                
            }
            else if (elemtype==2) { // prism
                nfe = 5;                          
                nve = 6;
                nle = 9;
                nqf = 3;
                ntf = 2;
                
                facetype.resize(nfe,0);                
                facetype[0] = 0;
                facetype[1] = 0;
                facetype[2] = 1;
                facetype[3] = 1;
                facetype[4] = 1;
                
                face.resize(nfe);                
                //face=[[0,2,1];[3,4,5];[1,2,5,4];[2,0,3,5];[0,1,4,3]];                
                face[0].resize(3);
                face[0][0] = 0;
                face[0][1] = 2;
                face[0][2] = 1;                
                
                face[1].resize(3);
                face[1][0] = 3;
                face[1][1] = 4;
                face[1][2] = 5;                
                
                face[3].resize(4);
                face[3][0] = 1;
                face[3][1] = 2;
                face[3][2] = 5;                
                face[3][3] = 4;                
                
                face[2].resize(4);
                face[2][0] = 2;
                face[2][1] = 0;
                face[2][2] = 3;
                face[2][3] = 5;                                
                
                face[4].resize(4);
                face[4][0] = 0;
                face[4][1] = 1;
                face[4][2] = 4;                
                face[4][3] = 3;                         
                //face=[[0,2,1];[3,4,5];[1,2,5,4];[2,0,3,5];[0,1,4,3]];                                
            }
            else if (elemtype==3) { // pyramid
                nfe = 5;
                nve = 5;
                nle = 8;
                nqf = 1;
                ntf = 4;
                
                facetype.resize(nfe,0);                
                facetype[0] = 1;
                facetype[1] = 0;
                facetype[2] = 0;
                facetype[3] = 0;
                facetype[4] = 0;
                
                face.resize(nfe);
                //face=[[0,3,2,1];[0,1,4];[1,2,4];[2,3,4];[3,0,4]];                
                face[0].resize(4);
                face[0][0] = 0;
                face[0][1] = 3;
                face[0][2] = 2;
                face[0][3] = 1;
                
                face[1].resize(3);
                face[1][0] = 0;
                face[1][1] = 1;
                face[1][2] = 4;
                
                face[2].resize(3);
                face[2][0] = 1;
                face[2][1] = 2;
                face[2][2] = 4;
                
                face[3].resize(3);
                face[3][0] = 2;
                face[3][1] = 3;
                face[3][2] = 4;                
                
                face[4].resize(3);
                face[4][0] = 3;
                face[4][1] = 0;
                face[4][2] = 4;                
            }
            break;
        default:
            mexErrMsgTxt("Only can handle dim=1, dim=2 or dim=3\n");
    }
    
    elementinfo eleminfo;
    eleminfo.nfe = nfe;
    eleminfo.nle = nle;
    eleminfo.nve = nve;
    eleminfo.nqf = nqf;
    eleminfo.ntf = ntf;
    eleminfo.facetype = facetype;
    eleminfo.face = face;
    return eleminfo;
}

vector< elementinfo > mkeleminfos(Int dim)
{
    Int nelemtype;
    switch (dim) {
        case 1: // line element            
            nelemtype = 1;
            break;            
        case 2:
            nelemtype = 2;
            break;
        case 3:
            nelemtype = 4;
            break;
        default:
            mexErrMsgTxt("Only can handle dim=1, dim=2 or dim=3\n");
    }
    
    vector< elementinfo > eleminfos(nelemtype, elementinfo());
    
    for (Int i=0; i<nelemtype; i++)
        eleminfos[i] = mkeleminfo(dim, i);
    
    return eleminfos;    
}

vector< vector<Int> > mkf(vector< vector<Int> > &face, vector<Int> &t)
{
    
    Int i, j, nfe;
    nfe = face.size();
    
    vector< vector<Int> > f = face;    
    for (i=0; i<nfe; i++) 
        for (j = 0; j<f[i].size(); j++)
            f[i][j] = t[face[i][j]];    
    
    return f;
}
 
void mexFunction(int nlhs,mxArray* plhs[], int nrhs, const mxArray* prhs[])
{    
    Int i, j, k, m, n, nve, nfe, nvf, ne, same;
    Int dim, nvemax, nfemax, nvfmax, ntemax, elemtype, nelemtype;
                    
    double* t = mxGetPr(prhs[0]);
    double* elementtype = mxGetPr(prhs[1]);
    double* nd = mxGetPr(prhs[2]);
    dim = (Int) nd[0];     
    ne = (Int) nd[1];      /* number of elements */
    nvemax  = (Int) nd[2]; /* number of vertices per element */        
    nfemax = (Int) nd[3];  /* number of faces per element */        
    nvfmax = (Int) nd[4];  /* number of vertices per face */        
    
    vector< elementinfo> eleminfos = mkeleminfos(dim);    
        
    nelemtype = eleminfos.size();       
    vector< Int> nfes(nelemtype,0), nves(nelemtype,0);    
    for (i=0; i<nelemtype; i++) {
        nfes[i] = eleminfos[i].nfe;
        nves[i] = eleminfos[i].nve;
    }        
    
    vector< Int > ti(nvemax, 0);        
    vector< vector< vector<Int> > > faces(ne, vector< vector<Int> >());               
    for (i=0; i<ne; i++)
    {        
        elemtype = (Int) elementtype[i];       
        nve = nves[elemtype];        
        
        for (j=0; j<nve; j++)
            ti[j] = (Int) t[i*nvemax+j];        
        
        faces[i] = mkf(eleminfos[elemtype].face, ti);                
        
        nfe = nfes[elemtype];
        for(j = 0; j<nfe; j++)
            sort(faces[i][j].begin(), faces[i][j].end());
    }
                
    Int nmax = (1+nvfmax+2);     
    mwSize szQ[2]; 
    szQ[0]=nfemax; szQ[1] = ne; 
    plhs[0]=mxCreateNumericArray(2,szQ,mxDOUBLE_CLASS,mxREAL);  
    plhs[1]=mxCreateNumericArray(2,szQ,mxDOUBLE_CLASS,mxREAL);  
    szQ[0]=nmax; szQ[1] = nfemax*ne; 
    plhs[2]=mxCreateNumericArray(2,szQ,mxDOUBLE_CLASS,mxREAL);  
        
    double *t2t, *t2f, *f;
    t2t = mxGetPr(plhs[0]);
    t2f = mxGetPr(plhs[1]);
    f = mxGetPr(plhs[2]);
    for (i=0; i<ne*nfemax; i++) {
        t2t[i] = -1;
        t2f[i] = -1;        
    }
    for (i=0; i<ne*nfemax*nmax; i++)
        f[i] = -1;
    
    Int nf = 0;
    for (i=0; i<ne; i++) { // for each element i      
        nfe = nfes[elementtype[i]]; 
        nve = nves[elementtype[i]]; 
        for (k=0; k<nfe; k++) { // for each face of element i          
            if ((t2t[i*nfemax+k] > i) || (t2t[i*nfemax+k] < 0)) {
                same = 0;
                for (j=i+1; j<ne; j++) { // for each element j      
                    nfe = nfes[elementtype[j]];                        
                    for (n=0; n<nfe; n++) // for each face of element j                                      
                        if (IsSameFace(faces[i][k],faces[j][n])) {
                            same = 1;
                            break;                
                        }
                    if (same==1)
                        break;
                }
                if (same==1) {// interior face
                    t2t[i*nfemax+k] = j;
                    t2t[j*nfemax+n] = i;                    
                    t2f[i*nfemax+k] = nf;
                    t2f[j*nfemax+n] = nf;                    
                    f[nf*nmax+nmax-2] = i;
                    f[nf*nmax+nmax-1] = j;
                }
                else {// boundary face
                    t2f[i*nfemax+k] = nf;    
                    f[nf*nmax+nmax-2] = i;
                }
                f[nf*nmax] = faces[i][k].size();
                for (m = 0; m<f[nf*nmax]; m++)                     
                    f[nf*nmax+m+1] =  t[i*nvemax+eleminfos[elementtype[i]].face[k][m]];                    
                    //f[nf*nmax+m+1] = faces[i][k][m];                                                
                nf = nf + 1;
            }
        }            
    }        
}

// void mkt2f(vector<Int> &t2f, vector<Int> &t2t, vector<Int> &f, vector<Int> &t, vector<Int> &elementtype, Int dim)
// {
//     
//     vector< elementinfo> eleminfos = mkeleminfos(dim);    
//     
//     Int i, j, k, m, n, nve, nfe, nvf, ne, same;
//     Int nvemax, nfemax, nvfmax, ntemax, elemtype, nelemtype;
//     ne = elementtype.size();        
//     nelemtype = eleminfos.size();   
//     
//     vector< Int> nfes(nelemtype,0), nves(nelemtype,0);    
//     for (i=0; i<nelemtype; i++) {
//         nfes[i] = eleminfos[i].nfe;
//         nves[i] = eleminfos[i].nve;
//     }        
//     nfemax = *max_element(nfes.begin(),nfes.end());    
//     ntemax = *max_element(elementtype.begin(),elementtype.end());    
//     nvemax = round(((double) t.size())/ne);            
//     
//     vector< Int > ti(nvemax, 0);        
//     vector< vector< vector<Int> > > faces(ne, vector< vector<Int> >());               
//     for (i=0; i<ne; i++)
//     {        
//         elemtype = elementtype[i];       
//         nve = nves[elemtype];        
//         
//         for (j=0; j<nve; j++)
//             ti[j] = t[i*nvemax+j];        
//         
//         faces[i] = mkf(eleminfos[elemtype].face, ti);                
//         
//         nfe = nfes[elemtype];
//         for(j = 0; j<nfe; j++)
//             sort(faces[i][j].begin(), faces[i][j].end());
//     }
//             
//     nvfmax = dim;
//     if ((ntemax>0) && (dim==3))
//         nvfmax = 4;            
//     
//     Int nmax = (1+nvfmax+2);
//     t2t.resize(nfemax*ne,-1);
//     t2f.resize(nfemax*ne,-1);
//     f.resize(nmax*nfemax*ne,-1);        
//     Int nf = 0;
//     for (i=0; i<ne; i++) { // for each element i      
//         nfe = nfes[elementtype[i]]; 
//         nve = nves[elementtype[i]]; 
//         for (k=0; k<nfe; k++) { // for each face of element i          
//             if ((t2t[i*nfemax+k] > i) || (t2t[i*nfemax+k] < 0)) {
//                 same = 0;
//                 for (j=i+1; j<ne; j++) { // for each element j      
//                     nfe = nfes[elementtype[j]];                        
//                     for (n=0; n<nfe; n++) // for each face of element j                                      
//                         if (IsSameFace(faces[i][k],faces[j][n])) {
//                             same = 1;
//                             break;                
//                         }
//                     if (same==1)
//                         break;
//                 }
//                 if (same==1) {// interior face
//                     t2t[i*nfemax+k] = j;
//                     t2t[j*nfemax+n] = i;                    
//                     t2f[i*nfemax+k] = nf;
//                     t2f[j*nfemax+n] = nf;                    
//                     f[nf*nmax+nmax-2] = i;
//                     f[nf*nmax+nmax-1] = j;
//                 }
//                 else {// boundary face
//                     t2f[i*nfemax+k] = nf;    
//                     f[nf*nmax+nmax-2] = i;
//                 }
//                 f[nf*nmax] = (Int) faces[i][k].size();
//                 for (m = 0; m<f[nf*nmax]; m++)                     
//                     f[nf*nmax+m+1] = t[i*nve+eleminfos[elementtype[i]].face[k][m]];                    
//                     //f[nf*nmax+m+1] = faces[i][k][m];                                                
//                 nf = nf + 1;
//             }
//         }            
//     }    
//     f.resize(nmax*nf);              
// }

