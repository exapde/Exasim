#ifndef __CPUIMPL_H__
#define __CPUIMPL_H__

void cpuArraySetValue(dstype *udg, dstype a, int N)
{
  for (int i=0; i<N; i++)
    udg[i] = a;
}

void cpuArrayInsert(dstype* u, const dstype* un, const int I, const int J, const int K, 
        const int i1, const int i2, const int j1, const int j2, const int k1, const int k2)
{        
    int ni = i2-i1;
    int nj = j2-j1;
    int nk = k2-k1;    
    int M = ni*nj;
    int N = M*nk;
    int Q = I*J;    
    for (int idx=0; idx<N; idx++) {
        int l = idx%M;        
        int i = l%ni+i1;
        int j = (l-i)/ni+j1;
        int k = (idx-l)/M+k1;
        u[i+I*j+Q*k] = un[idx];        
    }
}


dstype cpuArrayMin(dstype *a, int n)
{
    dstype b = a[0];
    for (int i=1; i<n; i++)
        if (a[i]<b)
            b = a[i];    
    return b;
}

dstype cpuArrayMax(dstype *a, int n)
{
    dstype b = a[0];
    for (int i=1; i<n; i++)
        if (a[i]>b)
            b = a[i];    
    return b;
}

void cpuElemGeom(dstype *Xx, dstype *jac, dstype *Jg, int ne, int ng, int nd)
{        
    dstype *Jg11, *Jg12, *Jg13, *Jg21, *Jg22, *Jg23, *Jg31, *Jg32, *Jg33;
    dstype *Xx11, *Xx12, *Xx13, *Xx21, *Xx22, *Xx23, *Xx31, *Xx32, *Xx33;
    int ngv = ng*ne;
     
    if (nd==1) {
        for (int i=0; i<ngv; i++) {
            jac[i] = Jg[i];
            Xx[i] = 1.0;
        }
    }
    else if (nd==2) {
        Jg11 = &Jg[0];
        Jg12 = &Jg[ngv];
        Jg21 = &Jg[2*ngv];
        Jg22 = &Jg[3*ngv];
        Xx11 = &Xx[0];
        Xx21 = &Xx[ngv];
        Xx12 = &Xx[2*ngv];
        Xx22 = &Xx[3*ngv];
        for (int i=0; i<ngv; i++) {
            jac[i] = Jg11[i]*Jg22[i] - Jg12[i]*Jg21[i];
            Xx11[i] = Jg22[i]; // dxi/dx
            Xx21[i] = -Jg21[i]; //dxi/dy 
            Xx12[i] = -Jg12[i]; //deta/dx
            Xx22[i] = Jg11[i]; //deta/dy
        }        
    }
    else if (nd==3) {
        Jg11 = &Jg[0];
        Jg12 = &Jg[ngv];
        Jg13 = &Jg[2*ngv];
        Jg21 = &Jg[3*ngv];
        Jg22 = &Jg[4*ngv];
        Jg23 = &Jg[5*ngv];
        Jg31 = &Jg[6*ngv];
        Jg32 = &Jg[7*ngv];
        Jg33 = &Jg[8*ngv];
        Xx11 = &Xx[0];
        Xx21 = &Xx[ngv];
        Xx31 = &Xx[2*ngv];
        Xx12 = &Xx[3*ngv];
        Xx22 = &Xx[4*ngv];
        Xx32 = &Xx[5*ngv];
        Xx13 = &Xx[6*ngv];
        Xx23 = &Xx[7*ngv];
        Xx33 = &Xx[8*ngv];
        for (int i=0; i<ngv; i++) {
            jac[i] = Jg11[i]*Jg22[i]*Jg33[i] - Jg11[i]*Jg32[i]*Jg23[i] +
                     Jg21[i]*Jg32[i]*Jg13[i] - Jg21[i]*Jg12[i]*Jg33[i] +
                     Jg31[i]*Jg12[i]*Jg23[i] - Jg31[i]*Jg22[i]*Jg13[i];
            Xx11[i] = Jg22[i]*Jg33[i] - Jg23[i]*Jg32[i];
            Xx21[i] = Jg23[i]*Jg31[i] - Jg21[i]*Jg33[i];
            Xx31[i] = Jg21[i]*Jg32[i] - Jg22[i]*Jg31[i];
            Xx12[i] = Jg13[i]*Jg32[i] - Jg12[i]*Jg33[i];
            Xx22[i] = Jg11[i]*Jg33[i] - Jg13[i]*Jg31[i];
            Xx32[i] = Jg12[i]*Jg31[i] - Jg11[i]*Jg32[i];
            Xx13[i] = Jg12[i]*Jg23[i] - Jg13[i]*Jg22[i];
            Xx23[i] = Jg13[i]*Jg21[i] - Jg11[i]*Jg23[i];
            Xx33[i] = Jg11[i]*Jg22[i] - Jg12[i]*Jg21[i];
        }        
    }
}

void cpuGetElemNodes(dstype* unView, const dstype* uView, const int np, const int nc, const int nc1, const int nc2, const int e1, const int e2) 
{
    int nn = np * (e2 - e1);
    int ncu = nc2 - nc1;
    int N = nn * ncu;
    int K = np * nc;

    for (int idx=0; idx<N; idx++) {
        int i = idx % nn;  // [0, np*(e2-e1)]
        int j = idx / nn;  // [0, ncu]
        int k = i % np;    // [0, np]
        int e = i / np + e1;
        unView[idx] = uView[k + (j + nc1) * np + e * K];
    }
}

void cpuApplyGivensRotation(dstype *H, dstype *s, dstype *cs, dstype *sn,  int i)
{        
    dstype temp;    
    for (int k=0; k<i; k++) {
        temp       =  cs[k]*H[k] + sn[k]*H[k+1];
        H[k+1] = -sn[k]*H[k] + cs[k]*H[k+1];
        H[k]   = temp;
    }

    if (H[i+1] == 0.0) {
        cs[i] = 1.0;
        sn[i] = 0.0;
    } 
    else if (fabs(H[i+1]) > fabs(H[i])) {
        temp = H[i] / H[i+1];
        sn[i] = 1.0 / sqrt( 1.0 + temp*temp );
        cs[i] = temp * sn[i];
    } 
    else {
        temp = H[i+1] / H[i];
        cs[i] = 1.0 / sqrt( 1.0 + temp*temp );
        sn[i] = temp * cs[i];
    }
   
    temp   = cs[i]*s[i];                       
    s[i+1] = -sn[i]*s[i];
    s[i]   = temp;
    H[i] = cs[i]*H[i] + sn[i]*H[i+1];
    H[i+1] = 0.0;
}

void cpuBackSolve(dstype *y, dstype *H, dstype *s, int i, int n)
{
    for (int j=i; j>=0; j--)
        y[j] = s[j];    
    
    for (int j=i; j>=0; j--) {
        y[j] =  y[j]/H[j+n*j]; 
        for (int k=j-1; k>=0; k--)
            y[k] = y[k] - H[k+n*j]*y[j]; 
    }
}

void faceindex(int *in1, int *in2, int *facecon, int npf, int ncu, int npe, int nc, int f1, int f2)
{    
    int nf = f2-f1;
    int ndf = npf*nf;
    int N = ndf*ncu;    
    for (int idx = 0; idx<N; idx++)
    {
        int i = idx%ndf;
        int j = (idx-i)/ndf;
        int m = npf*f1+i;
        int k1 = facecon[2*m];
        int k2 = facecon[2*m+1];
        int m1 = k1%npe;
        int m2 = k2%npe;
        int n1 = (k1-m1)/npe;
        int n2 = (k2-m2)/npe;          
        in1[idx] = m1+j*npe+n1*npe*nc;
        in2[idx] = m2+j*npe+n2*npe*nc;
    }                            
}

void faceperm(int *ind1, int *ind2, int *indpts, int *facecon, int *fblks, int npf, int ncu, int npe, int nc, int nbf)
{
    int N = 0;
    for (int j=0; j<nbf; j++) {
        int f1 = fblks[3*j]-1;
        int f2 = fblks[3*j+1];      
        int nf = f2-f1;
        int ndf = npf*nf;        
        faceindex(&ind1[N], &ind2[N], facecon, npf, ncu, npe, nc, f1, f2);
        indpts[j] = N;
        N = N + ndf*ncu;    
    }        
    indpts[nbf] = N;
}

void faceindex1(int *in1, int *facecon, int npf, int ncu, int npe, int nc, int f1, int f2)
{    
    int nf = f2-f1;
    int ndf = npf*nf;
    int N = ndf*ncu;    
    for (int idx = 0; idx<N; idx++)
    {
        int i = idx%ndf;
        int j = (idx-i)/ndf;
        int m = npf*f1+i;
        int k1 = facecon[2*m];
        int m1 = k1%npe;
        int n1 = (k1-m1)/npe;
        in1[idx] = m1+j*npe+n1*npe*nc;
    }                            
}

void faceperm1(int *ind1, int *indpts, int *facecon, int *fblks, int npf, int ncu, int npe, int nc, int nbf)
{
    int N = 0;
    for (int j=0; j<nbf; j++) {
        int f1 = fblks[3*j]-1;
        int f2 = fblks[3*j+1];      
        int nf = f2-f1;
        int ndf = npf*nf;        
        faceindex1(&ind1[N], facecon, npf, ncu, npe, nc, f1, f2);
        indpts[j] = N;
        N = N + ndf*ncu;    
    }        
    indpts[nbf] = N;
}

void elemindex(int *ind, int npe, int nc, int ncu, int e1, int e2)
{        
    int nn = npe*(e2-e1);
    int N = nn*ncu;
    for (int idx = 0; idx<N; idx++)
    {
        int i = idx%nn;   // [0, npe*ne]
        int j = (idx-i)/nn; // [0, ncu]
        int k = i%npe;  // [0, npe]
        int e = (i-k)/npe+e1;        
        ind[idx] = k+j*npe+e*npe*nc;        
    }        
}

void elemperm(int *ind, int *indpts, int *eblks, int npe, int nc, int ncu, int nbe)
{
    int N=0;
    for (int j=0; j<nbe; j++) {
        int e1 = eblks[3*j]-1;
        int e2 = eblks[3*j+1];    
        int nn = npe*(e2-e1);
        elemindex(&ind[N], npe, nc, ncu, e1, e2);
        indpts[j] = N;
        N = N + nn*ncu;
    }      
    indpts[nbe] = N;
}

void getboundaryfaces(int *numbf, int *boufaces, int *bf, int *eblks, int nbe, int nfe, int maxbc, int nboufaces)
{
    for (int i=0; i<1+maxbc*nbe; i++) numbf[i] = 0;
    for (int i=0; i<nboufaces; i++) boufaces[i] = 0;

    for (int j=0; j<nbe; j++) { // loop over each chunk
        int e1 = eblks[3*j]-1;
        int e2 = eblks[3*j+1];            
        for (int e=e1; e<e2; e++) {  // loop over each element in a chunk           
            for (int k=0; k<nfe; k++) { // loop over each local face
                int ib = bf[k + nfe*e];                
                if (ib > 0) numbf[ib + maxbc*j] += 1; // boundary face
            }
        }          
    }                     

    // accumulative sum of numbf
    for (int j=0; j<nbe; j++) { // loop over each chunk
        for (int k=0; k<maxbc; k++) { // loop over each boundary condition
            int n = k + maxbc*j;              
            numbf[n+1] += numbf[n];
        }
    }      
          
    for (int j=0; j<nbe; j++) { // loop over each chunk
        for (int k=0; k<maxbc; k++) { // loop over each boundary condition
            int n = k + maxbc*j;
            int start = numbf[n];
            int nfaces = numbf[n+1] - start; // number of boundary faces for condition k             
            if (nfaces > 0) { // if the number of faces is not zero
              int idx = 0;
              int e1 = eblks[3*j]-1;
              int e2 = eblks[3*j+1];            
              for (int e=e1; e<e2; e++) {  // loop over each element in a chunk                    
                  for (int l=0; l<nfe; l++) { // loop over each local face
                      int ib = bf[l + nfe*e];  // boundary condition ib
                      if (ib == k+1) { // if boundary condition ib match 
                        boufaces[start + idx] = l + nfe*(e-e1); // store the boundary face (l,e)
                        idx++;
                      }
                  }
              }                                  
            }
        }
    }      
}

int getinterfacefaces(const int* f2e, int ne1, int nf)
{
    int n = 0; 
    for (int j=0; j<nf; j++) { // loop over each face
        int e1 = f2e[4*j+0];   // 1st element sharing the face 
        int e2 = f2e[4*j+2];   // 2nd element sharing the face  
        if ((e1 >= ne1) || (e2 >= ne1)) { // interface faces
          n += 1;
        }
    }    
    return n;
}

int getinterfacefaces(int *interface, const int* f2e, int ne1, int nf)
{
    int n = 0; 
    for (int j=0; j<nf; j++) { // loop over each face
        int e1 = f2e[4*j+0];   // 1st element sharing the face 
        int e2 = f2e[4*j+2];   // 2nd element sharing the face  
        if ((e1 >= ne1) || (e2 >= ne1)) { // interface faces
          interface[n] = j;
          n += 1;
        }
    }    
    return n;
}

#endif  

