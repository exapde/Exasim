#ifndef __GETPOLY
#define __GETPOLY

void MGS(cublasHandle_t handle, dstype *V, dstype *H, Int N, Int m, Int backend)
{
    for (Int k = 0; k < m; k++) {
        PDOT(handle, N, &V[m*N], inc1, &V[k*N], inc1, &H[k], backend);
        ArrayAXPY(handle, &V[m*N], &V[k*N], -H[k], N, backend);
    }
    PDOT(handle, N, &V[m*N], inc1, &V[m*N], inc1, &H[m], backend);
    H[m] = sqrt(H[m]);    
    ArrayMultiplyScalar(handle, &V[m*N], one/H[m], N, backend);
}

void MGS(cublasHandle_t handle, dstype *V, dstype *H, Int N, Int m, Int L, Int backend)
{
    if (L==0) {
      MGS(handle, V, H, N, m, backend);
      return;
    }
      
    double h1, h2;
    Int M = N-L;
    for (Int k = 0; k < m; k++) {
        PDOT(handle, L, &V[m*N], inc1, &V[k*N], inc1, &h1, backend);
        PDOT(handle, M, &V[m*N+L], inc1, &V[k*N+L], inc1, &h2, backend);
        H[k] = h2 + 0.5*h1;
        ArrayAXPY(handle, &V[m*N], &V[k*N], -H[k], N, backend);
    }
    PDOT(handle, L, &V[m*N], inc1, &V[m*N], inc1, &h1, backend);
    PDOT(handle, M, &V[m*N+L], inc1, &V[m*N+L], inc1, &h2, backend);
    H[m] = h2 + 0.5*h1;
    H[m] = sqrt(H[m]);    
    ArrayMultiplyScalar(handle, &V[m*N], one/H[m], N, backend);
}

void makeH(CDiscretization &disc, CPreconditioner& prec, sysstruct &sys, 
        dstype *H, dstype *r, Int N, Int m, Int backend)
{
    int m1 = m + 1;
    for(int i=0; i<m*m1; i++)
        H[i] = 0.0;
    
    dstype normr = PNORM(disc.common.cublasHandle, N, r, backend);    
          
    ArrayAX(disc.common.cublasHandle, sys.v, r, 1.0/normr, N, backend);    
    for (int j=0; j<m; j++) {                    
        disc.evalMatVec(&sys.v[(j+1)*N], &sys.v[j*N], sys.u, sys.b, backend);      
        prec.ApplyPreconditioner(&sys.v[(j+1)*N], sys, disc, backend);            
        MGS(disc.common.cublasHandle, sys.v, &H[m1*j], N, j+1, backend);
    }
}

void ArrayTranspose(dstype  *A, dstype  *B, int m, int n, int m1)
{  
    for (int i=0; i<m; i++)
        for (int j=0; j<n; j++)
            A[j+n*i] = B[i + m1*j];     
}

void MatMul(dstype *c, dstype *a, dstype *b, int r1, int c1, int c2) 
{
    int i, j, k;

    for(j = 0; j < c2; j++)        
        for(i = 0; i < r1; i++)
            c[i + r1*j] = 0.0;        
    
    for(j = 0; j < c2; j++)
        for(i = 0; i < r1; i++)        
            for(k = 0; k < c1; k++)            
                c[i + r1*j] += a[i + r1*k] * b[k + c1*j];            
}

void eigH(dstype  *HT, dstype  *H, dstype *em, dstype *work, int *ipiv, int m)
{
    int m1 = m+1;
    
    ArrayTranspose(HT, H, m, m, m1);    
    cpuComputeInverse(HT, work, ipiv, m);
    MatMul(work, HT, em, m, m, 1); 
            
    dstype Beta2 = H[m1*m - 1]*H[m1*m - 1];
    for (int j=0; j<m; j++)
        for (int i=0; i<m; i++)        
            HT[i + m*j] = H[i + m1*j] + Beta2*work[i]*em[j]; 
}

dstype complexabs(dstype a, dstype b) {
    return sqrt(a*a + b*b);
}

void LejaSort(dstype *sr, dstype *si, dstype *lr, dstype *li, dstype *product, int n)
{
    sr[0] = lr[0];
    si[0] = li[0];
    
    int j = 1;
    for (int i=1; i<n; i++) {
        if (complexabs(lr[i], li[i]) > complexabs(sr[0], si[0])) {
            sr[0] = lr[i];
            si[0] = li[i];
        }
    }
        
    if (fabs(si[0]) > 1e-10) {
        sr[1] = sr[0];
        si[1] = -si[0];
        j = 2;
    }
    else 
        si[0] = 0.0;
    
    dstype productj;    
    while (j < n) {
        for (int i=0; i<n; i++) {
            product[i] = 1;
            for (int ii=0; ii<j; ii++) {
                product[i] = product[i] + log(complexabs(lr[i]-sr[ii], li[i]-si[ii]));
            }
        }
        sr[j] = lr[0];
        si[j] = li[0];
        productj = product[0];
        for (int i=1; i<n; i++) {
            if( product[i] > productj ) {
                sr[j] = lr[i];
                si[j] = li[i];
                productj = product[i];
            }
        }
    
        if (fabs(si[j]) > 1e-10) {
            sr[j+1] = sr[j];
            si[j+1] = -si[j];
            j += 1;
        }                
        else
            si[j] = 0.0;
            
        j += 1;        
    }            
}

void getPoly(CDiscretization &disc, CPreconditioner& prec, sysstruct &sys, 
        dstype  *lam, dstype *r, int *ipiv, int N, int m, int backend)
{
    dstype *Hm = &lam[0];    
    dstype *Ht = &lam[m+m*m];      
    dstype *work = &lam[m+2*m*m];    
    dstype *em = &lam[2*m+2*m*m];        
    dstype *wr = &lam[0];
    dstype *wi = &lam[m];    
    dstype *lamr = &lam[2*m];
    dstype *lami = &lam[3*m];    
    
    int info;
    int lwork = 5*m;
    char chv = 'N';    
    
    for (int i=0; i < m; i++)
        em[i] = 0.0;
    em[m-1] = 1.0;           
    
    makeH(disc, prec, sys, Hm, r, N, m, backend);
    eigH(Ht, Hm, em, work, ipiv, m);        
    DGEEV(&chv, &chv, &m, Ht, &m, wr, wi, work, &m, work, &m, work, &lwork, &info );               
    LejaSort(lamr, lami, wr, wi, work, m);
}

void makeH(CDiscretization &disc, CPreconditioner& prec, sysstruct &sys, 
        dstype *H, dstype *r, Int N, Int m, Int spatialScheme, Int backend)
{
    int m1 = m + 1;
    for(int i=0; i<m*m1; i++)
        H[i] = 0.0;
    
    dstype normr = PNORM(disc.common.cublasHandle, N, disc.common.ndofuhatinterface, r, backend);        
    ArrayAX(disc.common.cublasHandle, sys.v, r, 1.0/normr, N, backend);    
    for (int j=0; j<m; j++) {                    
        disc.evalMatVec(&sys.v[(j+1)*N], &sys.v[j*N], sys.u, sys.b, spatialScheme, backend);      
        prec.ApplyPreconditioner(&sys.v[(j+1)*N], sys, disc, spatialScheme, backend);            
        MGS(disc.common.cublasHandle, sys.v, &H[m1*j], N, j+1, disc.common.ndofuhatinterface, backend);
    }
}

void getPoly(CDiscretization &disc, CPreconditioner& prec, sysstruct &sys, 
        dstype  *lam, dstype *r, int *ipiv, int N, int m, int spatialScheme, int backend)
{
    dstype *Hm = &lam[0];    
    dstype *Ht = &lam[m+m*m];      
    dstype *work = &lam[m+2*m*m];    
    dstype *em = &lam[2*m+2*m*m];        
    dstype *wr = &lam[0];
    dstype *wi = &lam[m];    
    dstype *lamr = &lam[2*m];
    dstype *lami = &lam[3*m];    
    
    int info;
    int lwork = 5*m;
    char chv = 'N';    
    
    for (int i=0; i < m; i++)
        em[i] = 0.0;
    em[m-1] = 1.0;                   
        
    makeH(disc, prec, sys, Hm, r, N, m, spatialScheme, backend);        
    eigH(Ht, Hm, em, work, ipiv, m);                
    DGEEV(&chv, &chv, &m, Ht, &m, wr, wi, work, &m, work, &m, work, &lwork, &info );            
    LejaSort(lamr, lami, wr, wi, work, m);    
}
    
#endif

