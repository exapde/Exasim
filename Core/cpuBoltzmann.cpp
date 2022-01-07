#ifndef __CPUBOLTZMANN
#define __CPUBOLTZMANN

template <typename T> void cpugetxg(T* xb, T* xdg, int* fb, int nbf, int ngf, int nfe, int ne, int ncx)
{    
    int M = nbf*ngf;
    int N = nbf*ngf*ncx;
    #pragma omp parallel for
    for (int idx=0; idx<N; idx++){
        int m = idx%M;     // [1, nbf*ngf]
        int i = (idx-m)/M; // [1, ncx]
        int j = m%ngf;     // [1, ngf]
        int k = (m-j)/ngf; // [1, nbf]
        int e = fb[0+4*k];
        int f = fb[1+4*k];
        xb[idx] = xdg[j+ngf*f+ngf*nfe*e+ngf*nfe*ne*i];               
    }        
}
template void cpugetxg(double*, double*, int*, int, int, int, int, int);
template void cpugetxg(float*, float*, int*, int, int, int, int, int);

template <typename T> void cpugetug(T* ug, T* udg, int* fb, int* perm, int nbf, int npf, int npe, int ne, int ncu)
{    
    int M = nbf*npf;
    int N = nbf*npf*ncu;
    #pragma omp parallel for
    for (int idx=0; idx<N; idx++){
        int m = idx%M;     // [1, nbf*npf]
        int i = (idx-m)/M; // [1, ncu]
        int j = m%npf;     // [1, npf]
        int k = (m-j)/npf; // [1, nbf]
        int e = fb[0+4*k];
        int f = fb[1+4*k];
        ug[idx] = udg[perm[j+npf*f]+npe*e+npe*ne*i];         
    }        
}
template void cpugetug(double*, double*, int*, int*, int, int, int, int, int);
template void cpugetug(float*, float*, int*, int*, int, int, int, int, int);

template <typename T> void cpugetfg(T* fg, T* fdg, int* fb, int* perm, int nbf, int npf, int npe, int ne)
{    
    int N = nbf*npf;
    #pragma omp parallel for
    for (int idx=0; idx<N; idx++){
        int j = idx%npf;     // [1, npf]
        int k = (idx-j)/npf; // [1, nbf]        
        int e = fb[0+4*k];
        int f = fb[1+4*k];
        int n = fb[2+4*k];
        fg[idx] = fdg[perm[j+npf*f]+npe*e+npe*ne*n];
    }    
}
template void cpugetfg(double*, double*, int*, int*, int, int, int, int);
template void cpugetfg(float*, float*, int*, int*, int, int, int, int);

template <typename T> void cpuputfg(T* fdg, T* fg, T scalar, int* fc, int M, int N)
{    
    #pragma omp parallel for
    for (int i=0; i<M; i++) {   
        int K = fc[(N+1)+i+1]-fc[(N+1)+i]; 
        for (int j=0; j<K; j++)
            fdg[fc[i]] += scalar*fg[fc[2*(N+1)+fc[(N+1)+i]+j]];
    }            
}
template void cpuputfg(double*, double*, double, int*, int, int);
template void cpuputfg(float*, float*, float, int*, int, int);

template <typename T> void cpugetc(T* c, T* pv, int* fb, int nbf, int nd)
{    
    int N = nbf*nd;
    #pragma omp parallel for
    for (int idx=0; idx<N; idx++){
        int j = idx%nd;     // [1, nd]
        int k = (idx-j)/nd; // [1, nbf]        
        int n = fb[2+4*k];
        c[j+nd*k] = pv[j+nd*n]; 
    }           
}
template void cpugetc(double*, double*, int*, int, int);
template void cpugetc(float*, float*, int*, int, int);

template <typename T> void cpugetnjc(T *njc, T* DetJf, T* nlf, T *c, int ngf, int nbf, int nd)
{
    int N = nbf*ngf;
    #pragma omp parallel for
    for (int idx=0; idx<N; idx++){
        int j = idx%ngf;     // [1, ngf]
        int i = (idx-j)/ngf; // [1, nbf]      
        njc[idx] = 0.0;
        for (int k=0; k<nd; k++) 
            njc[idx] += DetJf[idx]*nlf[j+ngf*i+ngf*nbf*k]*c[k+nd*i];        
    }        
}
template void cpugetnjc(double*, double*, double*, double*, int, int, int);
template void cpugetnjc(float*, float*, float*, float*, int, int, int);

template <typename T> void cpuapplyinverse(T* F, T* R, T *f, T* Aeinv, T* Beinv, int* t2t,
        int* nvsweep, int* evsweep1, int* evsweep2, int* linflowi,  int* ninflowi,  
        int* finflowi, int npe, int nfe, int ne, int nev, int s)
{    
    //int nev = nvsweep[s+1]-nvsweep[s];
    #pragma omp parallel for
    for (int i=0; i<nev; i++) {
        int e = evsweep1[nvsweep[s]+i];
        int n = evsweep2[nvsweep[s]+i];
        for (int j=0; j<npe; j++)
            f[j] = R[npe*ne*n+npe*e+j];

        int k = e+(ne+1)*n;            
        int n1 = linflowi[n];                      
        int n2 = npe*npe*(n1+ninflowi[k]); 
        int nfaces = ninflowi[k+1]-ninflowi[k];
        for (int m=0; m<nfaces; m++) {
            int l = finflowi[n1+ninflowi[k]+m];
            int e1 = t2t[l+nfe*e];
            for (int j=0; j<npe; j++)
                for (int p=0; p<npe; p++)
                    f[j] = f[j]-Beinv[n2+npe*npe*m+npe*p+j]*F[p+npe*e1+npe*ne*n];                
        }
        
        for (int j=0; j<npe; j++) {
            F[npe*ne*n+npe*e+j] = 0.0;
            for (int p=0; p<npe; p++)
                F[npe*ne*n+npe*e+j] += Aeinv[npe*npe*ne*n+npe*npe*e+npe*p+j]*f[p];            
        }        
    }                
}
template void cpuapplyinverse(double*, double*, double*, double*, double*, 
        int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);
template void cpuapplyinverse(float*, float*, float*, float*, float*, 
        int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);

#endif