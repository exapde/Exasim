#ifndef __CPUELEMFACENODE
#define __CPUELEMFACENODE

template <typename T> void cpuGetElemFaceNodes(T *un, T *u, int *perm, int npe, int npf, int nfe, int nc, int e1, int e2)
{        
    int ndf = npf*nfe;
    int nn = ndf*(e2-e1);    
    int N = nn*nc;
    #pragma omp parallel for
    for (int idx = 0; idx<N; idx++)
    {
        int i = idx%nn; // [0, npf*nfe*(e2-e1)]
        int j = (idx-i)/nn; // [0, nc]
        int k = i%ndf;   //[0, ndf]
        int e = (i-k)/ndf+e1; // [e1, e2]
        int m = k%npf; //[0, npf]
        int n = (k-m)/npf; // [0, nfe]
        un[idx] = u[perm[m+n*npf]+j*npe+e*npe*nc];        
    }        
    // un: npf*nfe*ne*nc
}

template <typename T> void cpuGetElemNodes(T *un, T *u, int np, int nc, int nc1, int nc2, int e1, int e2)
{        
    int nn = np*(e2-e1);
    int ncu = nc2-nc1;
    int N = nn*ncu;
    #pragma omp parallel for
    for (int idx = 0; idx<N; idx++)
    {
        int i = idx%nn;
        int j = (idx-i)/nn;
        int k = i%np;
        int e = (i-k)/np+e1;        
        un[idx] = u[k+(j+nc1)*np+e*np*nc];        
    }        
}

template <typename T> void cpuPutElemNodes(T *u, T *un, int np, int nc, int nc1, int nc2, int e1, int e2)
{
    int nn = np*(e2-e1);
    int ncu = nc2-nc1;
    int N = nn*ncu;
    #pragma omp parallel for
    for (int idx = 0; idx<N; idx++)
    {
        int i = idx%nn;
        int j = (idx-i)/nn;
        int k = i%np;
        int e = (i-k)/np+e1;        
        u[k+(j+nc1)*np+e*np*nc] = un[idx];        
    }            
}

template <typename T> void cpuGetElemNodes2(T *un, T *u, int np, int nc, int nc1, int nc2, int e1, int e2)
{        
    int ncu = nc2-nc1;
    int M = np*ncu;
    int N = np*ncu*(e2-e1);
    #pragma omp parallel for
    for (int idx=0; idx<N; idx++)
    {
        int i = idx%M;        
        int k = i%np;
        int j = (i-k)/np+nc1;
        int e = (idx-i)/M+e1;
        un[idx] = u[k+j*np+e*np*nc];        
    }            
}

template <typename T> void cpuPutElemNodes2(T *u, T *un, int np, int nc, int nc1, int nc2, int e1, int e2)
{        
    int ncu = nc2-nc1;
    int M = np*ncu;
    int N = np*ncu*(e2-e1);
    #pragma omp parallel for
    for (int idx=0; idx<N; idx++)
    {
        int i = idx%M;        
        int k = i%np;
        int j = (i-k)/np+nc1;
        int e = (idx-i)/M+e1;
        u[k+j*np+e*np*nc] = un[idx];        
    }            
}

template <typename T> void cpuGetFaceNodes(T *uh, T *udg, int *facecon, int npf, int ncu, int npe, int nc, int f1, int f2, int opts)
{    
    int nf = f2-f1;
    int ndf = npf*nf;
    int N = ndf*ncu;
    
    if (opts==0) {
        #pragma omp parallel for
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
            uh[idx] = 0.5*(udg[m1+j*npe+n1*npe*nc]+udg[m2+j*npe+n2*npe*nc]);
        }                        
    }
    else if (opts==1) {
        #pragma omp parallel for
        for (int idx = 0; idx<N; idx++)    
        {
            int i = idx%ndf;
            int j = (idx-i)/ndf;
            int m = npf*f1+i;
            int k1 = facecon[2*m];
            int m1 = k1%npe;
            int n1 = (k1-m1)/npe;
            uh[idx] = udg[m1+j*npe+n1*npe*nc];
        }                        
    }
    else if (opts==2) {
        #pragma omp parallel for
        for (int idx = 0; idx<N; idx++)       
        {
            int i = idx%ndf;
            int j = (idx-i)/ndf;
            int m = npf*f1+i;
            int k2 = facecon[2*m+1];
            int m2 = k2%npe;
            int n2 = (k2-m2)/npe;            
            uh[idx] = udg[m2+j*npe+n2*npe*nc];
        }                        
    }
}


template <typename T> void cpuPutFaceNodes(T *udg, T *uh, int *facecon, int npf, int ncu, int npe, int nc, int f1, int f2, int opts)
{    
    int nf = f2-f1;
    int ndf = npf*nf;
    int N = ndf*ncu;
    
    if (opts==0) {
        #pragma omp parallel for        
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
            udg[m1+j*npe+n1*npe*nc] = udg[m1+j*npe+n1*npe*nc] - uh[idx];
            udg[m2+j*npe+n2*npe*nc] = udg[m2+j*npe+n2*npe*nc] + uh[idx];            
        }                        
    }
    else {
        #pragma omp parallel for
        for (int idx = 0; idx<N; idx++)    
        {
            int i = idx%ndf;
            int j = (idx-i)/ndf;
            int m = npf*f1+i;
            int k1 = facecon[2*m];
            int m1 = k1%npe;
            int n1 = (k1-m1)/npe;            
            udg[m1+j*npe+n1*npe*nc] = udg[m1+j*npe+n1*npe*nc] - uh[idx];
        }                        
    }
}

template <typename T> void cpuPutFaceNodes(T *udg, T *uh, int *rowe2f1, int *cole2f1, int *ent2ind1,
        int *rowe2f2, int *cole2f2, int *ent2ind2, int npf, int npe, int nc, int e1, int e2, int opts)
{    
    int ne = e2-e1;
    int K = npf*nc;
    int M = npe*nc;
    int N = M*ne;        
    int I = M*e1;
    if (opts==0) {
        #pragma omp parallel for
        for (int idx = 0; idx<N; idx++) {
            int j = idx%M;              //[1, npe*nc]
            int k = (idx-j)/M+e1;       //[1, ne]      
            int l = j%npe;              //[1, npe]
            int m = (j-l)/npe;          //[1, nc] 
            int q, p, s;
            
            int i = ent2ind1[l+npe*k];
            int e = (i > 0) ? i : 0;
            int n = rowe2f1[i+1] - rowe2f1[e];            
            for (j=0; j<n; j++) {
                q = cole2f1[rowe2f1[i]+j];
                p = q%npf;              // [1, npf]
                s = (q-p)/npf;          // [1, nf]    
                udg[I+idx] = udg[I+idx] - uh[p+npf*m+K*s]; 
            }            
            
            i = ent2ind2[l+npe*k];
            e = (i > 0) ? i : 0;
            n = rowe2f2[i+1] - rowe2f2[e];
            for (j=0; j<n; j++) {
                q = cole2f2[rowe2f2[i]+j];
                p = q%npf;              // [1, npf]
                s = (q-p)/npf;          // [1, nf]      
                udg[I+idx] = udg[I+idx] + uh[p+npf*m+K*s]; 
            }                        
        }        
    }
    else {
        #pragma omp parallel for
        for (int idx = 0; idx<N; idx++) {
            int j = idx%M;              //[1, npe*nc]
            int k = (idx-j)/M+e1;       //[1, ne]      
            int l = j%npe;              //[1, npe]
            int m = (j-l)/npe;          //[1, nc] 
            
            int i = ent2ind1[l+npe*k];
            int e = (i > 0) ? i : 0;
            int n = rowe2f1[i+1] - rowe2f1[e];
            for (j=0; j<n; j++) {
                int q = cole2f1[rowe2f1[i]+j];
                int p = q%npf;          // [1, npf]
                int s = (q-p)/npf;          // [1, nf]    
                udg[I+idx] = udg[I+idx] - uh[p+npf*m+K*s]; 
            }            
        }
    }
}

template void cpuGetElemFaceNodes(double*, double*, int*, int, int, int, int, int, int);
template void cpuGetElemNodes(double*, double*, int, int, int, int, int, int);
template void cpuPutElemNodes(double*, double*, int, int, int, int, int, int);
template void cpuGetElemNodes2(double*, double*, int, int, int, int, int, int);
template void cpuPutElemNodes2(double*, double*, int, int, int, int, int, int);
template void cpuGetFaceNodes(double*, double*, int*, int, int, int, int, int, int, int);
template void cpuPutFaceNodes(double*, double*, int*, int, int, int, int, int, int, int);
template void cpuPutFaceNodes(double*, double*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);

template void cpuGetElemFaceNodes(float*, float*, int*, int, int, int, int, int, int);
template void cpuGetElemNodes(float*, float*, int, int, int, int, int, int);
template void cpuPutElemNodes(float*, float*, int, int, int, int, int, int);
template void cpuGetElemNodes2(float*, float*, int, int, int, int, int, int);
template void cpuPutElemNodes2(float*, float*, int, int, int, int, int, int);
template void cpuGetFaceNodes(float*, float*, int*, int, int, int, int, int, int, int);
template void cpuPutFaceNodes(float*, float*, int*, int, int, int, int, int, int, int);
template void cpuPutFaceNodes(float*, float*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);

#endif


