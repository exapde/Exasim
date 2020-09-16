#ifndef __CPUARRAYPERMUTE
#define __CPUARRAYPERMUTE

template <typename T> void cpuKron(T *C, T *A, T *B, int M1, int N1, int M2, int N2)
{            
    int M = M1*M2;
    int N = N1*N2;
    int P = M*N;    
    #pragma omp parallel for
    for (int idx=0; idx<P; idx++)     
    {
        int i = idx%M;
        int j = (idx-i)/M;
        int ib = i%M2;
        int ia = (i-ib)/M2;
        int jb = j%N2;
        int ja = (j-jb)/N2;
        C[idx] = A[ia+M1*ja]*B[ib+M2*jb];        
    }
}

void cpuIndexPermute12(int *index, int I1, int I2, int I3)
{                     
    int M = I1*I2;
    int N = M*I3;
    #pragma omp parallel for
    for (int idx=0; idx<N; idx++)     
    {
        int l = idx%M;
        int i = l%I1;
        int j = (l-i)/I1;
        int k = (idx-l)/M;        
        index[idx] = j+I2*i+M*k;
    }
}

void cpuIndexPermute13(int *index, int I1, int I2, int I3, int I4)
{                
    int M = I1*I2;
    int N = M*I3;
    int P = N*I4;
    #pragma omp parallel for
    for (int idx=0; idx<P; idx++)     
    {
        int n = idx%N;
        int i4 = (idx-n)/N;        
        int m = n%M;
        int i3 = (n-m)/M;      
        int i1 = m%I1;
        int i2 = (m-i1)/I1;      
        index[idx] = i3+I3*i2+I3*I2*i1+N*i4;
    }
}

void cpuIndexPermute23(int *index, int I1, int I2, int I3, int I4)
{                
    int M = I1*I2;
    int N = M*I3;
    int P = N*I4;
    #pragma omp parallel for
    for (int idx=0; idx<P; idx++)     
    {
        int n = idx%N;
        int i4 = (idx-n)/N;        
        int m = n%M;
        int i3 = (n-m)/M;      
        int i1 = m%I1;
        int i2 = (m-i1)/I1;      
        index[idx] = i1+I1*i3+I1*I3*i2+N*i4;        
    }
}


template <typename T> void cpuPermute(T *B, T *A, int *index, int N)
{    
    #pragma omp parallel for
    for (int idx=0; idx<N; idx++)          
        B[index[idx]] = A[idx];            
}

template <typename T> void cpuPermuteSharedMem(T *B, T *A, int *index, int N, int BLOCKDIM)
{
    T Ashared[BLOCKDIM];     
    #pragma omp parallel for
    for (int idx=0; idx<N; idx++)          
    {
        int tx = idx%BLOCKDIM;
        Ashared[index[tx]] = A[idx];         
        B[idx] = Ashared[tx];        
    }
}

template <typename T> void cpuPermute12(T *B, T *A, int I1, int I2, int I3)
{            
    int M = I1*I2;
    int N = M*I3;    
    #pragma omp parallel for
    for (int idx=0; idx<N; idx++)          
    {
        int l = idx%M;
        int i = l%I1;
        int j = (l-i)/I1;
        int k = (idx-l)/M;        
        B[j+I2*i+M*k] = A[idx];
    }
}

template <typename T> void cpuPermute13(T *B, T *A, int I1, int I2, int I3, int I4)
{                
    int M = I1*I2;
    int N = M*I3;
    int P = N*I4;
    #pragma omp parallel for
    for (int idx=0; idx<P; idx++) 
    {
        int n = idx%N;
        int i4 = (idx-n)/N;        
        int m = n%M;
        int i3 = (n-m)/M;      
        int i1 = m%I1;
        int i2 = (m-i1)/I1;      
        B[i3+I3*i2+I3*I2*i1+N*i4] = A[idx];        
    }
}

template <typename T> void cpuPermute23(T *B, T *A, int I1, int I2, int I3, int I4)
{                
    int M = I1*I2;
    int N = M*I3;
    int P = N*I4;
    #pragma omp parallel for
    for (int idx=0; idx<P; idx++)  
    {
        int n = idx%N;
        int i4 = (idx-n)/N;        
        int m = n%M;
        int i3 = (n-m)/M;      
        int i1 = m%I1;
        int i2 = (m-i1)/I1;              
        B[i1+I1*i3+I1*I3*i2+N*i4] = A[idx];        
    }
}


//void cpuIndexPermute12(int*, int, int, int);
//void cpuIndexPermute13(int*, int, int, int, int);
//void cpuIndexPermute23(int*, int, int, int, int);

template void cpuKron(double*, double*, double*, int, int, int, int);
template void cpuKron(float*, float*, float*, int, int, int, int);

template void cpuPermute(double*, double*, int*, int);
template void cpuPermute(float*, float*, int*, int);
template void cpuPermuteSharedMem(double*, double*, int*, int, int);
template void cpuPermuteSharedMem(float*, float*, int*, int, int);

template void cpuPermute12(double*, double*, int, int, int);
template void cpuPermute12(float*, float*, int, int, int);
template void cpuPermute13(double*, double*, int, int, int, int);
template void cpuPermute13(float*, float*, int, int, int, int);
template void cpuPermute23(double*, double*, int, int, int, int);
template void cpuPermute23(float*, float*, int, int, int, int);

#endif
