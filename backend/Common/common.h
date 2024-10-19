#ifndef __COMMON_H__
#define __COMMON_H__

#define SCOPY scopy_
#define SSCAL sscal_
#define SAXPY saxpy_
#define SDOT sdot_
#define SGEMV sgemv_
#define SGEMM sgemm_
#define SGETRF sgetrf_
#define SGETRI sgetri_
#define SGEEV sgeev_

#define DCOPY dcopy_
#define DSCAL dscal_
#define DAXPY daxpy_
#define DDOT ddot_
#define DGEMV dgemv_
#define DGEMM dgemm_
#define DGEEV dgeev_
#define DGETRF dgetrf_
#define DGETRI dgetri_

#ifdef USE_FLOAT
typedef float dstype;
#else
typedef double dstype; //  double is default precision 
#endif

#ifdef USE_LONG
typedef long Int;
#else
typedef int Int; 
#endif

#ifndef HAVE_CUDA    
#define cublasHandle_t int
#define cudaEvent_t int
#endif

// #ifdef HAVE_ENZYME                
// template <typename... Args>
// void __enzyme_autodiff(void*, Args... args);
// void __enzyme_fwddiff(void*, ...);
// int enzyme_const, enzyme_dup;
// #endif

typedef Kokkos::View<int*, Kokkos::HostSpace> view_1ih;
typedef Kokkos::View<dstype*, Kokkos::HostSpace> view_1dh;
typedef Kokkos::View<int*> view_1i;
typedef Kokkos::View<dstype*> view_1d;

#ifdef HAVE_MPP
#include <mutation++.h>
#endif

#define MKL_INT int

#define CPUFREE(x)                                                           \
{                                                                         \
    if (x != nullptr) {                                                      \
        free(x);                                                          \
        x = nullptr;                                                         \
    }                                                                     \
}

extern "C" {
    double DNRM2(Int*,double*,Int*);
    double DDOT(Int*,double*,Int*,double*,Int*);
    void DCOPY(Int*,double*,Int*,double*,Int*);    
    void DSCAL(Int*,double*,double*,Int*);
    void DAXPY(Int*,double*,double*,Int*,double*,Int*);
    void DGEMV(char*,Int*,Int*,double*,double*,Int*,double*,Int*,double*,double*,Int*);  
    void DGEMM(char*,char*,Int*,Int*,Int*,double*,double*,Int*,
             double*,Int*,double*,double*,Int*);        
    void DGETRF(Int*,Int*,double*,Int*,Int*,Int*);
    void DGETRI(Int*,double*,Int*,Int*,double*,Int*,Int*);
    void DTRSM(char *, char*, char*, char *, Int *, Int *, double*, double*, Int*,
             double*, Int*);
    void DGEEV( char* jobvl, char* jobvr, int* n, double* a,
                int* lda, double* wr, double* wi, double* vl, int* ldvl,
                double* vr, int* ldvr, double* work, int* lwork, int* info );    
    
    float SNRM2(Int*,float*,Int*);  
    float SDOT(Int*,float*,Int*,float*,Int*);
    void SCOPY(Int*,float*,Int*,float*,Int*);
    void SSCAL(Int*,double*,double*,Int*);
    void SAXPY(Int*,float*,float*,Int*,float*,Int*);
    void SGEMM(char*,char*,Int*,Int*,Int*,float*,float*,Int*,
             float*,Int*,float*,float*,Int*);  
    void SGEMV(char*,Int*,Int*,float*,float*,Int*,float*,Int*,float*,float*,Int*);      
    void SGETRF(Int*,Int*,float*,Int*,Int*,Int*);    
    void SGETRI(Int*,float*,Int*,Int*,float*,Int*,Int*);
    void STRSM(char *, char*, char*, char *, Int *, Int*, float*, float*, Int*,
             float*, Int*);        
    void SGEEV( char* jobvl, char* jobvr, int* n, double* a,
                int* lda, double* wr, double* wi, double* vl, int* ldvl,
                double* vr, int* ldvr, double* work, int* lwork, int* info );    
}

// global variables for BLAS  
dstype one = 1.0;
dstype minusone = -1.0;
dstype zero = 0.0;
char chn = 'N';
char cht = 'T';
char chl = 'L';
char chu = 'U';
char chr = 'R';
char chv = 'V';
Int inc1 = 1;

// global variables for CUBLAS  
// dstype *cublasOne;
// dstype *cublasMinusone;
// dstype *cublasZero;
dstype cublasOne[1] = {one};
dstype cublasMinusone[1] = {minusone};
dstype cublasZero[1] = {zero};

#ifdef HAVE_CUDA       
   #define CUDA_SYNC cudaDeviceSynchronize();  
#else 
   #define CUDA_SYNC
#endif                      

#ifdef TIMING    
    #define INIT_TIMING auto begin = chrono::high_resolution_clock::now(); auto end = chrono::high_resolution_clock::now();
#else
    #define INIT_TIMING
#endif

#ifdef TIMING
   #define TIMING_START  begin = chrono::high_resolution_clock::now();   
#else 
   #define TIMING_START     
#endif       

#ifdef TIMING
   #define TIMING_END    end = chrono::high_resolution_clock::now();   
#else 
   #define TIMING_END     
#endif       

#ifdef TIMING       
   #define TIMING_GET(num) common.timing[num] += chrono::duration_cast<chrono::nanoseconds>(end-begin).count()/1e6;        
#else 
   #define TIMING_GET(num)  
#endif                      

#ifdef TIMING
   #define START_TIMING {CUDA_SYNC; TIMING_START;}       
#else 
   #define START_TIMING
#endif       

#ifdef TIMING
   #define END_TIMING(num) {CUDA_SYNC; TIMING_END; TIMING_GET(num)}   
#else 
   #define END_TIMING(num)
#endif       

#ifdef TIMING       
   #define TIMING_GET1(num) disc.common.timing[num] += chrono::duration_cast<chrono::nanoseconds>(end-begin).count()/1e6;        
#else 
   #define TIMING_GET1(num)  
#endif                      

#ifdef TIMING
   #define END_TIMING_DISC(num) {CUDA_SYNC; TIMING_END; TIMING_GET1(num)}   
#else 
   #define END_TIMING_DISC(num)
#endif       
                
#ifdef HAVE_CUDA     

#ifdef USE_FLOAT
#define cublasNRM2 cublasSnorm2
#define cublasDOT cublasSdot
#define cublasAXPY cublasSaxpy
#define cublasGEMV cublasSgemv
#define cublasGEMM cublasSgemm
#define cublasGEMVBatched cublasSgemvBatched
#define cublasGEMMBatched cublasSgemmBatched
#define cublasGEMVStridedBatched cublasSgemvStridedBatched
#define cublasGEMMStridedBatched cublasSgemmStridedBatched
#define cublasGETRF cublasSgetrf
#define cublasGETRI cublasSgetri
#define cublasGETRFBatched cublasSgetrfBatched
#define cublasGETRIBatched cublasSgetriBatched
#define cublasTRSM cublasStrsm 
#else
#define cublasNRM2 cublasDnorm2
#define cublasDOT cublasDdot
#define cublasAXPY cublasDaxpy
#define cublasGEMV cublasDgemv
#define cublasGEMM cublasDgemm
#define cublasGEMVBatched cublasDgemvBatched
#define cublasGEMMBatched cublasDgemmBatched
#define cublasGEMVStridedBatched cublasDgemvStridedBatched
#define cublasGEMMStridedBatched cublasDgemmStridedBatched
#define cublasGETRF cublasDgetrf
#define cublasGETRI cublasDgetri
#define cublasGETRFBatched cublasDgetrfBatched
#define cublasGETRIBatched cublasDgetriBatched
#define cublasTRSM cublasDtrsm 
#endif

#define CHECK(call)                                                            \
{                                                                              \
    const cudaError_t error = call;                                            \
    if (error != cudaSuccess)                                                  \
    {                                                                          \
        fprintf(stderr, "Error: %s:%d, ", __FILE__, __LINE__);                 \
        fprintf(stderr, "code: %d, reason: %s\n", error,                       \
                cudaGetErrorString(error));                                    \
        exit(1);                                                               \
    }                                                                          \
}

#define CHECK_CUBLAS(call)                                                     \
{                                                                              \
    cublasStatus_t err;                                                        \
    if ((err = (call)) != CUBLAS_STATUS_SUCCESS)                               \
    {                                                                          \
        fprintf(stderr, "Got CUBLAS error %d at %s:%d\n", err, __FILE__,       \
                __LINE__);                                                     \
        exit(1);                                                               \
    }                                                                          \
}

#define GPUFREE(x)                                                       \
{                                                                         \
    if (x != nullptr) {                                                      \
        cudaTemplateFree(x);                                              \
        x = nullptr;                                                         \
    }                                                                     \
}

template <typename T> static void cudaTemplateMalloc(T **d_data, Int n)
{
    // allocate the memory on the GPU            
    CHECK( cudaMalloc( (void**)d_data, n * sizeof(T) ) );
}

template <typename T> static void cudaTemplateMallocManaged(T **d_data, Int n)
{
    // allocate unified memory 
    CHECK( cudaMallocManaged( (void**)d_data, n * sizeof(T) ) );        
}

template <typename T> static void cudaTemplateHostAlloc(T **h_data, Int n, unsigned int flags)
{
    // allocate zero-copy memory on host    
    CHECK(cudaHostAlloc((void **)h_data, n * sizeof(T), flags));                
}

template <typename T> static void cudaTemplateHostAllocMappedMemory(T **h_data, Int n)
{
    // allocate zero-copy memory on host    
    CHECK(cudaHostAlloc((void **)h_data, n * sizeof(T), cudaHostAllocMapped));                
}

template <typename T> static void cudaTemplateHostAllocPinnedMemory(T **h_data, Int n)
{
    // allocate pinned memory on host        
    CHECK(cudaHostAlloc((void **)h_data, n * sizeof(T), cudaHostAllocDefault));                
}

template <typename T> static void cudaTemplateFree(T *d_data)
{
    // free the memory on the GPU            
    CHECK( cudaFree( d_data ) );    
}

template <typename T> static void cudaCopytoDevice(T *d_data, T *h_data, Int n)
{
    // copy data from CPU to GPU
    CHECK( cudaMemcpy( d_data, h_data, n * sizeof(T), cudaMemcpyHostToDevice ) );    
}

template <typename T> static void cudaCopytoHost(T *h_data, T *d_data, Int n)
{
    // copy data from GPU to CPU
    CHECK( cudaMemcpy( h_data, d_data, n * sizeof(T), cudaMemcpyDeviceToHost ) );    
}

#endif

template <typename T> static void TemplateMalloc(T **data, Int n, Int backend)
{
    if (backend == 0)               
        *data = (T *) malloc(n*sizeof(T));      

#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        // allocate the memory on the GPU            
        CHECK( cudaMalloc( (void**)data, n * sizeof(T) ) );
#endif                 
}

template <typename T> static void TemplateFree(T *data,  Int backend)
{
    if (backend == 0)  CPUFREE(data);
        
#ifdef HAVE_CUDA            
    if (backend == 2)  GPUFREE(data);
#endif                  
}

template <typename T> static void TemplateCopytoDevice(T *d_data, T *h_data, Int n, Int backend)
{
    if (backend == 0)  {
        for (Int i=0; i<n; i++)
            d_data[i] = h_data[i];
    }
    
#ifdef HAVE_CUDA            
    // copy data from CPU to GPU
    if (backend == 2) CHECK( cudaMemcpy( d_data, h_data, n * sizeof(T), cudaMemcpyHostToDevice ) );            
#endif    
}

template <typename T> static void TemplateCopytoHost(T *h_data, T *d_data, Int n, Int backend)
{
    if (backend == 0)  {
        for (Int i=0; i<n; i++)
            h_data[i] = d_data[i];
    }

#ifdef HAVE_CUDA
    // copy data from GPU to CPU
    if (backend == 2) CHECK( cudaMemcpy( h_data, d_data, n * sizeof(T), cudaMemcpyDeviceToHost ) );
#endif    
}

struct appstruct {              
    Int *lsize=nullptr;
    Int *nsize=nullptr;  // data size
    Int *ndims=nullptr;  // dimensions
    Int *flag=nullptr;   // flag parameters
    Int *problem=nullptr;// problem parameters    
    Int *comm=nullptr;   // communication parameters 
    Int *porder=nullptr; // polymnomial degrees
    Int *stgib=nullptr;
    Int *vindx=nullptr;
    
    dstype *uinf=nullptr;    // boundary data
    dstype *dt=nullptr;      // time steps       
    dstype *dae_dt=nullptr;  // dual time steps       
    dstype *factor=nullptr;  // factors      
    dstype *physicsparam=nullptr; // physical parameters
    dstype *solversparam=nullptr; // solvers parameters
    dstype *tau=nullptr; // stabilization parameters
    dstype *stgdata=nullptr; 
    dstype *stgparam=nullptr; 
    
    //dstype time=nullptr;     /* current time */
    dstype *fc_u=nullptr;    /* factor when discretizing the time derivative of the U equation. Allow scalar field for local time stepping in steady problems? */
    dstype *fc_q=nullptr;    /* factor when discretizing the time derivative of the Q equation. Allow scalar field for local time stepping in steady problems? */
    dstype *fc_w=nullptr;    /* factor when discretizing the time derivative of the P equation. Allow scalar field for local time stepping in steady problems? */    
    
    dstype *dtcoef_u=nullptr;    /* factor when discretizing the time derivative of the U equation. Allow scalar field for local time stepping in steady problems? */
    dstype *dtcoef_q=nullptr;    /* factor when discretizing the time derivative of the Q equation. Allow scalar field for local time stepping in steady problems? */
    dstype *dtcoef_w=nullptr;    /* factor when discretizing the time derivative of the P equation. Allow scalar field for local time stepping in steady problems? */    
    
    Int szflag=0, szproblem=0, szcomm=0, szporder=0, szstgib=0, szvindx=0;
    Int szuinf=0, szdt=0, szdae_dt=0, szfactor=0, szphysicsparam=0, szsolversparam=0;
    Int sztau=0, szstgdata=0, szstgparam=0, szfc_u=0, szfc_q=0, szfc_w=0;
    Int szdtcoef_u=0, szdtcoef_q=0, szdtcoef_w=0;

    int sizeofint() {
      int sz = szflag + szproblem + szcomm + szporder + szstgib + szvindx;
      return sz;
    }
    int sizeoffloat() {
      int sz = szuinf+szdt+szdae_dt+szfactor+szphysicsparam+szsolversparam+
               sztau+szstgdata+szstgparam+szfc_u+szfc_q+szfc_w+szdtcoef_u+
               szdtcoef_q+szdtcoef_w;
      return sz;        
    }

    void printinfo()
    {    
      printf("--------------- App Struct Information ----------------\n");
      printf("size of flag: %d\n", szflag);
      printf("size of problem: %d\n", szproblem);
      printf("size of comm: %d\n", szcomm);
      printf("size of porder: %d\n", szporder);
      printf("size of stgib: %d\n", szstgib);
      printf("size of vindx: %d\n", szvindx);
      printf("size of uinf: %d\n", szuinf);
      printf("size of dt: %d\n", szdt);
      printf("size of dae_dt: %d\n", szdae_dt);
      printf("size of factor: %d\n", szfactor);
      printf("size of physicsparam: %d\n", szphysicsparam);
      printf("size of solversparam: %d\n", szsolversparam);
      printf("size of tau: %d\n", sztau);
      printf("size of stgdata: %d\n", szstgdata);
      printf("size of stgparam: %d\n", szstgparam);
      printf("size of fc_u: %d\n", szfc_u);
      printf("size of fc_q: %d\n", szfc_q);
      printf("size of fc_w: %d\n", szfc_w);
      printf("size of dtcoef_u: %d\n", szdtcoef_u);
      printf("size of dtcoef_q: %d\n", szdtcoef_q);
      printf("size of dtcoef_w: %d\n", szdtcoef_w);
      printf("size of int: %d\n", sizeofint());
      printf("size of float: %d\n", sizeoffloat());
    }

    #ifdef HAVE_MPP
    Mutation::Mixture *mix=nullptr;
    #endif
    // custom destructor
    void freememory(Int hostmemory)
    {
       if (hostmemory==1) {
            CPUFREE(lsize);
            CPUFREE(nsize);
            CPUFREE(ndims);   
            CPUFREE(comm);   
            CPUFREE(porder);               
            CPUFREE(flag);    
            CPUFREE(problem);
            CPUFREE(stgib);
            CPUFREE(vindx);
            CPUFREE(uinf);
            CPUFREE(dt);
            CPUFREE(dae_dt);
            CPUFREE(factor);
            CPUFREE(physicsparam);
            CPUFREE(solversparam);
            CPUFREE(tau);
            CPUFREE(stgdata);
            CPUFREE(stgparam);
            CPUFREE(fc_u);
            CPUFREE(fc_q);
            CPUFREE(fc_w);
            CPUFREE(dtcoef_u);
            CPUFREE(dtcoef_q);
            CPUFREE(dtcoef_w);
        }            
#ifdef HAVE_CUDA      
       else {
            GPUFREE(lsize);
            GPUFREE(nsize);
            GPUFREE(ndims);    
            GPUFREE(comm);   
            GPUFREE(porder);   
            GPUFREE(flag);    
            GPUFREE(problem);
            GPUFREE(stgib);
            GPUFREE(vindx);
            GPUFREE(uinf);
            GPUFREE(dt);
            GPUFREE(dae_dt);
            GPUFREE(factor);
            GPUFREE(physicsparam);
            GPUFREE(solversparam);
            GPUFREE(tau);
            GPUFREE(stgdata);
            GPUFREE(stgparam);
            GPUFREE(fc_u);
            GPUFREE(fc_q);
            GPUFREE(fc_w);
            GPUFREE(dtcoef_u);
            GPUFREE(dtcoef_q);
            GPUFREE(dtcoef_w);
       }
#endif       
    }
};

struct masterstruct {       
    
    Int *lsize=nullptr;
    Int *nsize=nullptr;  // data size
    Int *ndims=nullptr;  // dimensions
    
    dstype *shapegwdotshapeg=nullptr;
    dstype *shapfgwdotshapfg=nullptr;
    dstype *shapegt=nullptr; // element shape functions at Gauss points (transpose)
    dstype *shapegw=nullptr; // element shape functions at Gauss points multiplied by Gauss weights
    dstype *shapfgt=nullptr; // face shape functions at Gauss points (transpose)
    dstype *shapfgw=nullptr; // face shape functions at Gauss points multiplied by Gauss weights    
    dstype *shapent=nullptr; // element shape functions at nodes (transpose)
    dstype *shapen=nullptr;  // element shape functions at nodes        
    dstype *shapfnt=nullptr; // element shape functions at nodes (transpose)
    dstype *shapfn=nullptr;  // element shape functions at nodes        
    dstype *xpe=nullptr; // nodal points on master element
    dstype *gpe=nullptr; // gauss points on master element
    dstype *gwe=nullptr; // gauss weighs on master element
    dstype *xpf=nullptr; // nodal points on master face
    dstype *gpf=nullptr; // gauss points on master face
    dstype *gwf=nullptr; // gauss weighs on master face
    
    dstype *shap1dgt=nullptr; 
    dstype *shap1dgw=nullptr; 
    dstype *shap1dnt=nullptr; 
    dstype *shap1dnl=nullptr; 
    dstype *xp1d=nullptr; // node points on 1D element
    dstype *gp1d=nullptr; // gauss points on 1D element
    dstype *gw1d=nullptr; // gauss weights on 1D element    
    
    Int szshapegwdotshapeg=0, szshapfgwdotshapfg=0, szshapegt=0, szshapegw=0;
    Int szshapfgt=0, szshapfgw=0, szshapent=0, szshapen=0, szshapfnt=0;
    Int szshapfn=0, szxpe=0, szgpe=0, szgwe=0, szxpf=0, szgpf=0, szgwf=0;
    Int szshap1dgt=0, szshap1dgw=0, szshap1dnt=0, szshap1dnl=0, szxp1d=0;
    Int szgp1d=0, szgw1d=0;

    int sizeofint() { return 0;} 

    int sizeoffloat()
    {
      int sz = szshapegwdotshapeg+szshapfgt + szshapfgw + szshapent + szshapen +
               szshapfnt + szshapfn + szxpe + szgpe + szgwe + szxpf + szgpf + 
               szgwf + szshap1dgt + szshap1dgw + szshap1dnt + szshap1dnl + 
               szxp1d + szgp1d + szgw1d;
      return sz;         
    }

    void printinfo()
    {
      printf("--------------- Master Struct Information ----------------\n");
      printf("size of shapegwdotshapeg: %d\n", szshapegwdotshapeg);
      printf("size of shapfgwdotshapfg: %d\n", szshapfgwdotshapfg);
      printf("size of shapegt: %d\n", szshapegt);
      printf("size of shapegw: %d\n", szshapegw);
      printf("size of shapfgt: %d\n", szshapfgt);
      printf("size of shapfgw: %d\n", szshapfgw);
      printf("size of shapent: %d\n", szshapent);
      printf("size of shapen: %d\n", szshapen);
      printf("size of shapfnt: %d\n", szshapfnt);
      printf("size of shapfn: %d\n", szshapfn);
      printf("size of xpe: %d\n", szxpe);
      printf("size of gpe: %d\n", szgpe);
      printf("size of gwe: %d\n", szgwe);
      printf("size of xpf: %d\n", szxpf);
      printf("size of gpf: %d\n", szgpf);
      printf("size of gwf: %d\n", szgwf);
      printf("size of shap1dgt: %d\n", szshap1dgt);
      printf("size of shap1dgw: %d\n", szshap1dgw);
      printf("size of shap1dnt: %d\n", szshap1dnt);
      printf("size of shap1dnl: %d\n", szshap1dnl);
      printf("size of xp1d: %d\n", szxp1d);
      printf("size of gp1d: %d\n", szgp1d);
      printf("size of gw1d: %d\n", szgw1d);
      printf("size of int: %d\n", sizeofint());
      printf("size of float: %d\n", sizeoffloat());
    }

    void freememory(Int hostmemory)
    {
       if (hostmemory==1) {
            CPUFREE(lsize);
            CPUFREE(nsize);
            CPUFREE(ndims);    
            CPUFREE(shapegt); // element shape functions at Gauss points (transpose)
            CPUFREE(shapegw); // element shape functions at Gauss points multiplied by Gauss weights
            CPUFREE(shapfgt); // face shape functions at Gauss points (transpose)
            CPUFREE(shapfgw); // face shape functions at Gauss points multiplied by Gauss weights    
            CPUFREE(shapent); // element shape functions at nodes (transpose)
            CPUFREE(shapen);  // element shape functions at nodes        
            CPUFREE(shapfnt); // element shape functions at nodes (transpose)
            CPUFREE(shapfn);  // element shape functions at nodes        
            CPUFREE(xpe); // nodal points on master element
            CPUFREE(gpe); // gauss points on master element
            CPUFREE(gwe); // gauss weighs on master element
            CPUFREE(xpf); // nodal points on master face
            CPUFREE(gpf); // gauss points on master face
            CPUFREE(gwf); // gauss weighs on master face            
            CPUFREE(shap1dgt); 
            CPUFREE(shap1dgw); 
            CPUFREE(shap1dnt); 
            CPUFREE(shap1dnl); 
            CPUFREE(xp1d); 
            CPUFREE(gp1d); 
            CPUFREE(gw1d);            
            CPUFREE(shapegwdotshapeg);
            CPUFREE(shapfgwdotshapfg);             
        }            
#ifdef HAVE_CUDA      
       else {
            GPUFREE(lsize);
            GPUFREE(nsize);
            GPUFREE(ndims);    
            GPUFREE(shapegt); // element shape functions at Gauss points (transpose)
            GPUFREE(shapegw); // element shape functions at Gauss points multiplied by Gauss weights
            GPUFREE(shapfgt); // face shape functions at Gauss points (transpose)
            GPUFREE(shapfgw); // face shape functions at Gauss points multiplied by Gauss weights    
            GPUFREE(shapent); // element shape functions at nodes (transpose)
            GPUFREE(shapen);  // element shape functions at nodes        
            GPUFREE(shapfnt); // element shape functions at nodes (transpose)
            GPUFREE(shapfn);  // element shape functions at nodes        
            GPUFREE(xpe); // nodal points on master element
            GPUFREE(gpe); // gauss points on master element
            GPUFREE(gwe); // gauss weighs on master element
            GPUFREE(xpf); // nodal points on master face
            GPUFREE(gpf); // gauss points on master face
            GPUFREE(gwf); // gauss weighs on master face
            GPUFREE(shap1dgt); 
            GPUFREE(shap1dgw); 
            GPUFREE(shap1dnt); 
            GPUFREE(shap1dnl); 
            GPUFREE(xp1d); 
            GPUFREE(gp1d); 
            GPUFREE(gw1d);             
            GPUFREE(shapegwdotshapeg);
            GPUFREE(shapfgwdotshapfg);                        
       }
#endif       
    }    
};
  
struct meshstruct {       
    Int *lsize=nullptr;
    Int *nsize=nullptr;  // data size
    Int *ndims=nullptr;  // dimensions
        
    Int *facecon=nullptr;    // face-to-element connectivities 
    Int *f2e=nullptr;        // face-to-element connectivities
    Int *elemcon=nullptr;    // element-to-face connectivities
    Int *perm=nullptr;       // indices of element nodes on faces
    Int *bf=nullptr;         // boundary faces  
    Int *boufaces=nullptr;   // boundary faces
    //Int *interfacefaces=nullptr;   // interface faces between two subdomains
    Int *eblks=nullptr;    // element blocks
    Int *fblks=nullptr;    // face blocks    
    Int *nbsd=nullptr;
    Int *elemsend=nullptr;
    Int *elemrecv=nullptr;
    Int *elemsendpts=nullptr;
    Int *elemrecvpts=nullptr;
    Int *elempart=nullptr;
    Int *elempartpts=nullptr;
    Int *cgelcon=nullptr;
    Int *rowent2elem=nullptr;
    Int *cgent2dgent=nullptr;
    Int *colent2elem=nullptr;
    Int *rowe2f1=nullptr;
    Int *cole2f1=nullptr;
    Int *ent2ind1=nullptr;
    Int *rowe2f2=nullptr;
    Int *cole2f2=nullptr;
    Int *ent2ind2=nullptr;
    
    Int *findxdg1=nullptr; 
    //Int *findxdg2=nullptr; 
    Int *findxdgp=nullptr; 
    Int *findudg1=nullptr; 
    Int *findudg2=nullptr; 
    Int *findudgp=nullptr;     
    Int *eindudg1=nullptr;     
    Int *eindudgp=nullptr;     
    Int *elemsendind=nullptr;
    Int *elemrecvind=nullptr;
    Int *elemsendodg=nullptr;
    Int *elemrecvodg=nullptr;
    Int *elemsendudg=nullptr;
    Int *elemrecvudg=nullptr;
    //Int *index=nullptr; 
    
    Int szfacecon=0, szf2e=0, szelemcon=0, szperm=0, szbf=0, szboufaces=0; 
    Int szeblks=0, szfblks=0, sznbsd=0, szelemsend=0;
    Int szelemrecv=0, szelemsendpts=0, szelemrecvpts=0, szelempart=0;
    Int szelempartpts=0, szcgelcon=0, szrowent2elem=0, szcgent2dgent=0;
    Int szcolent2elem=0, szrowe2f1=0, szcole2f1=0, szent2ind1=0, szrowe2f2=0;
    Int szcole2f2=0, szent2ind2=0, szfindxdg1=0, szfindxdg2=0, szfindxdgp=0; 
    Int szfindudg1=0, szfindudg2=0, szfindudgp=0, szeindudg1=0, szeindudgp=0;
    Int szelemsendind=0, szelemrecvind=0, szelemsendodg=0, szelemrecvodg=0;
    Int szelemsendudg=0, szelemrecvudg=0, szindex=0;

    int sizeoffloat() {return 0;}
    int sizeofint() {
      int sz = szeblks+szfblks + sznbsd + szelemsend + szelemrecv + 
               szelemsendpts + szelemrecvpts + szelempart + szelempartpts + 
               szcgelcon + szrowent2elem + szcgent2dgent + szcolent2elem + 
               szrowe2f1 + szcole2f1 + szent2ind1 + szrowe2f2 + szcole2f2 + 
               szent2ind2 + szfindxdg1 + szfindxdgp + szfindudg1 + szfindudg2 + 
               szfindudgp + szeindudg1 + szeindudgp + szelemsendind + szelemrecvind + 
               szelemsendodg + szelemrecvodg + szelemsendudg + szelemrecvudg;
      return sz;        
    }

    void printinfo()
    {
      printf("--------------- Mesh Struct Information ----------------\n");
      printf("size of facecon: %d\n", szfacecon);
      printf("size of f2e: %d\n", szf2e);
      printf("size of elemcon: %d\n", szelemcon);
      printf("size of perm: %d\n", szperm);
      printf("size of bf: %d\n", szbf);
      printf("size of boufaces: %d\n", szboufaces);
      //printf("size of interfacefaces: %d\n", szinterfacefaces);
      printf("size of eblks: %d\n", szeblks);
      printf("size of fblks: %d\n", szfblks);
      printf("size of nbsd: %d\n", sznbsd);
      printf("size of elemsend: %d\n", szelemsend);
      printf("size of elemrecv: %d\n", szelemrecv);
      printf("size of elemsendpts: %d\n", szelemsendpts);
      printf("size of elemrecvpts: %d\n", szelemrecvpts);
      printf("size of elempart: %d\n", szelempart);
      printf("size of elempartpts: %d\n", szelempartpts);
      printf("size of cgelcon: %d\n", szcgelcon);
      printf("size of rowent2elem: %d\n", szrowent2elem);
      printf("size of cgent2dgent: %d\n", szcgent2dgent);
      printf("size of colent2elem: %d\n", szcolent2elem);
      printf("size of rowe2f1: %d\n", szrowe2f1);
      printf("size of cole2f1: %d\n", szcole2f1);
      printf("size of ent2ind1: %d\n", szent2ind1);
      printf("size of rowe2f2: %d\n", szrowe2f2);
      printf("size of cole2f2: %d\n", szcole2f2);
      printf("size of ent2ind2: %d\n", szent2ind2);
      printf("size of findxdg1: %d\n", szfindxdg1);
      //printf("size of findxdg2: %d\n", szfindxdg2);
      printf("size of findxdgp: %d\n", szfindxdgp);
      printf("size of findudg1: %d\n", szfindudg1);
      printf("size of findudg2: %d\n", szfindudg2);
      printf("size of findudgp: %d\n", szfindudgp);
      printf("size of eindudg1: %d\n", szeindudg1);
      printf("size of eindudgp: %d\n", szeindudgp);
      printf("size of elemsendind: %d\n", szelemsendind);
      printf("size of elemrecvind: %d\n", szelemrecvind);
      printf("size of elemsendodg: %d\n", szelemsendodg);
      printf("size of elemrecvodg: %d\n", szelemrecvodg);
      printf("size of elemsendudg: %d\n", szelemsendudg);
      printf("size of elemrecvudg: %d\n", szelemrecvudg);      
      printf("size of int: %d\n", sizeofint());
      printf("size of float: %d\n", sizeoffloat());
    }

    void freememory(Int hostmemory)
    {
       if (hostmemory==1) {
            CPUFREE(lsize);
            CPUFREE(nsize);
            CPUFREE(ndims);    
            CPUFREE(facecon);    // face-to-element connectivities 
            CPUFREE(f2e);    // face-to-element connectivities 
            CPUFREE(elemcon);    // element-to-face connectivities
            CPUFREE(perm);       // indices of element nodes on faces
            CPUFREE(boufaces);   // boundary faces
            //CPUFREE(interfacefaces);   // interface faces
            CPUFREE(eblks);    // element blocks
            CPUFREE(fblks);    // face blocks    
            CPUFREE(nbsd);
            CPUFREE(elemsend);
            CPUFREE(elemrecv);
            CPUFREE(elemsendpts);
            CPUFREE(elemrecvpts);            
            CPUFREE(elempart);
            CPUFREE(elempartpts);   
            CPUFREE(cgelcon);            
            CPUFREE(rowent2elem);
            CPUFREE(cgent2dgent);
            CPUFREE(colent2elem);
            CPUFREE(rowe2f1);
            CPUFREE(cole2f1);
            CPUFREE(ent2ind1);
            CPUFREE(rowe2f2);
            CPUFREE(cole2f2);
            CPUFREE(ent2ind2);
            
            CPUFREE(findxdg1);   
            //CPUFREE(findxdg2);   
            CPUFREE(findxdgp);   
            CPUFREE(findudg1);   
            CPUFREE(findudg2);   
            CPUFREE(findudgp);               
            CPUFREE(eindudg1);               
            CPUFREE(eindudgp);  
            CPUFREE(elemsendind);   
            CPUFREE(elemrecvind); 
            CPUFREE(elemsendodg);
            CPUFREE(elemrecvodg);
            CPUFREE(elemsendudg);
            CPUFREE(elemrecvudg);
            //CPUFREE(index); 
        }            
#ifdef HAVE_CUDA      
       else {
            GPUFREE(lsize);
            GPUFREE(nsize);
            GPUFREE(ndims);
            GPUFREE(facecon);    // face-to-element connectivities 
            GPUFREE(f2e);    // face-to-element connectivities 
            GPUFREE(elemcon);    // element-to-face connectivities
            GPUFREE(perm);       // indices of element nodes on faces
            GPUFREE(boufaces);   // boundary faces
            //GPUFREE(interfacefaces);   // interface faces
            GPUFREE(eblks);    // element blocks
            GPUFREE(fblks);    // face blocks    
            GPUFREE(nbsd);
            GPUFREE(elemsend);
            GPUFREE(elemrecv);
            GPUFREE(elemsendpts);
            GPUFREE(elemrecvpts);     
            GPUFREE(elempart);
            GPUFREE(elempartpts);   
            GPUFREE(cgelcon);            
            GPUFREE(rowent2elem);
            GPUFREE(cgent2dgent);
            GPUFREE(colent2elem);
            GPUFREE(rowe2f1);
            GPUFREE(cole2f1);
            GPUFREE(ent2ind1);
            GPUFREE(rowe2f2);
            GPUFREE(cole2f2);
            GPUFREE(ent2ind2);
            
            GPUFREE(findxdg1);   
            //GPUFREE(findxdg2);   
            GPUFREE(findxdgp);   
            GPUFREE(findudg1);   
            GPUFREE(findudg2);   
            GPUFREE(findudgp);               
            GPUFREE(eindudg1);               
            GPUFREE(eindudgp);   
            GPUFREE(elemsendind);   
            GPUFREE(elemrecvind); 
            GPUFREE(elemsendodg);
            GPUFREE(elemrecvodg); 
            GPUFREE(elemsendudg);
            GPUFREE(elemrecvudg);   
            //GPUFREE(index);   
       }
#endif       
    }        
};

struct solstruct {       
    Int *lsize=nullptr;
    Int *nsize=nullptr;  // data size
    Int *ndims=nullptr;  // dimensions
    
    dstype *xdg=nullptr; // spatial coordinates
    dstype *udg=nullptr; // solution (u, q, p) 
    dstype *sdg=nullptr; // source term due to the previous solution
    dstype *odg=nullptr; // auxilary term 
    dstype *wdg=nullptr; // dw/dt = u (wave problem)
    dstype *uh=nullptr;  // uhat  
    // dstype *dudgt=nullptr; // for line search
    dstype *udg0=nullptr;
    dstype *uh0=nullptr;

    #ifdef HAVE_ENZYME
        dstype *dudg=nullptr; // solution (du, dq, dp) 
        dstype *dwdg=nullptr; // dw/dt = u (wave problem)
        dstype *duh=nullptr; // duhat
        dstype *dodg=nullptr;
        dstype *dodgg=nullptr;
        dstype *dog1=nullptr;
        dstype *dog2=nullptr;
    #endif
    dstype *elemg=nullptr;
    dstype *faceg=nullptr;
    dstype *elemfaceg=nullptr;        
    //dstype *udgg=nullptr;
    dstype *sdgg=nullptr;
    dstype *odgg=nullptr;
    dstype *og1=nullptr;
    dstype *og2=nullptr;    
    dstype *udgavg=nullptr; // time-average solution (u, q, p) 
    dstype *wsrc=nullptr;   // source term due to the time derivative for DAE equations  
    dstype *wdual=nullptr;   // source term due to the dual time derivative for DAE equations  
    dstype** udgarray;
    
    Int szxdg=0, szudg=0, szsdg=0, szodg=0, szwdg=0, szuh=0;
    Int szelemg=0, szfaceg=0, szelemfaceg=0, szsdgg=0, szodgg=0, szog1=0, szog2=0;
    Int szudgavg=0, szwsrc=0, szwdual=0;

    int sizeofint() {return 0;}
    int sizeoffloat() {
      int sz = szxdg + szudg + szsdg + szodg + szwdg + szuh + szelemg + szfaceg +
               szelemfaceg + szsdgg + szodgg + szog1 + szog2 + szudgavg + 
               szwsrc + szwdual;
      return sz;
    }

    void printinfo()
    {
      printf("--------------- Solution Struct Information ----------------\n");
      printf("size of xdg: %d\n", szxdg);
      printf("size of udg: %d\n", szudg);
      printf("size of sdg: %d\n", szsdg);
      printf("size of odg: %d\n", szodg);
      printf("size of wdg: %d\n", szwdg);
      printf("size of uh: %d\n", szuh);
      printf("size of elemg: %d\n", szelemg);
      printf("size of faceg: %d\n", szfaceg);
      printf("size of elemfaceg: %d\n", szelemfaceg);
      printf("size of sdgg: %d\n", szsdgg);
      printf("size of odgg: %d\n", szodgg);
      printf("size of og1: %d\n", szog1);
      printf("size of og2: %d\n", szog2);
      printf("size of udgavg: %d\n", szudgavg);
      printf("size of wsrc: %d\n", szwsrc);
      printf("size of wdual: %d\n", szwdual);     
      printf("size of int: %d\n", sizeofint());
      printf("size of float: %d\n", sizeoffloat());   
    } 

    void freememory(Int hostmemory)
    {
       if (hostmemory==1) {
            CPUFREE(lsize);
            CPUFREE(nsize);
            CPUFREE(ndims);    
            CPUFREE(xdg); // spatial coordinates
            CPUFREE(udg); // solution (u, q, p) 
            CPUFREE(sdg); // source term due to the previous solution
            CPUFREE(odg); // auxilary term 
            CPUFREE(wdg); // wave problem
#ifdef HAVE_ENZYME                   
            CPUFREE(dudg); // solution (u, q, p) 
            CPUFREE(dwdg); // wave problem
            CPUFREE(duh);
            CPUFREE(dodg);
            CPUFREE(dodgg);
            CPUFREE(dog1);
            CPUFREE(dog2);
#endif            
            CPUFREE(uh);  // uhat      
            CPUFREE(elemg); 
            CPUFREE(faceg); 
            CPUFREE(elemfaceg);
            //CPUFREE(udgg); 
            CPUFREE(sdgg); 
            CPUFREE(odgg); 
            CPUFREE(wsrc);
            CPUFREE(wdual);
            CPUFREE(og1);
            CPUFREE(og2);     
            CPUFREE(udgavg); // time-average solution (u, q, p) 
            //delete[] udgarray;
        }            
#ifdef HAVE_CUDA      
       else {
            GPUFREE(lsize);
            GPUFREE(nsize);
            GPUFREE(ndims);
            GPUFREE(xdg); // spatial coordinates
            GPUFREE(udg); // solution (u, q, p) 
            GPUFREE(sdg); // source term due to the previous solution
            GPUFREE(odg); // auxilary term 
            GPUFREE(wdg); // wave problem
#ifdef HAVE_ENZYME                   
            GPUFREE(dudg); // solution (u, q, p) 
            GPUFREE(dwdg); // wave problem
            GPUFREE(duh);
            GPUFREE(dodg);
            GPUFREE(dodgg);
            GPUFREE(dog1);
            GPUFREE(dog2);
#endif                        
            GPUFREE(uh);  // uhat          
            GPUFREE(elemg); 
            GPUFREE(faceg); 
            GPUFREE(elemfaceg);
            //GPUFREE(udgg); 
            GPUFREE(sdgg); 
            GPUFREE(odgg); 
            GPUFREE(wsrc); 
            GPUFREE(wdual); 
            GPUFREE(og1);
            GPUFREE(og2);        
            GPUFREE(udgavg); // time-average solution (u, q, p) 
            //delete[] udgarray;
       }
#endif       
    }            
};

struct resstruct {
    //dstype *R=nullptr;    // shared memory for all residual vectors
    dstype *Rqe=nullptr;  // element residual vector for q
    dstype *Rqf=nullptr;  // face residual vector for q   
    dstype *Rue=nullptr;  // element residual vector for u
    dstype *Ruf=nullptr;  // face residual vector for u
    dstype *Rq=nullptr;   // residual vector for q     
    dstype *Ru=nullptr;   // residual vector for u    
    dstype *Rp=nullptr;   // residual vector for p    
    dstype *Rh=nullptr;   // residual vector for uhat
    dstype *dudgt=nullptr;   // residual vector for uhat
    dstype *Rh=nullptr;   // residual vector for uhat    

    dstype *dRq=nullptr;   // residual vector for q     
    dstype *dRu=nullptr;   // residual vector for u        
    dstype *dRh=nullptr;   // residual vector for uhat
    dstype *dRqe=nullptr;  // element residual vector for q
    dstype *dRqf=nullptr;  // face residual vector for q   
    dstype *dRue=nullptr;  // element residual vector for u
    dstype *dRuf=nullptr;  // face residual vector for u    

    dstype *Mass=nullptr; // store the mass matrix
    dstype *Minv=nullptr; // store the inverse of the mass matrix
    dstype *Mass2=nullptr; // store the mass matrix
    dstype *Minv2=nullptr; // store the inverse of the mass matrix
    dstype *C=nullptr; // store the convection matrix
    dstype *E=nullptr; // store the convection matrix
    dstype *D=nullptr; // store the diffusion matrix
    dstype *B=nullptr; // store the diffusion matrix
    dstype *F=nullptr; // store the diffusion matrix
    dstype *G=nullptr; // store the diffusion matrix
    dstype *K=nullptr; // store the diffusion matrix
    dstype *H=nullptr; // store the diffusion matrix

    Int *ipiv=nullptr;    
        
    Int szipiv=0, szH=0, szK=0, szG=0, szF=0, szB=0, szD=0, szE=0, szC=0, szMass=0, szMinv=0, szMass2=0, szMinv2=0;
    Int szRq=0, szRu=0, szRh=0, szRuf=0, szRue=0, szRqf=0, szRqe=0;  

    int sizeofint() {return szipiv;}
    int sizeoffloat() {
      int sz = szH + szK + szG + szF + szB + szD + szE + szC + szMass + szMinv +
               szMass2 + szMinv2 + szRq + szRu + szRh + szRuf + szRue + szRqf + 
               szRqe;        
      return sz;
    }

    void printinfo()
    {
      printf("--------------- Residual Struct Information ----------------\n");
      printf("size of ipiv: %d\n", szipiv);
      printf("size of Rq: %d\n", szRq);
      printf("size of Ru: %d\n", szRu);
      printf("size of Rh: %d\n", szRh);
      // printf("size of dRq: %d\n", szdRq);
      // printf("size of dRu: %d\n", szdRu);
      // printf("size of dRh: %d\n", szdRh);
      printf("size of Mass: %d\n", szMass);
      printf("size of Minv: %d\n", szMinv);
      printf("size of Mass2: %d\n", szMass2);
      printf("size of Minv2: %d\n", szMinv2);
      printf("size of C: %d\n", szC);
      printf("size of E: %d\n", szE);
      printf("size of D: %d\n", szD);
      printf("size of B: %d\n", szB);
      printf("size of F: %d\n", szF);
      printf("size of G: %d\n", szG);
      printf("size of K: %d\n", szK);
      printf("size of H: %d\n", szH);
      printf("size of int: %d\n", sizeofint());
      printf("size of float: %d\n", sizeoffloat());
    }

    void freememory(Int hostmemory)
    {
       if (hostmemory==1) {            
            CPUFREE(Rq);    
            CPUFREE(Ru);    
            CPUFREE(Rh);    
            CPUFREE(dudgt);
#ifdef HAVE_ENZYME                   
            CPUFREE(dRq);   
            CPUFREE(dRu);   
            CPUFREE(dRh);   
#endif                                                
            CPUFREE(Mass);      
            CPUFREE(Minv);      
            CPUFREE(Mass2);      
            CPUFREE(Minv2);      
            CPUFREE(C);      
            CPUFREE(E);      
            CPUFREE(D);
            CPUFREE(B);
            CPUFREE(F);
            CPUFREE(G);
            CPUFREE(K);
            CPUFREE(H);
            CPUFREE(ipiv);
        }            
#ifdef HAVE_CUDA      
       else {
            GPUFREE(Rq);   
            GPUFREE(Ru);   
            GPUFREE(Rh); 
            GPUFREE(dudgt); 
#ifdef HAVE_ENZYME                   
            GPUFREE(dRq);   
            GPUFREE(dRu);   
            GPUFREE(dRh);   
#endif                                    
            GPUFREE(Mass);   
            GPUFREE(Minv);   
            GPUFREE(Mass2);   
            GPUFREE(Minv2);   
            GPUFREE(C);   
            GPUFREE(E);   
            GPUFREE(D);
            GPUFREE(B);
            GPUFREE(F);
            GPUFREE(G);
            GPUFREE(K);
            GPUFREE(H);
            GPUFREE(ipiv);
       }
//        else {
//             cudaFreeHost(Rq);  
//             cudaFreeHost(Ru);  
//             cudaFreeHost(Rh);  
//             cudaFreeHost(Mass);            
//             cudaFreeHost(Minv);            
//        }
#endif       
    }                
};

struct tempstruct {
    dstype *tempn=nullptr; 
    dstype *tempg=nullptr;
    dstype *buffrecv=nullptr;
    dstype *buffsend=nullptr;
    
    int sztempn=0, sztempg = 0, szbuffrecv=0, szbuffsend=0;

    int sizeofint() {return 0;}
    int sizeoffloat() 
    {
      int sz = sztempn + sztempg + szbuffrecv + szbuffsend;
      return sz;
    }

    void printinfo()
    {
      printf("--------------- Temp Struct Information ----------------\n");
      printf("size of tempn: %d\n", sztempn);
      printf("size of tempg: %d\n", sztempg);
      printf("size of buffrecv: %d\n", szbuffrecv);
      printf("size of buffsend: %d\n", szbuffsend);
      printf("size of int: %d\n", sizeofint());
      printf("size of float: %d\n", sizeoffloat());
    }

    void freememory(Int hostmemory)
    {
       if (hostmemory==1) {
            CPUFREE(tempn); 
            CPUFREE(tempg); 
            CPUFREE(buffrecv); 
            CPUFREE(buffsend); 
        }            
#ifdef HAVE_CUDA      
       else {
            GPUFREE(tempn);
            GPUFREE(tempg);
            GPUFREE(buffrecv); 
            GPUFREE(buffsend); 
       }
#endif       
    }                     
};

struct sysstruct {    
    Int cpuMemory;
    Int *ipiv=nullptr;

    dstype *x=nullptr; 
    dstype *u=nullptr;
    dstype *r=nullptr;
    dstype *b=nullptr;
    dstype *v=nullptr;
    dstype *randvect=nullptr;
    dstype *q=nullptr;
    dstype *p=nullptr;
    
    // unified memory for GMRES solver
    dstype *tempmem=nullptr;
    dstype *lam=nullptr;
    //dstype *normcu=nullptr;    
    
    // store PTC matrix
    dstype *PTCmatrix=nullptr;    

    // for DIRK schemes
    dstype *utmp=nullptr;
    dstype *wtmp=nullptr;
    
    // previous solutions for time-dependent problems
    dstype *udgprev=nullptr; // 
    dstype *udgprev1=nullptr; // 
    dstype *udgprev2=nullptr; // 
    dstype *udgprev3=nullptr; //         
    dstype *wprev=nullptr; // 
    dstype *wprev1=nullptr; // 
    dstype *wprev2=nullptr; // 
    dstype *wprev3=nullptr; // 
        
    Int szipiv = 0;    
    Int szx=0, szu=0, szr=0, szb=0, szv=0, szq=0, szp=0;
    Int szrandvect=0, sztempmem=0, szlam=0, szPTCmatrix=0, szutmp=0, szwtmp=0;    
    Int szudgprev = 0, szudgprev1 = 0, szudgprev2 = 0, szudgprev3 = 0;
    Int szwprev = 0, szwprev1 = 0, szwprev2 = 0, szwprev3 = 0;
    
    int sizeofint() { return szipiv; }
    int sizeoffloat() {
      int sz = szx + szu + szr + szb + szv + szq + szp +
              szrandvect + sztempmem + szPTCmatrix + szutmp + szwtmp + 
              szudgprev + szudgprev1 + szudgprev2 + szudgprev3 + szwprev + 
              szwprev1 + szwprev2 + szwprev3;
      return sz;
    }

    void printinfo()
    {
      printf("--------------- Sys Struct Information ----------------\n");
      printf("size of ipiv: %d\n", szipiv);
      printf("size of x: %d\n", szx);
      printf("size of u: %d\n", szu);
      printf("size of r: %d\n", szr);
      printf("size of b: %d\n", szb);
      printf("size of v: %d\n", szv);
      printf("size of q: %d\n", szq);
      printf("size of p: %d\n", szp);
      printf("size of randvect: %d\n", szrandvect);
      printf("size of tempmem: %d\n", sztempmem);
      printf("size of lam: %d\n", szlam);
      printf("size of PTCmatrix: %d\n", szPTCmatrix);
      printf("size of utmp: %d\n", szutmp);
      printf("size of wtmp: %d\n", szwtmp);
      printf("size of udgprev: %d\n", szudgprev);
      printf("size of udgprev1: %d\n", szudgprev1);
      printf("size of udgprev2: %d\n", szudgprev2);
      printf("size of udgprev3: %d\n", szudgprev3);
      printf("size of wprev: %d\n", szwprev);
      printf("size of wprev1: %d\n", szwprev1);
      printf("size of wprev2: %d\n", szwprev2);
      printf("size of wprev3: %d\n", szwprev3);
      printf("size of int: %d\n", sizeofint());
      printf("size of float: %d\n", sizeoffloat());
    }

    void freememory(Int hostmemory)
    {
       CPUFREE(lam);  
       //CPUFREE(normcu);  
       CPUFREE(ipiv);  
       if (hostmemory==1) {
            CPUFREE(x); 
            CPUFREE(u); 
            CPUFREE(r); 
            CPUFREE(b); 
            CPUFREE(v); 
            CPUFREE(q); 
            CPUFREE(p); 
            CPUFREE(randvect);
            CPUFREE(tempmem);     
            CPUFREE(utmp);            
            //CPUFREE(w);
            //CPUFREE(wsrc);             
            CPUFREE(wtmp);             
            CPUFREE(udgprev);  
            CPUFREE(udgprev1);  
            CPUFREE(udgprev2);  
            CPUFREE(udgprev3);  
            CPUFREE(wprev);  
            CPUFREE(wprev1);  
            CPUFREE(wprev2);  
            CPUFREE(wprev3);  
            CPUFREE(PTCmatrix);  
        }            
#ifdef HAVE_CUDA      
       else {
            GPUFREE(x);
            GPUFREE(u);
            GPUFREE(r);
            GPUFREE(b);
            GPUFREE(v);
            GPUFREE(q);
            GPUFREE(p);
            GPUFREE(randvect);
            cudaFreeHost(tempmem);      
            //GPUFREE(w);
            //GPUFREE(wsrc); 
            GPUFREE(utmp);
            GPUFREE(wtmp);             
            GPUFREE(udgprev);  
            GPUFREE(udgprev1);  
            GPUFREE(udgprev2);  
            GPUFREE(udgprev3);  
            GPUFREE(wprev);  
            GPUFREE(wprev1);  
            GPUFREE(wprev2);  
            GPUFREE(wprev3);  
            GPUFREE(PTCmatrix);  
       }
#endif       
    }                     
};

struct precondstruct {    
    Int cpuMemory;
    
    dstype *W=nullptr; 
    dstype *U=nullptr;
    // dstype *V=nullptr;
    // dstype *R=nullptr;
    // dstype *C=nullptr;
    // dstype *H=nullptr;
    
    // dstype *y;
    // dstype *z=nullptr;
    
    // dstype *Cmat;
    
    // dstype *RBcoef=nullptr;     
    Int *ipiv=nullptr;
    
    Int szipiv = 0, szW = 0, szU = 0;
    int sizeofint() { return szipiv; }
    int sizeoffloat() {
      int sz = szW + szU;
      return sz;
    }

    void printinfo()
    { 
      printf("--------------- Precond Struct Information ----------------\n");
      printf("size of ipiv: %d\n", szipiv);
      printf("size of W: %d\n", szW);
      printf("size of U: %d\n", szU);
      printf("size of int: %d\n", sizeofint());
      printf("size of float: %d\n", sizeoffloat());
    }

    void freememory(Int hostmemory)
    {
       if (hostmemory==1) {
            CPUFREE(W); 
            CPUFREE(U); 
            // CPUFREE(V); 
            // CPUFREE(R); 
            // CPUFREE(C); 
            //CPUFREE(H); 
            //CPUFREE(y); 
            //CPUFREE(z); 
            CPUFREE(ipiv); 
        }            
#ifdef HAVE_CUDA      
       else {
            GPUFREE(W);
            GPUFREE(U);
            // GPUFREE(V);
            // GPUFREE(R);
            // GPUFREE(C);
            //GPUFREE(H);
            //GPUFREE(y); 
            //GPUFREE(z); 
            GPUFREE(ipiv);
       }
#endif       
    }                     
};

struct commonstruct {     
    string filein;       // Name of binary file with input data
    string fileout;      // Name of binary file to write the solution            
    
    Int backend;   // 0: Serial; 1: OpenMP; 2: CUDA  
    Int cpuMemory; // 1: data in CPU memory; 0: data in GPU memory
    Int maxnbc;    // maximum number of boundary conditions
    
    Int mpiRank;  // MPI rank      
    Int mpiProcs;    // number of MPI ranks
    Int nomodels; // number of models
    Int enzyme=0;
    
    Int nc; // number of compoments of (u, q)
    Int ncu;// number of compoments of (u)
    Int ncq;// number of compoments of (q)
    Int ncw;// number of compoments of (w)
    Int nco;// number of compoments of (o)
    Int nch;// number of compoments of (uhat)
    Int ncx;// number of compoments of (xdg)
    Int ncs;// number of compoments of (sdg)
    Int nce;// number of compoments of (edg)
    
    Int nd; // spatial dimension    
    Int elemtype;
    Int nodetype;
    Int porder; // polynomial degree for solution approximation
    Int pgauss; // polynomial degree for Gauss quadrature
    Int npe;  // number of nodes on master element
    Int npf;  // number of nodes on master face       
    Int nge;  // number of gauss poInts on master element
    Int ngf;  // number of gauss poInts on master face                    
    Int np1d;
    Int ng1d;    
    Int ppdegree=10; // polynomial preconditioner degree
    
    Int ne; // number of elements
    Int nf; // number of faces
    Int nv; // number of vertices      
    Int nfe; // number of faces per element        
    Int nbe; // number of blocks for elements 
    Int neb; // maximum number of elements per block
    Int nbf; // number of blocks for faces   
    Int nfb; // maximum number of faces per block
    Int nbe0; // number of blocks for interior elements 
    Int nbe1; // number of blocks for interior+interface elements 
    Int nbe2; // number of blocks for interior+interface+exterior elements 
    Int nbf0; // number of blocks for interior faces
    Int nbf1; // number of blocks for interior+interface faces
    
    Int ndof; // number of degrees of freedom of u
    Int ndofq; // number of degrees of freedom of q
    Int ndofw; // number of degrees of freedom of w
    Int ndofuhat; // number of degrees of freedom of uhat
    Int ndofudg; // number of degrees of freedom of udg
    Int ndofsdg; // number of degrees of freedom of sdg
    Int ndofodg; // number of degrees of freedom of odg
    Int ndofedg; // number of degrees of freedom of edg
    Int ntau; // length of the stabilization
    
    Int ne0;  // number of interior elements
    Int ne1;  // number of interior+interface elements
    Int ne2;  // number of interior+interface+exterior elements
    Int nf0;  // number of interior faces
    Int ndof1; // number of degrees of freedom of u
    Int ndofq1; // number of degrees of freedom of q
    Int ndofw1; // number of degrees of freedom of w
    Int ndofudg1; // number of degrees of freedom of udg
    Int ndofsdg1; // number of degrees of freedom of sdg
    Int ndofodg1; // number of degrees of freedom of odg
    Int ndofedg1; // number of degrees of freedom of edg
    Int ndofucg;
        
    Int RBdim; // maximum dimension of the reduced basis space
    Int RBcurrentdim; // current dimension of the reduced basis space
    Int RBremovedind; // the vector to be removed from the RB space and replaced with new vector    
    Int Wdim; // maxium dimension of W
    Int Wcurrentdim;
    
    Int extUhat=0; // external uhat function flag
    Int extFhat=0; // external fhat function flag
    Int extStab=0; // external stabilization function flag
    Int curvedMesh;// curved mesh    
    Int debugMode; // 1: save data to binary files for debugging
    Int appname;   /* 0: Euler; 1: Compressible Navier-Stokes; etc. */
    Int tdep;      // 0: steady-state; 1: time-dependent;  
    Int wave;      // wave problem    
    Int linearProblem; // 0: nonlinear problem;  1: linear problem
    Int subproblem=0;
    Int saveSolFreq;   // number of time steps to save the solution
    Int saveSolOpt;    // option to save the solution
    Int timestepOffset=0; // timestep offset to restart the simulation 
    Int stgNmode=0;       // number of synthetic turbulence generation modes
    Int tdfunc;           // time-derivative function flag
    Int source;           // source function flag
    Int modelnumber;      // model number
    Int ibs;              // boundary index to save solution 
    Int saveSolBouFreq=0; // number of time steps to save the solution on the boundary
    Int compudgavg=1;     // compute time-averaged solution udg
    Int readudgavg=0;     // flag to read time-averaged solution udg from file
    Int saveResNorm=0;   
    
    Int spatialScheme;   /* 0: HDG; 1: EDG; 2: IEDG, HEDG */
    Int temporalScheme;  // 0: DIRK; 1: BDF; 2: ERK
    Int torder;    /* temporal accuracy order */
    Int tstages;    /* DIRK stages */
    Int tsteps;    /* number of time steps */
    Int dae_steps=0; /* number of dual time steps */
    Int currentstep;  // current time step
    Int currentstage; // curent time stage
    Int convStabMethod;  // Flag for convective stabilization tensor. 0: Constant tau, 1: Lax-Friedrichs; 2: Roe.
    Int diffStabMethod;  // Flag for diffusive stabilization tensor. 0: No diffusive stabilization.
    Int rotatingFrame;   // Flag for rotating frame. Options: 0: Velocities are respect to a non-rotating frame. 1: Velocities are respect to a rotating frame.
    Int viscosityModel;  // Flag for viscosity model. 0: Constant kinematic viscosity; 1: Sutherland's law. 2: Hack for viscosity for compressible HIT case in (Johnsen, 2010)
    Int SGSmodel;        // Flag for sub-grid scale (SGS) model. Only available for 3D solver.
                                        //  0: No SGS model. 1: Static Smagorinsky/Yoshizawa/Knight model. 
                                        //  2: Static WALE/Yoshizawa/Knight model. 3: Static Vreman/Yoshizawa/Knight model.
                                        //  4: Dynamic Smagorinsky/Yoshizawa/Knight model.        
    Int ALEflag;                    // Flag for Arbitrary Lagrangian-Eulerian (ALE) formulation. 0: No ALE; 1: Translation; 2: Translation + rotation; 3: Translation + rotation + deformation
    Int ncAV;                     // Number of components for artificial viscosity; formulation specified in symbolic code as pde.avfield function   
    Int AVsmoothingIter;           //Number of times artificial viscosity is smoothed
    Int frozenAVflag;    // Flag deciding if artificial viscosity is calculated once per non-linear solve or in every residual evluation
                                //   0: AV not frozen, evaluated as part of residual
                                //   1: AV frozen, evluated once per solve (default)   

    Int linearSolver;  /* 0: GMRES; 1: CG; etc. */      
    Int nonlinearSolver;
    Int linearSolverMaxIter;   
    Int linearSolverIter;   // current iteration
    Int nonlinearSolverMaxIter;                
    Int nonlinearSolverIter;                
    Int matvecOrder;    
    Int gmresRestart;    
    Int gmresOrthogMethod;        
    Int preconditioner; // 0: low-rank; 1: reduced basis
    Int precMatrixType; // 0: identity; 1: inverse mass matrix
    Int ptcMatrixType;  // 0: identity; 1: mass matrix
    Int runmode;
    Int ninterfacefaces; // number of interface faces
    
    dstype dtfactor; // factor asssociated with udg for temporal discretization
    dstype time; // current simulation time    
    dstype matvecTol;
    dstype linearSolverTol;
    dstype linearSolverTolFactor=1.0;
    dstype nonlinearSolverTol;
    dstype linearSolverRelError;
    dstype rampFactor;               // Ramp factor for artificial viscosity flux
    dstype PTCparam;
    dstype tau0=0.0;
    dstype dae_alpha=1.0;
    dstype dae_beta=0.0;
    dstype dae_gamma=0.0;
    dstype dae_epsilon=0.0;
            
    Int* eblks=nullptr; // element blocks
    Int* fblks=nullptr; // face blocks   
    Int* ncarray=nullptr;
    Int* nboufaces=nullptr;
    
    Int nstgib;
    Int nnbsd; // number of neighboring subdomains
    Int nelemsend;
    Int nelemrecv;
    Int nvindx;
    Int* nbsd=nullptr; // neighboring subdomains
    Int* elemsend=nullptr;
    Int* elemrecv=nullptr;       
    Int* elemsendpts=nullptr;
    Int* elemrecvpts=nullptr;        
    Int* stgib=nullptr;
    Int *vindx=nullptr;
    
    dstype timing[128];
    dstype* dt=nullptr;
    dstype* dae_dt=nullptr;
    dstype* DIRKcoeff_c=nullptr;
    dstype* DIRKcoeff_d=nullptr;
    dstype* DIRKcoeff_t=nullptr;
    dstype* BDFcoeff_c=nullptr;
    dstype* BDFcoeff_t=nullptr;
    
    cudaEvent_t eventHandle;
    cublasHandle_t cublasHandle;
    
#ifdef  HAVE_MPI
    MPI_Request * requests;
    MPI_Status * statuses;
#endif
    
    void printinfo()
    {
      printf("--------------- Common Struct Information ----------------\n");
      printf("backend: %d\n", backend);   
      printf("number of MPI ranks: %d\n", mpiProcs);   
      printf("number of models: %d\n", nomodels);               
      printf("number of compoments of (u, q): %d\n", nc);   
      printf("number of compoments of u: %d\n", ncu);   
      printf("number of compoments of q: %d\n", ncq);   
      printf("number of compoments of w: %d\n", ncw);   
      printf("number of compoments of v: %d\n", nco);   
      printf("number of compoments of uhat: %d\n", nch);   
      printf("number of compoments of x: %d\n", ncx);   
      printf("number of compoments of s: %d\n", ncs);   
      printf("number of compoments of outputs: %d\n", nce);    
      printf("spatial dimension: %d\n", nd);   
      printf("spatial scheme: %d\n", spatialScheme);        
      printf("element type: %d\n", elemtype);   
      printf("node type: %d\n", 1);   
      printf("polynomial degree: %d\n", porder);   
      printf("gauss quadrature degree: %d\n", pgauss); 
      printf("number of nodes on master element: %d\n", npe); 
      printf("number of gauss points on master element: %d\n", nge); 
      printf("number of nodes on master face: %d\n", npf); 
      printf("number of gauss points on master face: %d\n", ngf); 
      printf("temporal scheme: %d\n", temporalScheme);   
      printf("temporal order: %d\n", torder);   
      printf("number of DIRK stages: %d\n", tstages);   
      printf("number of time steps: %d\n", tsteps);   
      
      printf("total number of elements: %d\n", ne);   
      printf("number of interior elements: %d\n", ne0);   
      printf("number of interior+interface elements: %d\n", ne1);   
      printf("number of interior+interface+exterior elements: %d\n", ne2);   
      printf("total number of faces: %d\n", nf);   
      printf("number of interior faces: %d\n", nf0);   
      
      printf("number of faces per elements: %d\n", nfe);
      printf("number of blocks for elements: %d\n", nbe);
      printf("number of blocks for faces: %d\n", nbf);        
      printf("maximum number of faces per block: %d\n", nfb);
      printf("number of blocks for interior elements: %d\n", nbe0);
      printf("number of blocks for interior+interface elements: %d\n", nbe1);
      printf("number of blocks for interior+interface+exterior elements: %d\n", nbe2);
      printf("number of blocks for interior faces: %d\n", nbf0);
      printf("number of blocks for interior+interface faces: %d\n", nbf1);
      printf("number of interface faces: %d\n", ninterfacefaces);

      printf("number of degrees of freedom of u: %d\n", ndof);   
      printf("number of degrees of freedom of q: %d\n", ndofq);   
      printf("number of degrees of freedom of w: %d\n", ndofw);   
      printf("number of degrees of freedom of uhat: %d\n", ndofuhat);   
      printf("number of degrees of freedom of udg: %d\n", ndofudg);   
      printf("number of degrees of freedom of sdg: %d\n", ndofsdg);   
      printf("number of degrees of freedom of odg: %d\n", ndofodg);   
      printf("number of degrees of freedom of edg: %d\n", ndofedg);   
      printf("length of the stabilization: %d\n", ntau);   

      printf("maximum dimension of the reduced basis space: %d\n", RBdim);
      printf("current dimension of the reduced basis space: %d\n", RBcurrentdim);
      printf("the vector to be removed from the RB space and replaced with new vector: %d\n", RBremovedind);
     
      printf("external uhat function flag: %d\n", extUhat);
      printf("external fhat function flag: %d\n", extFhat);
      printf("external stabilization function flag: %d\n", extStab);
      printf("curved mesh flag: %d\n", curvedMesh);
      printf("debug mode flag: %d\n", debugMode);
      printf("time-dependent problem flag: %d\n", tdep);
      printf("wave problem flag: %d\n", wave);
      printf("linear problem flag: %d\n", linearProblem);
      printf("save solution frequency: %d\n", saveSolFreq);
      printf("save solution option: %d\n", saveSolOpt);
      printf("timestep offset to restart simulation: %d\n", timestepOffset);
      printf("time-derivative function flag: %d\n", tdfunc);
      printf("source function flag: %d\n", source);
      printf("model number: %d\n", modelnumber);
      printf("boundary index to save solution: %d\n", ibs);
      printf("save solution boundary frequency: %d\n", saveSolBouFreq);
      printf("compute time-averaged solution flag: %d\n", compudgavg);
      printf("read time-averaged solution flag: %d\n", readudgavg);
    
      printf("number of components of artificial viscosity: %d\n", ncAV);
      printf("number of artificial viscosity smoothing iterations: %d\n", AVsmoothingIter);
      printf("frozen artificial viscosity flag: %d\n", frozenAVflag);
      printf("linear solver type: %d\n", linearSolver);
      printf("nonlinear solver type: %d\n", nonlinearSolver);
      printf("maximum linear solver iterations: %d\n", linearSolverMaxIter);
      printf("current linear solver iteration: %d\n", linearSolverIter);
      printf("maximum nonlinear solver iterations: %d\n", nonlinearSolverMaxIter);
      printf("current nonlinear solver iteration: %d\n", nonlinearSolverIter);
      printf("matrix-vector multiplication order: %d\n", matvecOrder);
      printf("GMRES restart parameter: %d\n", gmresRestart);
      printf("GMRES orthogonalization method: %d\n", gmresOrthogMethod);
      printf("preconditioner type: %d\n", preconditioner);
      printf("preconditioner matrix type: %d\n", precMatrixType);
      printf("PTC matrix type: %d\n", ptcMatrixType);
      printf("run mode: %d\n", runmode);
      printf("time step factor: %f\n", dtfactor);
      printf("current simulation time: %f\n", time);
      printf("matrix-vector multiplication tolerance: %f\n", matvecTol);
      printf("linear solver tolerance: %f\n", linearSolverTol);
      printf("linear solver tolerance factor: %f\n", linearSolverTolFactor);
      printf("nonlinear solver tolerance: %f\n", nonlinearSolverTol);
      printf("linear solver relative error: %f\n", linearSolverRelError);
      printf("artificial viscosity ramp factor: %f\n", rampFactor);
      printf("PTC parameter: %f\n", PTCparam);
      printf("initial stabilization parameter: %f\n", tau0);
      printf("DAE alpha parameter: %f\n", dae_alpha);
      printf("DAE beta parameter: %f\n", dae_beta);
      printf("DAE gamma parameter: %f\n", dae_gamma);
      printf("DAE epsilon parameter: %f\n", dae_epsilon);
      
      printf("number of boundary conditions: %d\n", maxnbc);
      printf("number of neighboring subdomains: %d\n", nnbsd);      
      printf("number of elements to send: %d\n", nelemsend);
      printf("number of elements to receive: %d\n", nelemrecv);
      
      printf("eblks array: %d by %d\n", 3, nbe);
      for (int j=0; j<3; j++) {
        for (int i=0; i<nbe; i++)
          printf("%d  ", eblks[j+3*i]);
        printf("\n");  
      }

      printf("fblks array: %d by %d\n", 3, nbf);
      for (int j=0; j<3; j++) {
        for (int i=0; i<nbf; i++)
          printf("%d  ", fblks[j+3*i]);
        printf("\n");  
      }

      if (spatialScheme==1) {
        printf("nboufaces array: %d by %d\n", maxnbc, nbe);
        for (int j=0; j<maxnbc; j++) {
          for (int i=0; i<nbe; i++)
            printf("%d  ", nboufaces[1+j+maxnbc*i]);
          printf("\n");  
        }        
      }
      
      if (nnbsd > 1) {
        printf("nbsd array: %d by %d\n", 1, nnbsd);
        for (int i=0; i<nnbsd; i++)
          printf("%d  ", nbsd[i]);
        printf("\n");        

        printf("elemsendpts array: %d by %d\n", 1, nnbsd);
        for (int i=0; i<nnbsd; i++)
          printf("%d  ", elemsendpts[i]);
        printf("\n");        

        printf("elemrecvpts array: %d by %d\n", 1, nnbsd);
        for (int i=0; i<nnbsd; i++)
          printf("%d  ", elemrecvpts[i]);
        printf("\n");          

        printf("elemsend array: %d by %d\n", 1, nelemsend);
        for (int i=0; i<nelemsend; i++)
          printf("%d  ", elemsend[i]);
        printf("\n");        

        printf("elemrecv array: %d by %d\n", 1, nelemrecv);
        for (int i=0; i<nelemrecv; i++)
          printf("%d  ", elemrecv[i]);
        printf("\n");        
      }
    }
    
    void freememory()
    {
        CPUFREE(eblks); 
        CPUFREE(fblks);
        CPUFREE(nboufaces);
        CPUFREE(nbsd); 
        CPUFREE(elemsend); 
        CPUFREE(elemrecv); 
        CPUFREE(elemsendpts); 
        CPUFREE(elemrecvpts); 
        CPUFREE(stgib); 
        CPUFREE(vindx); 
        CPUFREE(dt); 
        CPUFREE(dae_dt); 
        CPUFREE(DIRKcoeff_c); 
        CPUFREE(DIRKcoeff_d); 
        CPUFREE(DIRKcoeff_t); 
        CPUFREE(BDFcoeff_c); 
        CPUFREE(BDFcoeff_t);         
    }                         
};

#endif  
