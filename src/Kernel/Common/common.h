#ifndef __COMMON_H__
#define __COMMON_H__

#define SDOT sdot_
#define SGEMV sgemv_
#define SGEMM sgemm_
#define SGETRF sgetrf_
#define SGETRI sgetri_

#define DDOT ddot_
#define DGEMV dgemv_
#define DGEMM dgemm_

#ifndef HAVE_CUDA    
#define DGETRF dgetrf_
#define DGETRI dgetri_
#endif

#ifdef USE_FLOAT
typedef float dstype;
#define cpuNRM2 SNRM2
#define cpuDOT SDOT
#define cpuAXPY SAXPY
#define cpuGEMV SGEMV
#define cpuGEMM SGEMM
#define cpuGETRF SGETRF
#define cpuGETRI SGETRI
#define cpuTRSM STRSM
#else
typedef double dstype; //  double is default precision 
#define cpuNRM2 DNRM2
#define cpuDOT DDOT
#define cpuAXPY DAXPY
#define cpuGEMV DGEMV
#define cpuGEMM DGEMM
#define cpuGETRF DGETRF
#define cpuGETRI DGETRI
#define cpuTRSM DTRSM
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

#define MKL_INT int

#define CPUFREE(x)                                                           \
{                                                                         \
    if (x != NULL) {                                                      \
        free(x);                                                          \
        x = NULL;                                                         \
    }                                                                     \
}

extern "C" {
    double DNRM2(Int*,double*,Int*);
    double DDOT(Int*,double*,Int*,double*,Int*);
    void DAXPY(Int*,double*,double*,Int*,double*,Int*);
    void DGEMV(char*,Int*,Int*,double*,double*,Int*,double*,Int*,double*,double*,Int*);  
    void DGEMM(char*,char*,Int*,Int*,Int*,double*,double*,Int*,
             double*,Int*,double*,double*,Int*);        
    void DGETRF(Int*,Int*,double*,Int*,Int*,Int*);
    void DGETRI(Int*,double*,Int*,Int*,double*,Int*,Int*);
    void DTRSM(char *, char*, char*, char *, Int *, Int *, double*, double*, Int*,
             double*, Int*);

    float SNRM2(Int*,float*,Int*);  
    float SDOT(Int*,float*,Int*,float*,Int*);
    void SAXPY(Int*,float*,float*,Int*,float*,Int*);
    void SGEMM(char*,char*,Int*,Int*,Int*,float*,float*,Int*,
             float*,Int*,float*,float*,Int*);  
    void SGEMV(char*,Int*,Int*,float*,float*,Int*,float*,Int*,float*,float*,Int*);      
    void SGETRF(Int*,Int*,float*,Int*,Int*,Int*);    
    void SGETRI(Int*,float*,Int*,Int*,float*,Int*,Int*);
    void STRSM(char *, char*, char*, char *, Int *, Int*, float*, float*, Int*,
             float*, Int*);    
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
    if (x != NULL) {                                                      \
        cudaTemplateFree(x);                                              \
        x = NULL;                                                         \
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
#ifdef HAVE_ONETHREAD         
    if (backend == 0)       // One thread CPU
        // allocate the memory on the CPU        
        *data = (T *) malloc(n*sizeof(T));      
#endif                   
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        // allocate the memory on the CPU        
        *data = (T *) malloc(n*sizeof(T));      
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        // allocate the memory on the GPU            
        CHECK( cudaMalloc( (void**)data, n * sizeof(T) ) );
#endif                  
}

struct appstruct {              
    Int *lsize=NULL;
    Int *nsize=NULL;  // data size
    Int *ndims=NULL;  // dimensions
    Int *flag=NULL;   // flag parameters
    Int *problem=NULL;// problem parameters    
    Int *comm=NULL;   // communication parameters 
    Int *porder=NULL; // polymnomial degrees
    Int *stgib=NULL;
    Int *vindx=NULL;
    
    dstype *uinf=NULL;    // boundary data
    dstype *dt=NULL;      // time steps       
    dstype *dae_dt=NULL;  // dual time steps       
    dstype *factor=NULL;  // factors      
    dstype *physicsparam=NULL; // physical parameters
    dstype *solversparam=NULL; // solvers parameters
    dstype *tau=NULL; // stabilization parameters
    dstype *stgdata=NULL; 
    dstype *stgparam=NULL; 
    
    //dstype time=NULL;     /* current time */
    dstype *fc_u=NULL;    /* factor when discretizing the time derivative of the U equation. Allow scalar field for local time stepping in steady problems? */
    dstype *fc_q=NULL;    /* factor when discretizing the time derivative of the Q equation. Allow scalar field for local time stepping in steady problems? */
    dstype *fc_p=NULL;    /* factor when discretizing the time derivative of the P equation. Allow scalar field for local time stepping in steady problems? */    
    
    dstype *dtcoef_u=NULL;    /* factor when discretizing the time derivative of the U equation. Allow scalar field for local time stepping in steady problems? */
    dstype *dtcoef_q=NULL;    /* factor when discretizing the time derivative of the Q equation. Allow scalar field for local time stepping in steady problems? */
    dstype *dtcoef_p=NULL;    /* factor when discretizing the time derivative of the P equation. Allow scalar field for local time stepping in steady problems? */    
    
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
            CPUFREE(fc_p);
            CPUFREE(dtcoef_u);
            CPUFREE(dtcoef_q);
            CPUFREE(dtcoef_p);
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
            GPUFREE(fc_p);
            GPUFREE(dtcoef_u);
            GPUFREE(dtcoef_q);
            GPUFREE(dtcoef_p);
       }
#endif       
    }
};


struct masterstruct {       
    
    Int *lsize=NULL;
    Int *nsize=NULL;  // data size
    Int *ndims=NULL;  // dimensions
    
    dstype *shapegt=NULL; // element shape functions at Gauss points (transpose)
    dstype *shapegw=NULL; // element shape functions at Gauss points multiplied by Gauss weights
    dstype *shapfgt=NULL; // face shape functions at Gauss points (transpose)
    dstype *shapfgw=NULL; // face shape functions at Gauss points multiplied by Gauss weights    
    dstype *shapent=NULL; // element shape functions at nodes (transpose)
    dstype *shapen=NULL;  // element shape functions at nodes        
    dstype *shapfnt=NULL; // element shape functions at nodes (transpose)
    dstype *shapfn=NULL;  // element shape functions at nodes        
    dstype *xpe=NULL; // nodal points on master element
    dstype *gpe=NULL; // gauss points on master element
    dstype *gwe=NULL; // gauss weighs on master element
    dstype *xpf=NULL; // nodal points on master face
    dstype *gpf=NULL; // gauss points on master face
    dstype *gwf=NULL; // gauss weighs on master face
    
    dstype *shap1dgt=NULL; 
    dstype *shap1dgw=NULL; 
    dstype *shap1dnt=NULL; 
    dstype *shap1dnl=NULL; 
    dstype *xp1d=NULL; // node points on 1D element
    dstype *gp1d=NULL; // gauss points on 1D element
    dstype *gw1d=NULL; // gauss weights on 1D element    
    
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
       }
#endif       
    }    
};
  
struct meshstruct {       
    Int *lsize=NULL;
    Int *nsize=NULL;  // data size
    Int *ndims=NULL;  // dimensions
        
    Int *facecon=NULL;    // face-to-element connectivities 
    Int *eblks=NULL;    // element blocks
    Int *fblks=NULL;    // face blocks    
    Int *nbsd=NULL;
    Int *elemsend=NULL;
    Int *elemrecv=NULL;
    Int *elemsendpts=NULL;
    Int *elemrecvpts=NULL;
    Int *elempart=NULL;
    Int *elempartpts=NULL;
    Int *cgelcon=NULL;
    Int *rowent2elem=NULL;
    Int *cgent2dgent=NULL;
    Int *colent2elem=NULL;
    Int *rowe2f1=NULL;
    Int *cole2f1=NULL;
    Int *ent2ind1=NULL;
    Int *rowe2f2=NULL;
    Int *cole2f2=NULL;
    Int *ent2ind2=NULL;
    
    Int *findxdg1=NULL; 
    Int *findxdg2=NULL; 
    Int *findxdgp=NULL; 
    Int *findudg1=NULL; 
    Int *findudg2=NULL; 
    Int *findudgp=NULL;     
    Int *eindudg1=NULL;     
    Int *eindudgp=NULL;     
    Int *elemsendind=NULL;
    Int *elemrecvind=NULL;
    Int *elemsendodg=NULL;
    Int *elemrecvodg=NULL;
    Int *elemsendudg=NULL;
    Int *elemrecvudg=NULL;
    Int *index=NULL; 
    
    void freememory(Int hostmemory)
    {
       if (hostmemory==1) {
            CPUFREE(lsize);
            CPUFREE(nsize);
            CPUFREE(ndims);    
            CPUFREE(facecon);    // face-to-element connectivities 
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
            CPUFREE(findxdg2);   
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
            CPUFREE(index); 
        }            
#ifdef HAVE_CUDA      
       else {
            GPUFREE(lsize);
            GPUFREE(nsize);
            GPUFREE(ndims);
            GPUFREE(facecon);    // face-to-element connectivities 
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
            GPUFREE(findxdg2);   
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
            GPUFREE(index);   
       }
#endif       
    }        
};

struct solstruct {       
    Int *lsize=NULL;
    Int *nsize=NULL;  // data size
    Int *ndims=NULL;  // dimensions
    
    dstype *xdg=NULL; // spatial coordinates
    dstype *udg=NULL; // solution (u, q, p) 
    dstype *sdg=NULL; // source term due to the previous solution
    dstype *odg=NULL; // auxilary term 
    dstype *wdg=NULL; // dw/dt = u (wave problem)
    dstype *dudg=NULL; // solution (du, dq, dp) 
    dstype *dwdg=NULL; // dw/dt = u (wave problem)
    dstype *uh=NULL;  // uhat  
    dstype *duh=NULL; // duhat
    dstype *elemg=NULL;
    dstype *faceg=NULL;    
    //dstype *udgg=NULL;
    dstype *sdgg=NULL;
    dstype *odgg=NULL;
    dstype *og1=NULL;
    dstype *og2=NULL;    
    dstype *udgavg=NULL; // time-average solution (u, q, p) 
    dstype *wsrc=NULL;   // source term due to the time derivative for DAE equations  
    dstype *wdual=NULL;   // source term due to the dual time derivative for DAE equations  
    dstype** udgarray;
    
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
#endif            
            CPUFREE(uh);  // uhat      
            CPUFREE(elemg); 
            CPUFREE(faceg); 
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
#endif                        
            GPUFREE(uh);  // uhat          
            GPUFREE(elemg); 
            GPUFREE(faceg); 
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
    //dstype *R=NULL;    // shared memory for all residual vectors
    dstype *Rqe=NULL;  // element residual vector for q
    dstype *Rqf=NULL;  // face residual vector for q   
    dstype *Rue=NULL;  // element residual vector for u
    dstype *Ruf=NULL;  // face residual vector for u
    dstype *Rpe=NULL;  // element residual vector for p
    dstype *Rpf=NULL;  // face residual vector for p
    dstype *Rq=NULL;   // residual vector for q     
    dstype *Ru=NULL;   // residual vector for u    
    dstype *Rp=NULL;   // residual vector for p    
    dstype *Rh=NULL;   // residual vector for uhat

    dstype *dRq=NULL;   // residual vector for q     
    dstype *dRu=NULL;   // residual vector for u        
    dstype *dRh=NULL;   // residual vector for uhat
    dstype *dRqe=NULL;  // element residual vector for q
    dstype *dRqf=NULL;  // face residual vector for q   
    dstype *dRue=NULL;  // element residual vector for u
    dstype *dRuf=NULL;  // face residual vector for u
    
    dstype *Mass=NULL; // store the mass matrix
    dstype *Minv=NULL; // store the inverse of the mass matrix
    
    void freememory(Int hostmemory)
    {
       if (hostmemory==1) {            
            CPUFREE(Rq);    
            CPUFREE(Ru);    
            CPUFREE(Rh);    
#ifdef HAVE_ENZYME                   
            CPUFREE(dRq);   
            CPUFREE(dRu);   
            CPUFREE(dRh);   
#endif                                                
            CPUFREE(Mass);      
            CPUFREE(Minv);      
        }            
#ifdef HAVE_CUDA      
       else {
            GPUFREE(Rq);   
            GPUFREE(Ru);   
            GPUFREE(Rh);   
#ifdef HAVE_ENZYME                   
            GPUFREE(dRq);   
            GPUFREE(dRu);   
            GPUFREE(dRh);   
#endif                                    
            GPUFREE(Mass);   
            GPUFREE(Minv);   
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
    dstype *tempn=NULL; 
    dstype *tempg=NULL;
    dstype *buffrecv=NULL;
    dstype *buffsend=NULL;
    
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
    
    dstype *x=NULL; 
    dstype *u=NULL;
    dstype *r=NULL;
    dstype *b=NULL;
    dstype *v=NULL;
     
    // unified memory for GMRES solver
    dstype *tempmem=NULL;
    
    // store PTC matrix
    dstype *PTCmatrix=NULL;
    dstype *PTCmat;

    // for DIRK schemes
    dstype *utmp=NULL;
    //dstype *w=NULL;
    //dstype *wsrc=NULL;     
    dstype *wtmp=NULL;
    
    // previous solutions for time-dependent problems
    dstype *udgprev=NULL; // 
    dstype *udgprev1=NULL; // 
    dstype *udgprev2=NULL; // 
    dstype *udgprev3=NULL; //         
    dstype *wprev=NULL; // 
    dstype *wprev1=NULL; // 
    dstype *wprev2=NULL; // 
    dstype *wprev3=NULL; // 
        
    void freememory(Int hostmemory)
    {
       if (hostmemory==1) {
            CPUFREE(x); 
            CPUFREE(u); 
            CPUFREE(r); 
            CPUFREE(b); 
            CPUFREE(v); 
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
    
    dstype *W=NULL; 
    dstype *U=NULL;
    dstype *V=NULL;
    dstype *R=NULL;
    dstype *C=NULL;
    dstype *H=NULL;
    
    dstype *y;
    dstype *z=NULL;
    
    dstype *Cmat;
    
    //dstype *RBcoef=NULL;     
    Int *ipiv=NULL;
    
    void freememory(Int hostmemory)
    {
       if (hostmemory==1) {
            CPUFREE(W); 
            CPUFREE(U); 
            CPUFREE(V); 
            CPUFREE(R); 
            CPUFREE(C); 
            CPUFREE(H); 
            //CPUFREE(y); 
            CPUFREE(z); 
            CPUFREE(ipiv); 
        }            
#ifdef HAVE_CUDA      
       else {
            GPUFREE(W);
            GPUFREE(U);
            GPUFREE(V);
            GPUFREE(R);
            GPUFREE(C);
            GPUFREE(H);
            //GPUFREE(y); 
            GPUFREE(z); 
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
    Int ndofudg; // number of degrees of freedom of udg
    Int ndofsdg; // number of degrees of freedom of sdg
    Int ndofodg; // number of degrees of freedom of odg
    Int ndofedg; // number of degrees of freedom of edg
    Int ntau; // length of the stabilization
    
    Int ne0;
    Int ne1;
    Int ne2;
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
    
    Int extUhat=0;
    Int extFhat=0;
    Int extStab=0;
    Int curvedMesh;    
    Int debugMode; // 1: save data to binary files for debugging
    Int appname;   /* 0: Euler; 1: Compressible Navier-Stokes; etc. */
    Int tdep;      // 0: steady-state; 1: time-dependent;  
    Int wave;      //     
    Int linearProblem; // 0: nonlinear problem;  1: linear problem
    Int saveSolFreq;
    Int saveSolOpt;
    Int timestepOffset=0;
    Int stgNmode=0;
    Int tdfunc;
    Int source;
    Int modelnumber;
    Int ibs;
    Int saveSolBouFreq=0;
    Int compudgavg=1;
    Int readudgavg=0;
    Int saveResNorm=1;
    
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
            
    Int* eblks=NULL; // element blocks
    Int* fblks=NULL; // face blocks   
    Int* ncarray=NULL;
    
    Int nstgib;
    Int nnbsd;
    Int nelemsend;
    Int nelemrecv;
    Int nvindx;
    Int* nbsd=NULL;
    Int* elemsend=NULL;
    Int* elemrecv=NULL;       
    Int* elemsendpts=NULL;
    Int* elemrecvpts=NULL;        
    Int* stgib=NULL;
    Int *vindx=NULL;
    
    dstype timing[128];
    dstype* dt=NULL;
    dstype* dae_dt=NULL;
    dstype* DIRKcoeff_c=NULL;
    dstype* DIRKcoeff_d=NULL;
    dstype* DIRKcoeff_t=NULL;
    dstype* BDFcoeff_c=NULL;
    dstype* BDFcoeff_t=NULL;
    
    cudaEvent_t eventHandle;
    cublasHandle_t cublasHandle;
    
#ifdef  HAVE_MPI
    MPI_Request * requests;
    MPI_Status * statuses;
#endif
    
    void freememory()
    {
        CPUFREE(eblks); 
        CPUFREE(fblks);
        CPUFREE(ncarray);    
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
