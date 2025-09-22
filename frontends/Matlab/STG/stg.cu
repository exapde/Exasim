#ifndef __GPUSTGCOMP
#define __GPUSTGCOMP

template <typename T>  
__global__ void gpuTemplateStgHomoTurb2D(T *up, T *xdg, T *stgdata, T *uc, T t, int M, int N)
{                                 
    // loop over each grid point
    int m = threadIdx.x + blockIdx.x * blockDim.x;
    while (m < M) {  

    // for m = 1:M    
        
        T un = 0.0;
        T vn = 0.0;

        // loop over wavenumbers
        for (int n=0; n<N; n++) {
            // wavenumber at mode n
            T kn  = stgdata[0*N+n];
            T amp = stgdata[1*N+n];

            // random numbers
            T phi = stgdata[2*N+n]; 
            T dx = stgdata[3*N+n];
            T dy = stgdata[4*N+n]; 
            T sigmax  = stgdata[6*N+n];
            T sigmay  = stgdata[7*N+n]; 
            T omega = stgdata[9*N+n];
            
            // angle
            T an = kn*(dx*(xdg[0*M+m]-uc[0]*t) + dy*(xdg[1*M+m]-uc[1]*t)) + phi + omega*t;
            // Fourier mode
            T bn = amp*cos(an);

            // fluctuating velocity field
            un = un + sigmax*bn;
            vn = vn + sigmay*bn;
        }

        up[0*M+m] = un;
        up[1*M+m] = vn;          

        m += blockDim.x * gridDim.x;     
    }        
}

template <typename T> void gpuStgHomoTurb2D(T *up, T *xdg, T *stgdata, T *uc, T t, int M, int N)
{
    int blockDim = 256;
    int gridDim = (M + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateStgHomoTurb2D<<<gridDim, blockDim>>>(up, xdg, stgdata, uc, t, M, N);
}

template <typename T>  
__global__ void gpuTemplateStgHomoTurb3D(T *up, T *xdg, T *stgdata, T *uc, T t, int M, int N)
{                                 
    // loop over each grid point
    int m = threadIdx.x + blockIdx.x * blockDim.x;
    while (m < M) {  
        
        T un = 0.0;
        T vn = 0.0;
        T wn = 0.0;

        // loop over wavenumbers
        for (int n=0; n<N; n++) {
            // wavenumber at mode n
            T kn  = stgdata[0*N+n];
            T amp = stgdata[1*N+n];

            // random numbers
            T phi = stgdata[2*N+n]; 
            T dx = stgdata[3*N+n];
            T dy = stgdata[4*N+n]; 
            T dz = stgdata[5*N+n];
            T sigmax  = stgdata[6*N+n];
            T sigmay  = stgdata[7*N+n]; 
            T sigmaz = stgdata[8*N+n];
            T omega = stgdata[9*N+n];
            
            // angle
            //T an = kn*(dx*(xdg[0*M+m]-uc[0]*t) + dy*(xdg[1*M+m]-uc[1]*t) + dz*(xdg[2*M+m]-uc[2]*t)) + phi + omega*t;
            T an = kn*(dx*(xdg[0*M+m]-uc[0]*t) + dy*(xdg[1*M+m]-uc[1]*t) + dz*(xdg[1*M+m]*xdg[2*M+m]-uc[2]*t)) + phi + omega*t;
            // Fourier mode
            T bn = amp*cos(an);

            // fluctuating velocity field
            un = un + sigmax*bn;
            vn = vn + sigmay*bn;
            wn = wn + sigmaz*bn;                    
        }

        up[0*M+m] = un;
        up[1*M+m] = vn;
        up[2*M+m] = wn;        

        m += blockDim.x * gridDim.x;        
    }        
}

template <typename T> void gpuStgHomoTurb3D(T *up, T *xdg, T *stgdata, T *uc, T t, int M, int N)
{
    int blockDim = 256;
    int gridDim = (M + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateStgHomoTurb3D<<<gridDim, blockDim>>>(up, xdg, stgdata, uc, t, M, N);
}

template <typename T>  
void gpuStgHomoTurb(T *up, T *xdg, T *stgdata, T *uc, T t, int M, int N, int nd)
{
	if (nd == 1) {
	}
	else if (nd == 2) {
        gpuStgHomoTurb2D(up, xdg, stgdata, uc, t, M, N);
    }
    else {
        gpuStgHomoTurb3D(up, xdg, stgdata, uc, t, M, N);
    }
}
template void gpuStgHomoTurb(double *, double *, double *, double *, double, int, int, int);
template void gpuStgHomoTurb(float *, float *, float *, float *,  float, int, int, int);

#endif

