#ifndef __GPUSTGCOMP
#define __GPUSTGCOMP

template <typename T>   
__global__ void gpuTemplateStg(T *U, T *xdg, T *Umean, T *Amean, T *kcute, T *waveno, T *randno, T *STGparam, 
        T time, int M, int N)
{        
    
    T gam = STGparam[3];
    T Minf = STGparam[4];  
    T Uinf = STGparam[7]/STGparam[6];
    T lemax = STGparam[11];
                 
    int m = threadIdx.x + blockIdx.x * blockDim.x;
    while (m < M) { 

        T sn = 0.0;
        T un = 0.0;
        T vn = 0.0;
        T wn = 0.0;

        // loop over wavenumbers
        for (int n=0; n<N; n++) {
            // wavenumber at mode n
            T kn  = waveno[0*N+n];
            T dkn = waveno[1*N+n];

            // random numbers
            T phi = randno[0*N+n]; 
            T dx = randno[1*N+n];
            T dy = randno[2*N+n]; 
            T dz = randno[3*N+n];
            T sigmax  = randno[4*N+n];
            T sigmay  = randno[5*N+n]; 
            T sigmaz = randno[6*N+n];

            // calculate feta
            T feta = exp(-1.0 * (12.0*kn/kcute[2*M+m]) * (12.0*kn/kcute[2*M+m]));

            // calculate fcut             
            T base = 4.0 * fmax(kn - 0.9*kcute[1*M+m], 0.0)/kcute[1*M+m];
            T fcut = exp( -1.0 * base * base * base);

            // calculate Ek
            T kke = kn/kcute[0*M+m];
            T kke2 = kke*kke;
            T Ek = feta*fcut*(kke2 * kke2)/pow((1.0 + 2.4 * kke2),(17.0/6.0));

            // wave amplitude
            T qn = Ek*dkn;
            // sum of N wave amplitudes
            sn = sn + qn;

            // position 
            T rx = (2.0 * 3.141592653589793 / (kn * lemax)) * (xdg[0*M+m] - Uinf*time);
            // angle
            T an = kn*(dx*rx + dy*xdg[1*M+m] + dz*xdg[2*M+m]) + phi;
            // Fourier mode
            T bn = sqrt(qn)*cos(an);

            // fluctuating velocity field
            un = un + sigmax*bn;
            vn = vn + sigmay*bn;
            wn = wn + sigmaz*bn;                    
        }

        // apply scaling
        T scale = 2.0 * sqrt(3.0 / 2.0);
        un = scale*un/sqrt(sn);
        vn = scale*vn/sqrt(sn);
        wn = scale*wn/sqrt(sn);

        // fluctuating velocity field
        T uprime = Amean[0*M+m]*un;
        T vprime = Amean[1*M+m]*un + Amean[2*M+m]*vn;
        T wprime = Amean[3*M+m]*un + Amean[4*M+m]*vn + Amean[5*M+m]*wn;

        // velocity field
        T u = Umean[1*M+m] + uprime;
        T v = Umean[2*M+m] + vprime;
        T w = Umean[3*M+m] + wprime;

        // density and temperature
        T gam1 = gam-1;
        T Tprime = - (gam1*Minf*Minf)*Umean[1*M+m]*uprime;
        T Rprime = - Tprime * (Umean[0*M+m]/Umean[4*M+m]);
        T r = Umean[0*M+m] + Rprime;    
        T Temp = Umean[4*M+m] + Tprime;   

        // conservative variables
        U[0*M+m] = r;
        U[1*M+m] = r * u;
        U[2*M+m] = r * v;
        U[3*M+m] = r * w;
        U[4*M+m] = r * Temp/((gam-1)*gam*Minf*Minf) + 0.5 * r * (u*u + v*v + w*w);   

        m += blockDim.x * gridDim.x;             
    }        
}

template <typename T> void gpuStg(T *U, T *xdg, T *Umean, T *Amean, T *kcute, T *waveno, T *randno, 
            T *STGparam, T time, int M, int N)
{
    int blockDim = 256;
    int gridDim = (M + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateStg<<<gridDim, blockDim>>>(U, xdg, Umean, Amean, kcute, waveno, randno, STGparam, time, M, N);
}

template void gpuStg(double *, double *, double *, double *, double *, double *, double *,
        double *, double, int, int);
template void gpuStg(float *, float *, float *, float *, float *, float *, float *, 
        float *, float, int, int);


template <typename T>   
__global__ void gpuTemplateStgHomoTurb(T *up, T *xdg, T *stgdata, T *uc, T t, int M, int N)
{        
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
            T an = kn*(dx*(xdg[0*M+m]-uc[0]*t) + dy*(xdg[1*M+m]-uc[1]*t) + dz*(xdg[2*M+m]-uc[2]*t)) + phi + omega*t;
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

template <typename T> void gpuStgHomoTurb(T *up, T *xdg, T *stgdata, T *uc, T t, int M, int N)
{
    int blockDim = 256;
    int gridDim = (M + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateStgHomoTurb<<<gridDim, blockDim>>>(up, xdg, stgdata, uc, t, M, N);
}

template void gpuStgHomoTurb(double *, double *, double *, double *, double, int, int);
template void gpuStgHomoTurb(float *, float *, float *, float *, float, int, int);


template <typename T>  
__global__ void gpuTemplateStgHomoTurb2D(T *up, T *xdg, T *stgdata, T *uc, T t, int M, int N)
{                                 
    // loop over each grid point
    int m = threadIdx.x + blockIdx.x * blockDim.x;
    while (m < M) {  
        
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
__global__ void gpuTemplateStgHomoTurb2D(T *up, T *xdg, T *stgdata, T *ui, T *uc, T *c, T t, int M, int N)
{                                 
    // loop over each grid point
    int m = threadIdx.x + blockIdx.x * blockDim.x;
    while (m < M) {  
                
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

        // velocity field
        up[0*M+m] = (ui[1]/ui[0]) + c[0]*un;
        up[1*M+m] = (ui[2]/ui[0]) + c[1]*vn;

        m += blockDim.x * gridDim.x;     
    }        
}

template <typename T> void gpuStgHomoTurb2D(T *up, T *xdg, T *stgdata, T *ui, T *uc, T *c, T t, int M, int N)
{
    int blockDim = 256;
    int gridDim = (M + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateStgHomoTurb2D<<<gridDim, blockDim>>>(up, xdg, stgdata, ui, uc, c, t, M, N);
}


template <typename T>  
__global__ void gpuTemplateStgHomoTurb2D(T *up, T *xdg, T *stgdata, T *ui, T *uc, T *c, T *param, T t, int M, int N)
{                               
    T gam = param[0];
    T Minf = param[3];  
    
    T rm = ui[0];
    T um = ui[1]/rm;
    T vm = ui[2]/rm;
    T velm2 = um*um + vm*vm;    
    T pm = (gam-1)*(ui[3] - 0.5*rm*velm2);
            
    // loop over each grid point
    int m = threadIdx.x + blockIdx.x * blockDim.x;
    while (m < M) {  
        
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

        // pressure
        T pn = pm + 2*um*(c[0]*un)/(gam*Minf*Minf);
        
        // velocity field
        un = um + c[0]*un;
        vn = vm + c[1]*vn;        
                
        // conservative variables
        up[0*M+m] = rm;
        up[1*M+m] = rm * un;
        up[2*M+m] = rm * vn;        
        up[3*M+m] = pn/(gam-1) + 0.5 * rm * (un*un + vn*vn);      

        m += blockDim.x * gridDim.x;                   
    }        
}

template <typename T> void gpuStgHomoTurb2D(T *up, T *xdg, T *stgdata, T *ui, T *uc, T *c, T *param, T t, int M, int N)
{
    int blockDim = 256;
    int gridDim = (M + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateStgHomoTurb2D<<<gridDim, blockDim>>>(up, xdg, stgdata, ui, uc, c, param, t, M, N);
}


template <typename T>  
__global__ void gpuTemplateStgHomoTurb2D2(T *up, T *xdg, T *stgdata, T *ui, T *uc, T *c, T *param, T t, int M, int N)
{                               
    T gam = param[0];
    //T Minf = param[3];  
    
    T rm = ui[0];
    T um = ui[1]/rm;
    T vm = ui[2]/rm;
    T velm2 = um*um + vm*vm;    
    T pm = (gam-1)*(ui[3] - 0.5*rm*velm2);
            
    // loop over each grid point
    int m = threadIdx.x + blockIdx.x * blockDim.x;
    while (m < M) {  
        
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

        // density
        T rn = rm - 2*rm*c[0]*un/um;
        
        // velocity field
        un = um + c[0]*un;
        vn = vm + c[1]*vn;        
                
        // conservative variables
        up[0*M+m] = rn;
        up[1*M+m] = rn * un;
        up[2*M+m] = rn * vn;        
        up[3*M+m] = pm/(gam-1) + 0.5 * rn * (un*un + vn*vn);                    

        m += blockDim.x * gridDim.x;                   
    }        
}

template <typename T> void gpuStgHomoTurb2D2(T *up, T *xdg, T *stgdata, T *ui, T *uc, T *c, T *param, T t, int M, int N)
{
    int blockDim = 256;
    int gridDim = (M + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateStgHomoTurb2D2<<<gridDim, blockDim>>>(up, xdg, stgdata, ui, uc, c, param, t, M, N);
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
__global__ void gpuTemplateStgHomoTurb3D(T *up, T *xdg, T *stgdata, T *ui, T *uc, T *c, T t, int M, int N)
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
            T an = kn*(dx*(xdg[0*M+m]-uc[0]*t) + dy*(xdg[1*M+m]-uc[1]*t) + dz*(xdg[2*M+m]-uc[2]*t)) + phi + omega*t;
            // Fourier mode
            T bn = amp*cos(an);

            // fluctuating velocity field
            un = un + sigmax*bn;
            vn = vn + sigmay*bn;
            wn = wn + sigmaz*bn;                    
        }

        // velocity field
        up[0*M+m] = (ui[1]/ui[0]) + c[0]*un;
        up[1*M+m] = (ui[2]/ui[0]) + c[1]*vn;
        up[2*M+m] = (ui[3]/ui[0]) + c[2]*wn;      

        m += blockDim.x * gridDim.x;        
    }        
}

template <typename T> void gpuStgHomoTurb3D(T *up, T *xdg, T *stgdata, T *ui, T *uc, T *c, T t, int M, int N)
{
    int blockDim = 256;
    int gridDim = (M + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateStgHomoTurb3D<<<gridDim, blockDim>>>(up, xdg, stgdata, ui, uc, c, t, M, N);
}


template <typename T>  
__global__ void gpuTemplateStgHomoTurb3D(T *up, T *xdg, T *stgdata, T *ui, T *uc, T *c, T *param, T t, int M, int N)
{                               
    T gam = param[0];
    T Minf = param[3];  
    
    T rm = ui[0];
    T um = ui[1]/rm;
    T vm = ui[2]/rm;
    T wm = ui[3]/rm;
    T velm2 = um*um + vm*vm + wm*wm;    
    T pm = (gam-1)*(ui[4] - 0.5*rm*velm2);
            
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
            T an = kn*(dx*(xdg[0*M+m]-uc[0]*t) + dy*(xdg[1*M+m]-uc[1]*t) + dz*(xdg[2*M+m]-uc[2]*t)) + phi + omega*t;
            // Fourier mode
            T bn = amp*cos(an);

            // fluctuating velocity field
            un = un + sigmax*bn;
            vn = vn + sigmay*bn;
            wn = wn + sigmaz*bn;                    
        }

        // pressure
        T pn = pm + 2*um*(c[0]*un)/(gam*Minf*Minf);
        
        // velocity field
        un = um + c[0]*un;
        vn = vm + c[1]*vn;
        wn = wm + c[2]*wn;        
                
        // conservative variables
        up[0*M+m] = rm;
        up[1*M+m] = rm * un;
        up[2*M+m] = rm * vn;
        up[3*M+m] = rm * wn;
        up[4*M+m] = pn/(gam-1) + 0.5 * rm * (un*un + vn*vn + wn*wn);       

        m += blockDim.x * gridDim.x;                     
    }        
}

template <typename T> void gpuStgHomoTurb3D(T *up, T *xdg, T *stgdata, T *ui, T *uc, T *c, T *param, T t, int M, int N)
{
    int blockDim = 256;
    int gridDim = (M + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateStgHomoTurb3D<<<gridDim, blockDim>>>(up, xdg, stgdata, ui, uc, c, param, t, M, N);
}


template <typename T>  
__global__ void gpuTemplateStgHomoTurb3D2(T *up, T *xdg, T *stgdata, T *ui, T *uc, T *c, T *param, T t, int M, int N)
{                               
    T gam = param[0];
    //T Minf = param[3];  
    
    T rm = ui[0];
    T um = ui[1]/rm;
    T vm = ui[2]/rm;
    T wm = ui[3]/rm;
    T velm2 = um*um + vm*vm + wm*wm;    
    T pm = (gam-1)*(ui[4] - 0.5*rm*velm2);
            
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
            T an = kn*(dx*(xdg[0*M+m]-uc[0]*t) + dy*(xdg[1*M+m]-uc[1]*t) + dz*(xdg[2*M+m]-uc[2]*t)) + phi + omega*t;
            // Fourier mode
            T bn = amp*cos(an);

            // fluctuating velocity field
            un = un + sigmax*bn;
            vn = vn + sigmay*bn;
            wn = wn + sigmaz*bn;                    
        }

        // density
        T rn = rm - 2*rm*c[0]*un/um;
        
        // velocity field
        un = um + c[0]*un;
        vn = vm + c[1]*vn;        
        wn = wm + c[2]*wn;        
                
        // conservative variables
        up[0*M+m] = rn;
        up[1*M+m] = rn * un;
        up[2*M+m] = rn * vn;        
        up[3*M+m] = rn * wn;
        up[4*M+m] = pm/(gam-1) + 0.5 * rn * (un*un + vn*vn + wn*wn);                            

        m += blockDim.x * gridDim.x;                     
    }        
}

template <typename T> void gpuStgHomoTurb3D2(T *up, T *xdg, T *stgdata, T *ui, T *uc, T *c, T *param, T t, int M, int N)
{
    int blockDim = 256;
    int gridDim = (M + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateStgHomoTurb3D2<<<gridDim, blockDim>>>(up, xdg, stgdata, ui, uc, c, param, t, M, N);
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

template <typename T>  
void gpuStgHomoTurb(T *up, T *xdg, T *stgdata, T *ui, T *uc, T *c, T t, int M, int N, int nd)
{
	if (nd == 1) {
	}
	else if (nd == 2) {
        gpuStgHomoTurb2D(up, xdg, stgdata, ui, uc, c, t, M, N);
    }
    else {
        gpuStgHomoTurb3D(up, xdg, stgdata, ui, uc, c, t, M, N);
    }
}
template void gpuStgHomoTurb(double *, double *, double *, double *, double *, double *, double, int, int, int);
template void gpuStgHomoTurb(float *, float *, float *, float *, float *, float *, float, int, int, int);

template <typename T>  
void gpuStgHomoTurb(T *up, T *xdg, T *stgdata, T *ui, T *uc, T *c, T *param, T t, int M, int N, int nd)
{
	if (nd == 1) {
	}
	else if (nd == 2) {
        gpuStgHomoTurb2D(up, xdg, stgdata, ui, uc, c, param, t, M, N);
    }
    else {
        gpuStgHomoTurb3D(up, xdg, stgdata, ui, uc, c, param, t, M, N);
    }
}
template void gpuStgHomoTurb(double *, double *, double *, double *, double *, double *, double *, double, int, int, int);
template void gpuStgHomoTurb(float *, float *, float *, float *, float *, float *, float *, float, int, int, int);

template <typename T>  
void gpuStgHomoTurb2(T *up, T *xdg, T *stgdata, T *ui, T *uc, T *c, T *param, T t, int M, int N, int nd)
{
	if (nd == 1) {
	}
	else if (nd == 2) {
        gpuStgHomoTurb2D2(up, xdg, stgdata, ui, uc, c, param, t, M, N);
    }
    else {
        gpuStgHomoTurb3D2(up, xdg, stgdata, ui, uc, c, param, t, M, N);
    }
}
template void gpuStgHomoTurb2(double *, double *, double *, double *, double *, double *, double *, double, int, int, int);
template void gpuStgHomoTurb2(float *, float *, float *, float *, float *, float *, float *, float, int, int, int);


#endif

