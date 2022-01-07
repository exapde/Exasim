#ifndef __GPUSPHERICALHARMONICS
#define __GPUSPHERICALHARMONICS

// core function
template <typename T> __global__ void gpuKernelSphericalHarmonicsBessel(T *Sr, T *Si, T *x, T *y, T *z, 
                T *x0, T *P, T *tmp, T *f, T *fac, T pi, int L, int K, int N)
{                        
    // L                     : the maximum degree of spherical harmonics
    // K                     : the maximum degree of spherical Bessel functions
    // N                     : length of x, y, z
    // x0  [(L+1)*K]         : zeros of pherical Bessel functions
    // P   [(L+1)*(L+2)/2]   : temporary storage for the recurrence formula
    // tmp [(L+1)*(L+2)/2]   : temporary storage for the recurrence formula
    // f   [L+1]             : temporary storage for the recurrence formula
    // fac                   : factorial look-up table
    // Sr  [N*K*(L+1)*(L+1)] : real part of spherical harmonics Bessel functions
    // Si  [N*K*(L+1)*(L+1)] : imag part of spherical harmonics Bessel functions
        
    // Compute spherical harmonics Bessel functions: S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z)
    
    // https://en.wikipedia.org/wiki/Spherical_harmonics
    // Spherical harmonics functions, Y_{lm}(x,y,z), for l = 0,1,...,L and m = -l,...,0,...,l        
    //  l = 0:                0
    //  l = 1:              1 2 3
    //  l = 2:            4 5 6 7 8
    //  l = 3:        9 10 11 12 13 14 15
    //  l = 4:    16 17 18 19 20 21 22 23 24
    //  ....
    //  l = L:  L*L L*L+1 ....   ....  (L+1)*(L+1)-1 
    // The total number of spherical harmonics is (L+1)*(L+1).
    
    // https://en.wikipedia.org/wiki/Bessel_function
    // Spherical Bessel functions, g_{lk}(x,y,z), for l = 0,1,...,L and k = 1,2,...,K
    //  l = 0:    0    1     2    ....   K-1
    //  l = 1:    K   K+1   K+2   ....  2K-1
    //  l = 2:   2K  2K+1  2K+2   ....  3K-1        
    //  l = 3:   3K  3K+1  3K+2   ....  4K-1        
    //  l = 4:   4K  4K+1  4K+2   ....  5K-1        
    //  ....
    //  l = L:   LK  LK+1  LK+2   .... (L+1)K-1        
    // The total number of spherical Bessel functions is K*(L+1).
        
    // Hence, the total number of spherical harmonics Bessel functions is K*(L+1)*(L+1).            
            
    int i = threadIdx.x + blockIdx.x * blockDim.x;    
    while (i < N) { // loop over each neighbor atom                   
        int k, j, jm;
        tmp[0] = 1.0;
        
        // Cartersian -> Spherical
        T r = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
        T the = acos(z[i]/r);
        T phi = atan2(y[i],x[i]);
        
        // Spherical harmonics for l = 0        
        int l = 0;
        T Ylmr = sqrt(1/(4*pi));
        T Ylmi = 0.0;                
        T xr, g;
        // Spherical Bessel for l = 0, see https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions:_jn,_yn 
        for (k=0; k<K; k++) {
            xr = x0[k*(L+1)+l]*r;                        
            g = cos(xr)/xr;            
            j = k*N + i;                      
            Sr[j] = g*Ylmr;  // real part                   
            Si[j] = g*Ylmi;  // imag part       
        }
                
        // Spherical harmonics for l = 1;
        l = 1;
        T costhe = cos(the);
        T a = -sin(the);        
        int m = 0;    
        P[0] = costhe;
        P[1] = a;
        
        T C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m]));
        Ylmr = C*P[0];
        Ylmi = 0.0;
        // Spherical Bessel for l = 1, see https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions:_jn,_yn
        for (k=0; k<K; k++) {
            xr = x0[k*(L+1)+l]*r;                        
            g = cos(xr)/(xr*xr) + sin(xr)/xr;
            j = ((l*l + l + m)*K + k)*N + i;                
            Sr[j] = g*Ylmr; // real part           
            Si[j] = g*Ylmi; // imag part                                    
        }        
        
        m = 1;
        C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m]));
        Ylmr = C*cos(phi)*P[1];
        Ylmi = C*sin(phi)*P[1];                
        // Spherical Bessel for l = 1, see https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions:_jn,_yn
        for (k=0; k<K; k++) {
            xr = x0[k*(L+1)+l]*r;                        
            g = cos(xr)/(xr*xr) + sin(xr)/xr; 
            j = ((l*l + l + m)*K + k)*N + i; // spherical harmonics Bessel functions for m > 0                
            Sr[j] = g*Ylmr; // real part           
            Si[j] = g*Ylmi; // imag part                                                            
            jm = ((l*l + l - m)*K + k)*N + i; // spherical harmonics Bessel functions for m < 0                
            Sr[jm] = -Sr[j]; // real part                      
            Si[jm] =  Si[j]; // imag part                                                                       
        }        
                
        for (l=2; l<=L; l++) {                                        
            // Compute associated Legendre polynomial using recurrence formula
            // https://en.wikipedia.org/wiki/Associated_Legendre_polynomials            
            T Pll = P[l-1];
            tmp[(l-1)] = Pll;
            P[(l-1)] = (2*(l-1)+1)*(costhe*Pll);
            P[l] = (2*(l-1)+1)*(a*Pll);             
            for (m=0; m<l-1; m++) {
                Pll = P[m];
                P[m] = ((2*(l-1)+1)*(costhe*Pll) - (l+m-1)*tmp[m])/(T (l-m));
                tmp[m] = Pll;                
            }
                            
            // Compute spherical harmonics Bessel function at level l for m = 0, 1, ..., l
            for (m=0; m<=l; m++) {
                C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m])); 
                Ylmr = C*(cos(m*phi)*P[m]);
                Ylmi = C*(sin(m*phi)*P[m]);                
                for (k=0; k<K; k++) {
                    // Compute the spherical Bessel functions using recurrence formula
                    xr = x0[k*(L+1)+l]*r;                                            
                    f[0] = -cos(xr)/xr;
                    f[1] = -cos(xr)/(xr*xr) - sin(xr)/xr;                    
                    for (int jk=2; jk<(l+1); jk++) 
                        f[jk] = ((2*(jk+1)-3)/xr)*f[jk-1] - f[jk-2];                                     
                    g = -f[l];                                                              
                    j = ((l*l + l + m)*K + k)*N + i;   
                    Sr[j] = g*Ylmr; // real part            
                    Si[j] = g*Ylmi; // imag part                                                            
                }                
            }        
            
            // Compute the spherical harmonics Bessel functions for m < 0 using symmetry properties
            C = -1.0;
            for (int m=1; m<=l; m++)  {                 
                for (k=0; k<K; k++) {
                    j =  ((l*l + l + m)*K + k)*N + i;  // spherical harmonics Bessel functions for m > 0 
                    jm = ((l*l + l - m)*K + k)*N + i;  // spherical harmonics Bessel functions for m < 0                 
                    Sr[jm] = C*Sr[j]; // real part
                    Si[jm] =-C*Si[j]; // imag part                                                                             
                }                     
                C = C*(-1.0);        
            }
        }
        
        i += blockDim.x * gridDim.x;
    }          
}
template <typename T> void gpuSphericalHarmonicsBessel(T *Sr, T *Si, T *x, T *y, T *z, 
                T *x0, T *P, T *tmp, T *f, T *fac, T pi, int L, int K, int N)
{        
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelSphericalHarmonicsBessel<<<gridDim, blockDim>>>(Sr, Si, x, y, z, 
            x0, P, tmp, f, fac, pi, L, K, N);
}
template void gpuSphericalHarmonicsBessel(double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, double, int, int, int);
template void gpuSphericalHarmonicsBessel(float*, float*, float*, float*, float*, float*, 
        float*, float*, float*, float*, float, int, int, int);

// core function
template <typename T> __global__ void gpuKernelSphericalHarmonicsBesselDeriv(T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, T *x, T *y, T *z, 
                T *x0, T *P, T *tmp, T *f, T *dP, T *dtmp, T *df, T *fac, T pi, int L, int K, int N)
{                        
    // L                      : the maximum degree of spherical harmonics
    // K                      : the maximum degree of spherical Bessel functions
    // N                      : length of x, y, z
    // x0  [(L+1)*K]          : zeros of pherical Bessel functions
    // P   [(L+1)*(L+2)/2]    : temporary storage for the recurrence formula
    // tmp [(L+1)*(L+2)/2]    : temporary storage for the recurrence formula
    // f   [L+1]              : temporary storage for the recurrence formula
    // dP   [(L+1)*(L+2)/2]   : temporary storage for the recurrence formula
    // dtmp [(L+1)*(L+2)/2]   : temporary storage for the recurrence formula
    // df   [L+1]             : temporary storage for the recurrence formula    
    // fac                    : factorial look-up table
    // Srx  [N*K*(L+1)*(L+1)] : x-derivative of real spherical harmonics Bessel functions
    // Six  [N*K*(L+1)*(L+1)] : x-derivative of imag spherical harmonics Bessel functions
    // Sry  [N*K*(L+1)*(L+1)] : y-derivative of real spherical harmonics Bessel functions
    // Siy  [N*K*(L+1)*(L+1)] : y-derivative of imag spherical harmonics Bessel functions
    // Srz  [N*K*(L+1)*(L+1)] : z-derivative of real spherical harmonics Bessel functions
    // Siz  [N*K*(L+1)*(L+1)] : z-derivative of imag spherical harmonics Bessel functions
    
    // Compute partial derivatives of spherical harmonics Bessel functions: S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z)
    
    // https://en.wikipedia.org/wiki/Spherical_harmonics
    // Spherical harmonics functions, Y_{lm}(x,y,z), for l = 0,1,...,L and m = -l,...,0,...,l        
    //  l = 0:                0
    //  l = 1:              1 2 3
    //  l = 2:            4 5 6 7 8
    //  l = 3:        9 10 11 12 13 14 15
    //  l = 4:    16 17 18 19 20 21 22 23 24
    //  ....
    //  l = L:  L*L L*L+1 ....   ....  (L+1)*(L+1)-1 
    // The total number of spherical harmonics is (L+1)*(L+1).
    
    // https://en.wikipedia.org/wiki/Bessel_function
    // Spherical Bessel functions, g_{lk}(x,y,z), for l = 0,1,...,L and k = 1,2,...,K
    //  l = 0:    0    1     2    ....   K-1
    //  l = 1:    K   K+1   K+2   ....  2K-1
    //  l = 2:   2K  2K+1  2K+2   ....  3K-1        
    //  l = 3:   3K  3K+1  3K+2   ....  4K-1        
    //  l = 4:   4K  4K+1  4K+2   ....  5K-1        
    //  ....
    //  l = L:   LK  LK+1  LK+2   .... (L+1)K-1        
    // The total number of spherical Bessel functions is K*(L+1).
            
    // Hence, the total number of spherical harmonics Bessel functions, S_{klm}(x,y,z), is K*(L+1)*(L+1).            
    
    int i = threadIdx.x + blockIdx.x * blockDim.x;    
    while (i < N) { // loop over each neighbor atom                   
        int k, j, jm;
        
        // Cartersian -> Spherical
        T r = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
        T the = acos(z[i]/r);
        T phi = atan2(y[i],x[i]);
                
        // partial derivatives of spherical coordinates
        T r2 = r*r;
        T rxy = sqrt(x[i]*x[i] + y[i]*y[i]);
        T rxy2 = rxy*rxy;
        T rr2 = rxy*r2;
        T Rx = x[i]/r;
        T Ry = y[i]/r;
        T Rz = z[i]/r;
        T Thex = x[i]*z[i]/rr2;
        T They = y[i]*z[i]/rr2;
        T Thez = -rxy/r2;
        T Phix = -y[i]/rxy2;
        T Phiy = x[i]/rxy2;
        T Phiz = 0.0;        
        
        // Spherical harmonics Y_{lm}(x,y,z) and their derivatives for l = 0        
        int l = 0;
        T Ylmr = sqrt(1/(4*pi));
        T Ylmi = 0.0;
        T YlmrThe = 0.0;
        T YlmrPhi = 0.0;
        T YlmiThe = 0.0;        
        T YlmiPhi = 0.0;        
        T Ylmrx = YlmrThe*Thex + YlmrPhi*Phix;
        T Ylmry = YlmrThe*They + YlmrPhi*Phiy;
        T Ylmrz = YlmrThe*Thez + YlmrPhi*Phiz;
        T Ylmix = YlmiThe*Thex + YlmiPhi*Phix;
        T Ylmiy = YlmiThe*They + YlmiPhi*Phiy;
        T Ylmiz = YlmiThe*Thez + YlmiPhi*Phiz;
        
        T xr, g, gR, gx, gy, gz;
        for (k=0; k<K; k++) {
            // https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions:_jn,_yn             
            // Spherical Bessel g_{lk}(x,y,z) and their derivatives for l = 0             
            xr = x0[k*(L+1)+l]*r;                        
            g = cos(xr)/xr;
            gR = x0[k*(L+1)+l]*(-cos(xr)/(xr*xr) - sin(xr)/xr);
            gx = gR*Rx;
            gy = gR*Ry;
            gz = gR*Rz;            
            
            // derivatives of S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z) for l = 0
            j = k*N + i;                
            Srx[j] = gx*Ylmr + g*Ylmrx;
            Sry[j] = gy*Ylmr + g*Ylmry;
            Srz[j] = gz*Ylmr + g*Ylmrz;
            Six[j] = gx*Ylmi + g*Ylmix;
            Siy[j] = gy*Ylmi + g*Ylmiy;
            Siz[j] = gz*Ylmi + g*Ylmiz;            
        }
                        
        l = 1;
        T costhe = cos(the);
        T a = -sin(the);        
        T dcosthe = a;
        T da = -costhe;        
        int m = 0;    
        P[0] = costhe;
        P[1] = a;
        dP[0] = dcosthe;    
        dP[1] = da;
        tmp[0] = 1.0;
        dtmp[0] = 0.0;
        
        // Spherical harmonics Y_{lm}(x,y,z) and their derivatives for l = 1 and m = 0  
        m = 0;
        T C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m]));
        Ylmr = C*P[0];
        Ylmi = 0.0;
        YlmrThe = C*a;    
        YlmiThe = 0.0;    
        YlmrPhi = 0.0;          
        YlmiPhi = 0.0;                  
        Ylmrx = YlmrThe*Thex + YlmrPhi*Phix;
        Ylmry = YlmrThe*They + YlmrPhi*Phiy;
        Ylmrz = YlmrThe*Thez + YlmrPhi*Phiz;
        Ylmix = YlmiThe*Thex + YlmiPhi*Phix;
        Ylmiy = YlmiThe*They + YlmiPhi*Phiy;
        Ylmiz = YlmiThe*Thez + YlmiPhi*Phiz;        
        for (k=0; k<K; k++) {
            // Spherical Bessel g_{lk}(x,y,z) and their derivatives for l = 1
            xr = x0[k*(L+1)+l]*r;                        
            g = cos(xr)/(xr*xr) + sin(xr)/xr;
            gR = x0[k*(L+1)+l]*(cos(xr)/xr - (2*cos(xr))/(xr*xr*xr) - (2*sin(xr))/(xr*xr));
            gx = gR*Rx;
            gy = gR*Ry;
            gz = gR*Rz;            
                        
            // derivatives of S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z) for l = 1 and m = 0
            j = ((l*l + l + m)*K + k)*N + i;                
            Srx[j] = gx*Ylmr + g*Ylmrx;
            Sry[j] = gy*Ylmr + g*Ylmry;
            Srz[j] = gz*Ylmr + g*Ylmrz;
            Six[j] = gx*Ylmi + g*Ylmix;
            Siy[j] = gy*Ylmi + g*Ylmiy;
            Siz[j] = gz*Ylmi + g*Ylmiz;            
        }        
        
        // Spherical harmonics Y_{lm}(x,y,z) and their derivatives for l = 1 and m = 1
        m = 1;
        C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m]));
        Ylmr = C*cos(phi)*P[1];
        Ylmi = C*sin(phi)*P[1];                
        YlmrThe = C*(cos(phi)*da);      
        YlmiThe = C*(sin(phi)*da);      
        YlmrPhi = -C*(sin(phi)*P[1]); 
        YlmiPhi = C*(cos(phi)*P[1]);  
        Ylmrx = YlmrThe*Thex + YlmrPhi*Phix;
        Ylmry = YlmrThe*They + YlmrPhi*Phiy;
        Ylmrz = YlmrThe*Thez + YlmrPhi*Phiz;
        Ylmix = YlmiThe*Thex + YlmiPhi*Phix;
        Ylmiy = YlmiThe*They + YlmiPhi*Phiy;
        Ylmiz = YlmiThe*Thez + YlmiPhi*Phiz;
        for (k=0; k<K; k++) {
            // Spherical Bessel g_{lk}(x,y,z) and their derivatives for l = 1
            xr = x0[k*(L+1)+l]*r;                        
            g = cos(xr)/(xr*xr) + sin(xr)/xr;
            gR = x0[k*(L+1)+l]*(cos(xr)/xr - (2*cos(xr))/(xr*xr*xr) - (2*sin(xr))/(xr*xr));
            gx = gR*Rx;
            gy = gR*Ry;
            gz = gR*Rz;            
            
            // derivatives of S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z) for l = 1 and m = 1
            j = ((l*l + l + m)*K + k)*N + i;                
            Srx[j] = gx*Ylmr + g*Ylmrx;
            Sry[j] = gy*Ylmr + g*Ylmry;
            Srz[j] = gz*Ylmr + g*Ylmrz;
            Six[j] = gx*Ylmi + g*Ylmix;
            Siy[j] = gy*Ylmi + g*Ylmiy;
            Siz[j] = gz*Ylmi + g*Ylmiz;        
            
            // derivatives of S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z) for l = 1 and m = -1
            jm = ((l*l + l - m)*K + k)*N + i;                
            Srx[jm] = -Srx[j];
            Sry[jm] = -Sry[j];
            Srz[jm] = -Srz[j];
            Six[jm] = Six[j];
            Siy[jm] = Siy[j];
            Siz[jm] = Siz[j];            
        }        
                
        for (l=2; l<=L; l++) {        
            // Compute associated Legendre polynomials P_m and their derivatives using recurrence formula
            // https://en.wikipedia.org/wiki/Associated_Legendre_polynomials                        
            T Pll = P[l-1];
            tmp[(l-1)] = Pll;
            P[(l-1)] = (2*(l-1)+1)*(costhe*Pll);
            P[l] = (2*(l-1)+1)*(a*Pll); 
            T dPll = dP[l-1];
            dtmp[l-1] = dPll;
            dP[l-1] = (2*(l-1)+1)*(costhe*dPll + dcosthe*Pll);   
            dP[l] = (2*(l-1)+1)*(a*dPll + da*Pll);             
            for (m=0; m<l-1; m++) {
                Pll = P[m];
                P[m] = ((2*(l-1)+1)*(costhe*Pll) - (l+m-1)*tmp[m])/(T (l-m));
                tmp[m] = Pll;
                dPll = dP[m];
                dP[m] = ((2*(l-1)+1)*(costhe*dPll+dcosthe*Pll) - (l+m-1)*dtmp[m])/(T (l-m));
                dtmp[m] = dPll;
            }
                            
            // Compute spherical harmonics Bessel functions and their derivatives at level l for m = 0, 1, ..., l
            for (m=0; m<=l; m++) {
                // Spherical harmonics Y_{lm}(x,y,z) and their derivatives for l and m = 0,1,..,l  
                C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m])); 
                Ylmr = C*(cos(m*phi)*P[m]);
                Ylmi = C*(sin(m*phi)*P[m]);
                YlmrThe = C*(cos(m*phi)*dP[m]); 
                YlmiThe = C*(sin(m*phi)*dP[m]); 
                YlmrPhi = -(m*C)*(sin(m*phi)*P[m]); 
                YlmiPhi = (m*C)*(cos(m*phi)*P[m]);          
                Ylmrx = YlmrThe*Thex + YlmrPhi*Phix;
                Ylmry = YlmrThe*They + YlmrPhi*Phiy;
                Ylmrz = YlmrThe*Thez + YlmrPhi*Phiz;
                Ylmix = YlmiThe*Thex + YlmiPhi*Phix;
                Ylmiy = YlmiThe*They + YlmiPhi*Phiy;
                Ylmiz = YlmiThe*Thez + YlmiPhi*Phiz; 

                for (k=0; k<K; k++) {
                    // Spherical Bessel g_{lk}(x,y,z) and their derivatives for l >=2
                    xr = x0[k*(L+1)+l]*r;                                            
                    f[0] = -cos(xr)/xr;
                    f[1] = -cos(xr)/(xr*xr) - sin(xr)/xr;
                    df[0] = x0[k*(L+1)+l]*(cos(xr)/(xr*xr) + sin(xr)/xr);
                    df[1] = x0[k*(L+1)+l]*((2*cos(xr))/(xr*xr*xr) - cos(xr)/xr  + (2*sin(xr))/(xr*xr));
                    for (int jk=2; jk<(l+1); jk++) {
                        f[jk] = ((2*(jk+1)-3)/xr)*f[jk-1] - f[jk-2];                 
                        df[jk] = ((2*(jk+1)-3)/xr)*df[jk-1] - x0[k*(L+1)+l]*((2*(jk+1)-3)/(xr*xr))*f[jk-1] - df[jk-2];        
                    }
                    g = -f[l];
                    gR = -df[l];                    
                    gx = gR*Rx;
                    gy = gR*Ry;
                    gz = gR*Rz; 
                    
                    // derivatives of S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z) for l >= 2 and m = 0,1,...,l
                    j = ((l*l + l + m)*K + k)*N + i;                
                    Srx[j] = gx*Ylmr + g*Ylmrx;
                    Sry[j] = gy*Ylmr + g*Ylmry;
                    Srz[j] = gz*Ylmr + g*Ylmrz;
                    Six[j] = gx*Ylmi + g*Ylmix;
                    Siy[j] = gy*Ylmi + g*Ylmiy;
                    Siz[j] = gz*Ylmi + g*Ylmiz;        
                }                
            }
            
            // derivatives of S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z) for l >= 2 and m = -1,-2,...,-l            
            C = -1.0;
            for (int m=0; m<l; m++) {  
                for (k=0; k<K; k++) {
                    j  = ((l*l + l + m + 1)*K + k)*N + i;  // spherical harmonics Bessel functions for m > 0 
                    jm = ((l*l + l - m - 1)*K + k)*N + i;  // spherical harmonics Bessel functions for m < 0                                     
                    Srx[jm] = C*Srx[j]; 
                    Six[jm] =-C*Six[j];
                    Sry[jm] = C*Sry[j]; 
                    Siy[jm] =-C*Siy[j];
                    Srz[jm] = C*Srz[j]; 
                    Siz[jm] =-C*Siz[j];                    
                }                      
                C = C*(-1.0);        
            }
        }
        i += blockDim.x * gridDim.x;
    }                        
}
template <typename T> void gpuSphericalHarmonicsBesselDeriv(T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, 
      T *x, T *y, T *z, T *x0, T *P, T *tmp, T *f, T *dP, T *dtmp, T *df, T *fac, T pi, int L, int K, int N)
{        
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelSphericalHarmonicsBesselDeriv<<<gridDim, blockDim>>>(Srx, Six, Sry, Siy, Srz, Siz, x, y, z, 
            x0, P, tmp, f, dP, dtmp, df, fac, pi, L, K, N);
}
template void gpuSphericalHarmonicsBesselDeriv(double*, double*, double*, double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, double*, double*, double*, double*, double, int, int, int);
template void gpuSphericalHarmonicsBesselDeriv(float*, float*, float*, float*, float*, float*, float*, float*, float*, 
        float*, float*, float*, float*, float*, float*, float*, float*, float, int, int, int);

// core function
template <typename T> __global__ void gpuKernelRadialSphericalHarmonicsSum(T *ar, T *ai, T *Sr, T *Si, 
        int *Nnb, int K, int N, int Na, int N2, int N3)
{                            
    // L                         : the maximum degree of spherical harmonics
    // K                         : the maximum degree of spherical Bessel functions
    // Na                        : number of atoms in the simulation domain 
    // Nnb [Na+1]                : a list containing the number of neighbors for each global atom    
    // N = Nnb[Na]-Nnb[0]      : total number of neighbors
    // Sr [N*K*(L+1)*(L+1)]      : real spherical harmonics Bessel functions
    // Si [N*K*(L+1)*(L+1)]      : imag spherical harmonics Bessel functions        
    // ar [Na*K*(L+1)*(L+1)]     : sum of real spherical harmonics Bessel functions
    // ai [Na*K*(L+1)*(L+1)]     : sum of imag spherical harmonics Bessel functions    
        
    // Sum the spherical harmonics Bessel functions over neighbors    
    
    int tid = threadIdx.x + blockIdx.x * blockDim.x;    
    while (tid < N3) { 
        int q = tid%N2;     // [1, Na*K]           
        int l = (tid-q)/N2; // spherical basis index [1, L2] 
        int n = q%Na;       // atom index n [1, Na] 
        int k = (q-n)/Na;   // radial basis index [1, K]                
        int j = (l*K + k)*Na + n; // index of atom n
        int m = (l*K + k)*N + (Nnb[n]-Nnb[0]); //  starting index of neighbors
        ar[j] = 0.0; 
        ai[j] = 0.0;
        int numnb = Nnb[n+1]-Nnb[n];  // number of neighbors for atom n
        for (int i=0; i<numnb; i++) { // loop over neighbors 
            ar[j] += Sr[m + i]; // sum over neighbors of atom n
            ai[j] += Si[m + i];                                                    
        }                    
        tid += blockDim.x * gridDim.x;
    }                    
}
template <typename T> void gpuRadialSphericalHarmonicsSum(T *ar, T *ai, T *Sr, T *Si, 
        int *Nnb, int Na, int L, int K)
{        
    int L2 = (L+1)*(L+1);
    int N = Nnb[Na]-Nnb[0]; // total number of neighbors 
    int N1 = Na;            // number of atoms in the simulation domain
    int N2 = K*Na;
    int N3 = L2*K*Na;                        
    int blockDim = 256;
    int gridDim = (N3 + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelRadialSphericalHarmonicsSum<<<gridDim, blockDim>>>(ar, ai, Sr, Si, Nnb, K, N, N1, N2, N3);
}
template void gpuRadialSphericalHarmonicsSum(double*, double*, double*, double*, int*, int, int, int);
template void gpuRadialSphericalHarmonicsSum(float*, float*, float*, float*, int *, int, int, int);


// core function
template <typename T> __global__ void gpuKernelRadialSphericalHarmonicsPower(T *p, T *ar, T *ai, int *indk, 
        int L, int K, int K2, int Na, int N2, int N3)
{   
    // L                      : the maximum degree of spherical harmonics
    // K                      : the maximum degree of spherical Bessel functions
    // Na                     : number of atoms in the simulation domain 
    // indk [K*(K+1)/2]       : store indices for symmetric tensor products of spherical Bessel functions
    // ar [Na*K*(L+1)*(L+1)]  : sum of real spherical harmonics Bessel functions
    // ai [Na*K*(L+1)*(L+1)]  : sum of imag spherical harmonics Bessel functions    
    // p [Na*(L+1)*K*(K+1)/2] : power spectrum components
    
    // Compute the power spectrum components for radial spherical harmonics
    
    int tid = threadIdx.x + blockIdx.x * blockDim.x;    
    while (tid < N3) { 
        int q = tid%N2;     // [1, Na*(L+1)]           
        int k = (tid-q)/N2; // spherical basis index [1, K2] 
        int n = q%Na;       // atom index n [1, Na] 
        int l = (q-n)/Na;   // radial basis index [1, (L+1)]                      
        int k2 = indk[k];
        int k1 = indk[K2+k];
        int l1 = (k*(L+1) + l)*Na + n; // global index of atom n
        p[l1] = (T) 0.0;
        for (int m=0; m<(2*l+1); m++) { // loop over magnetic quantum number
            int i1 = ((l*l + m)*K + k1)*Na + n;
            int j1 = ((l*l + m)*K + k2)*Na + n;
            p[l1] += ar[i1]*ar[j1] + ai[i1]*ai[j1]; // power spectrum           
        }                
        tid += blockDim.x * gridDim.x;
    }             
}
template <typename T> void gpuRadialSphericalHarmonicsPower(T *p, T *ar, T *ai, 
        int *indk, int Na, int L, int K)
{        
    int K2 = K*(K+1)/2;        
    int N1 = Na;    // number of atoms in the simulation domain
    int N2 = (L+1)*Na;
    int N3 = K2*(L+1)*Na;                        
    int blockDim = 256;
    int gridDim = (N3 + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelRadialSphericalHarmonicsPower<<<gridDim, blockDim>>>(p, ar, ai, indk, L, K, K2, N1, N2, N3);
}
template void gpuRadialSphericalHarmonicsPower(double*, double*, double*, int*, int, int, int);
template void gpuRadialSphericalHarmonicsPower(float*, float*, float*, int*, int, int, int);

// core function
template <typename T> __global__ void gpuKernelRadialSphericalHarmonicsPowerDeriv(T *px, T *py, T *pz, T *ar, T *ai, T *arx, T *aix, 
       T *ary, T *aiy, T *arz, T *aiz, int *indk, int *Nnb, int L, int K, int N, int K2, int Na, int N2, int N3)
{   
    // L                      : the maximum degree of spherical harmonics
    // K                      : the maximum degree of spherical Bessel functions    
    // N                      : total number of neighbors 
    // Na                     : number of global atoms in the simulation domain 
    // Nnb [Na+1]             : a list containing the number of neighbors for each global atom    
    // indk [K*(K+1)/2]       : store indices for symmetric tensor products of spherical Bessel functions
    // ar   [Na*K*(L+1)*(L+1)]: sum of real spherical harmonics Bessel functions
    // ai   [Na*K*(L+1)*(L+1)]: sum of imag spherical harmonics Bessel functions    
    // arx  [N*K*(L+1)*(L+1)] : x-derivative of real spherical harmonics Bessel functions
    // aix  [N*K*(L+1)*(L+1)] : x-derivative of imag spherical harmonics Bessel functions
    // ary  [N*K*(L+1)*(L+1)] : y-derivative of real spherical harmonics Bessel functions
    // aiy  [N*K*(L+1)*(L+1)] : y-derivative of imag spherical harmonics Bessel functions
    // arz  [N*K*(L+1)*(L+1)] : z-derivative of real spherical harmonics Bessel functions
    // aiz  [N*K*(L+1)*(L+1)] : z-derivative of imag spherical harmonics Bessel functions
    // p [(L+1)*K*(K+1)/2]    : power spectrum components
    // px [N*(L+1)*K*(K+1)/2] : x-derivative of power spectrum components
    // py [N*(L+1)*K*(K+1)/2] : y-derivative of power spectrum components
    // pz [N*(L+1)*K*(K+1)/2] : z-derivative of power spectrum components
    
    // Compute partial derivatives of the power spectrum components
    
    int tid = threadIdx.x + blockIdx.x * blockDim.x;    
    while (tid < N3) { 
        int q = tid%N2;     // [1, Na*(L+1)]           
        int k = (tid-q)/N2; // spherical basis index [1, K2] 
        int n = q%Na;       // atom index n [1, Na] 
        int l = (q-n)/Na;   // radial basis index [1, (L+1)]                      
        int k2 = indk[k];
        int k1 = indk[K2+k];

        int numnb = Nnb[n+1]-Nnb[n];  // number of neighbors for atom n
        for (int i=0; i<numnb; i++) {// loop over each neighbor of atom n                         
            int j = (k*(L+1) + l)*N + Nnb[n] + i; // index of px, py, pz
            px[j] = (T) 0.0;
            py[j] = (T) 0.0;
            pz[j] = (T) 0.0;
            for (int m=0; m<(2*l+1); m++) {
                int i1 = ((l*l + m)*K + k1)*Na + n; // index of ar and ai 
                int j1 = ((l*l + m)*K + k2)*Na + n; // index of ar and ai                                         
                int i2 = ((l*l + m)*K + k1)*N + Nnb[n] + i;  // index of arx and aix     
                int j2 = ((l*l + m)*K + k2)*N + Nnb[n] + i;  // index of arx and aix                                                    
                px[j] += ar[i1]*arx[j2] + arx[i2]*ar[j1] + ai[i1]*aix[j2] + aix[i2]*ai[j1];
                py[j] += ar[i1]*ary[j2] + ary[i2]*ar[j1] + ai[i1]*aiy[j2] + aiy[i2]*ai[j1];
                pz[j] += ar[i1]*arz[j2] + arz[i2]*ar[j1] + ai[i1]*aiz[j2] + aiz[i2]*ai[j1];                    
            }
        }
        tid += blockDim.x * gridDim.x;
    }                         
}
template <typename T> void gpuRadialSphericalHarmonicsPowerDeriv(T *px, T *py, T *pz, T *ar, T *ai, 
        T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, int *indk, int *Nnb, int Na, int L, int K)
{        
    int N = Nnb[Na]-Nnb[0]; // total number of neighbors 
    int K2 = K*(K+1)/2;        
    int N1 = Na;    // number of atoms in the simulation domain
    int N2 = (L+1)*Na;
    int N3 = K2*(L+1)*Na;                        
    int blockDim = 256;
    int gridDim = (N3 + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelRadialSphericalHarmonicsPowerDeriv<<<gridDim, blockDim>>>(px, py, pz, ar, ai, arx, aix, 
            ary, aiy, arz, aiz, indk, Nnb, L, K, N, K2, N1, N2, N3);
}
template void gpuRadialSphericalHarmonicsPowerDeriv(double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, double*, int*, int*, int, int, int);
template void gpuRadialSphericalHarmonicsPowerDeriv(float*, float*, float*, float*, float*, float*,  
        float*, float*, float*, float*, float*, int*, int*, int, int, int);

// core function
template <typename T> __global__ void gpuKernelRadialSphericalHarmonicsBispectrum(T *b, T *ar, T *ai, T *cg, int *indk, 
        int *indl, int *indm, int *rowm, int Nub, int Ncg, int Na, int L, int K, int K2, int N2, int N3)
{   
    // L                      : the maximum degree of spherical harmonics
    // K                      : the maximum degree of spherical Bessel functions
    // Na                     : number of atoms in the simulation domain    
    // Ncg                    : the total number of non-zero Clebsch-Gordan coefficients 
    // Nub                    : the number of non-zero unique spherical harmonics bispectrum components
    // cg   [Ncg*1]           : store non-zero Clebsch-Gordan coefficients
    // indl [Nub*3]           : store the indices of tuple (l2,l1,l)        
    // indk [K*(K+1)/2]       : store indices for symmetric tensor products of spherical Bessel functions
    // indm [Ncg*3]           : store the indices of tuple (m2, m1, m) for non-zero Clebsch-Gordan coefficients
    // rowm [Nub+1]           : store the number of indices of tuple (m2, m1, m) for each tuple (l2, l1, l)    
    // ar  [Na*K*(L+1)*(L+1)] : sum of real spherical harmonics Bessel functions
    // ai  [Na*K*(L+1)*(L+1)] : sum of imag spherical harmonics Bessel functions
    // b [Na*Nub*K*(K+1)/2]   : bispectrum components    
    
    // Compute the bispectrum components for radial spherical harmonics
    
    int tid = threadIdx.x + blockIdx.x * blockDim.x;    
    while (tid < N3) { 
        int q = tid%N2;     // [1, Na*Nub]           
        int k = (tid-q)/N2; // spherical basis index [1, K2] 
        int n = q%Na;       // atom index n [1, Na] 
        int i = (q-n)/Na;   // radial basis index [1, Nub]                      
        int k2 = indk[k];
        int k1 = indk[K2+k];
        int l2 = indl[i];
        int l1 = indl[Nub+i];
        int l = indl[2*Nub+i];                 
        int nm = rowm[i+1]-rowm[i];
        int indn = (k*Nub + i)*Na + n; // global index of atom n
        b[indn] = (T) 0.0;
        for (int j = 0; j<nm; j++) {
            int m2 = indm[rowm[i]+j];
            int m1 = indm[Ncg+rowm[i]+j];
            int m = indm[2*Ncg + rowm[i]+j];                  
            T a1, b1, a2, b2, a3, b3;                                
            a1 = ar[((l*l + l + m)*K + k1)*Na + n];
            b1 = ai[((l*l + l + m)*K + k1)*Na + n];
            a2 = ar[((l1*l1 + l1 + m1)*K + k2)*Na + n];
            b2 = ai[((l1*l1 + l1 + m1)*K + k2)*Na + n];
            a3 = ar[((l2*l2 + l2 + m2)*K + k2)*Na + n];
            b3 = ai[((l2*l2 + l2 + m2)*K + k2)*Na + n];
            b[indn] += cg[rowm[i]+j]*(a1*a2*a3 + a2*b1*b3 + a3*b1*b2 - a1*b2*b3);                                              
        }
        tid += blockDim.x * gridDim.x;
    }                         
}
template <typename T> void gpuRadialSphericalHarmonicsBispectrum(T *b, T *ar, T *ai, T *cg, int *indk, 
        int *indl, int *indm, int *rowm, int Nub, int Ncg, int Na, int L, int K)
{        
    int K2 = K*(K+1)/2;        
    int N1 = Na;    // number of atoms in the simulation domain
    int N2 = Nub*Na;
    int N3 = K2*Nub*Na;                        
    int blockDim = 256;
    int gridDim = (N3 + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelRadialSphericalHarmonicsBispectrum<<<gridDim, blockDim>>>(b, ar, ai, cg, indk, indl, indm, rowm,
            Nub, Ncg, Na, L, K, K2, N2, N3);
}
template void gpuRadialSphericalHarmonicsBispectrum(double*, double*, double*, double*, int*, int*, int*, int*, int, int, int, int, int);
template void gpuRadialSphericalHarmonicsBispectrum(float*, float*, float*, float*, int*, int*, int*, int*, int, int, int, int, int);

// core function
template <typename T> __global__ void gpuKernelRadialSphericalHarmonicsBispectrumDeriv(T *bx, T *by, T *bz, 
        T *ar, T *ai, T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, T*cg, int *indk, int *indl,
        int *indm, int *rowm, int *Nnb, int Na, int Nub, int Ncg, int K, int N, int K2, int N2, int N3)
{   
    // K                      : the maximum degree of spherical Bessel functions
    // N                      : total number of neighbors 
    // Na                     : number of global atoms in the simulation domain     
    // Ncg                    : the total number of non-zero Clebsch-Gordan coefficients 
    // Nub                    : the number of non-zero unique spherical harmonics bispectrum components
    // Nnb [Na+1]             : a list containing the number of neighbors for each global atom        
    // indl [Nub*3]           : store the indices of tuple (l2,l1,l)        
    // indk [K*(K+1)/2]       : store indices for symmetric tensor products of spherical Bessel functions
    // indm [Ncg*3]           : store the indices of tuple (m2, m1, m) for non-zero Clebsch-Gordan coefficients
    // rowm [Nub+1]           : store the number of indices of tuple (m2, m1, m) for each tuple (l2, l1, l)    
    // cg   [Ncg*1]           : store non-zero Clebsch-Gordan coefficients
    // ar   [Na*K*(L+1)*(L+1)]: real part of spherical harmonics Bessel functions
    // ai   [Na*K*(L+1)*(L+1)]: imag part of spherical harmonics Bessel functions    
    // arx  [N*K*(L+1)*(L+1)] : x-derivative of real spherical harmonics Bessel functions
    // aix  [N*K*(L+1)*(L+1)] : x-derivative of imag spherical harmonics Bessel functions
    // ary  [N*K*(L+1)*(L+1)] : y-derivative of real spherical harmonics Bessel functions
    // aiy  [N*K*(L+1)*(L+1)] : y-derivative of imag spherical harmonics Bessel functions
    // arz  [N*K*(L+1)*(L+1)] : z-derivative of real spherical harmonics Bessel functions
    // aiz  [N*K*(L+1)*(L+1)] : z-derivative of imag spherical harmonics Bessel functions
    // b    [Na*Nub*K*(K+1)/2]: bispectrum components        
    // bx   [N*Nub*K*(K+1)/2] : x-derivative of power spectrum components
    // by   [N*Nub*K*(K+1)/2] : y-derivative of power spectrum components
    // bz   [N*Nub*K*(K+1)/2] : z-derivative of power spectrum components
        
    // Compute partial derivatives of the bispectrum components
    
    int tid = threadIdx.x + blockIdx.x * blockDim.x;    
    while (tid < N3) { 
        int s = tid%N2;     // [1, Na*Nub]           
        int k = (tid-s)/N2; // spherical basis index [1, K2] 
        int n = s%Na;       // atom index n [1, Na] 
        int i = (s-n)/Na;   // radial basis index [1, Nub]                      
        int k2 = indk[k];
        int k1 = indk[K2+k];
        int l2 = indl[i];
        int l1 = indl[Nub+i];
        int l = indl[2*Nub+i];                 
        int nm = rowm[i+1]-rowm[i];
        int numnb = Nnb[n+1]-Nnb[n];  // number of neighbors for atom n
        for (int q=0; q<numnb; q++) {// loop over each neighbor of atom n                         
            int ii = (k*Nub + i)*N + Nnb[n] + q; // index of bx, by, bz                                        
            bx[ii] = (T) 0.0;
            by[ii] = (T) 0.0;
            bz[ii] = (T) 0.0;                                    
            for (int j = 0; j<nm; j++) {
                int m2 = indm[rowm[i]+j];
                int m1 = indm[Ncg+rowm[i]+j];
                int m = indm[2*Ncg + rowm[i]+j];                            
                int n1 = (l*l + l + m)*K + k1;    
                int n2 = (l1*l1 + l1 + m1)*K + k2;
                int n3 = (l2*l2 + l2 + m2)*K + k2;                        
                int lm1 = n1*Na + n; // index of ar and ai 
                int lm2 = n2*Na + n; // index of ar and ai                                         
                int lm3 = n3*Na + n; // index of ar and ai                                                                 
                int mlk1 = n1*N + Nnb[n] + q; // index of arx and aix      
                int mlk2 = n2*N + Nnb[n] + q; // index of arx and aix            
                int mlk3 = n3*N + Nnb[n] + q; // index of arx and aix                                                  

                T c = cg[rowm[i]+j];                    
                T a1, b1, a2, b2, a3, b3;                                
                T a1x, b1x, a2x, b2x, a3x, b3x;
                T a1y, b1y, a2y, b2y, a3y, b3y;
                T a1z, b1z, a2z, b2z, a3z, b3z;
                a1 = ar[lm1];
                b1 = ai[lm1];
                a2 = ar[lm2];
                b2 = ai[lm2];
                a3 = ar[lm3];
                b3 = ai[lm3];
                a1x = arx[mlk1];
                a1y = ary[mlk1];
                a1z = arz[mlk1];
                b1x = aix[mlk1];
                b1y = aiy[mlk1];
                b1z = aiz[mlk1];
                a2x = arx[mlk2];
                a2y = ary[mlk2];
                a2z = arz[mlk2];
                b2x = aix[mlk2];
                b2y = aiy[mlk2];
                b2z = aiz[mlk2];
                a3x = arx[mlk3];
                a3y = ary[mlk3];
                a3z = arz[mlk3];
                b3x = aix[mlk3];
                b3y = aiy[mlk3];
                b3z = aiz[mlk3];

                T t1 = a1x*a2*a3 + a1*a2x*a3 + a1*a2*a3x;                
                T t2 = a2x*b1*b3 + a2*b1x*b3 + a2*b1*b3x;
                T t3 = a3x*b1*b2 + a3*b1x*b2 + a3*b1*b2x;
                T t4 = a1x*b2*b3 + a1*b2x*b3 + a1*b2*b3x;
                bx[ii] += c*(t1 + t2 + t3 - t4);

                t1 = a1y*a2*a3 + a1*a2y*a3 + a1*a2*a3y;                
                t2 = a2y*b1*b3 + a2*b1y*b3 + a2*b1*b3y;
                t3 = a3y*b1*b2 + a3*b1y*b2 + a3*b1*b2y;
                t4 = a1y*b2*b3 + a1*b2y*b3 + a1*b2*b3y;
                by[ii] += c*(t1 + t2 + t3 - t4);

                t1 = a1z*a2*a3 + a1*a2z*a3 + a1*a2*a3z;                
                t2 = a2z*b1*b3 + a2*b1z*b3 + a2*b1*b3z;
                t3 = a3z*b1*b2 + a3*b1z*b2 + a3*b1*b2z;
                t4 = a1z*b2*b3 + a1*b2z*b3 + a1*b2*b3z;
                bz[ii] += c*(t1 + t2 + t3 - t4);                    
            }
        }               
        tid += blockDim.x * gridDim.x;
    }                                     
}
template <typename T> void gpuRadialSphericalHarmonicsBispectrumDeriv(T *bx, T *by, T *bz, 
        T *ar, T *ai, T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, T*cg, int *indk, int *indl,
        int *indm, int *rowm, int *Nnb, int Na, int Nub, int Ncg, int K)
{        
    int N = Nnb[Na]-Nnb[0]; // total number of neighbors 
    int K2 = K*(K+1)/2;        
    int N1 = Na;    // number of atoms in the simulation domain
    int N2 = Nub*Na;
    int N3 = K2*Nub*Na;                        
    int blockDim = 256;
    int gridDim = (N3 + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelRadialSphericalHarmonicsBispectrumDeriv<<<gridDim, blockDim>>>(bx, by, bz, ar, ai, arx, aix, 
            ary, aiy, arz, aiz, cg, indk, indl, indm, rowm, Nnb, Na, Nub, Ncg, K, N, K2, N2, N3);
}
template void gpuRadialSphericalHarmonicsBispectrumDeriv(double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, int, int, int, int);
template void gpuRadialSphericalHarmonicsBispectrumDeriv(float*, float*, float*, float*, float*, float*,  
        float*, float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, int, int, int, int);

// core function
template <typename T> __global__ void gpuKernelRadialSphericalHarmonicsBasis(T *d, T *c, int *atomtype, 
        int Ntype, int Na, int Nbf, int N)
{       
    // Ntype                : number of atom types 
    // Na                   : number of atoms in the simulation domain 
    // Nbf                  : number of basis functions per atom type    
    // atomtype [Ntype]     : a list containing the types of atoms
    // c        [Na*Nbf]    : spectrum components for all atoms
    // d        [Nbf*Ntype] : basis functions based on atom types
    
    int tid = threadIdx.x + blockIdx.x * blockDim.x;    
    while (tid < N) { 
        int m = tid%Nbf;     // [1, Nbf]           
        int t = (tid-m)/Nbf; // [1, Ntype]        
        int k = Nbf*t + m;       // index of the basis function            
        d[k] = (T) 0.0;
        for (int n=0; n<Na; n++) // for each atom n        
        {
            int tn = atomtype[n]; // type of atom n      
            if (tn==t) { // type of atom n is equal to t
                d[k] += c[m*Nbf + n];
            }
        }                        
        tid += blockDim.x * gridDim.x;
    }
}
template <typename T> void gpuRadialSphericalHarmonicsBasis(T *d, T *c, int *atomtype, 
        int Ntype, int Na, int Nbf)
{        
    int N = Nbf*Ntype;
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelRadialSphericalHarmonicsBasis<<<gridDim, blockDim>>>(d, c, atomtype, Ntype, Na, Nbf, N);
}
template void gpuRadialSphericalHarmonicsBasis(double*, double*, int*, int, int, int);
template void gpuRadialSphericalHarmonicsBasis(float*, float*, int*, int, int, int);

// core function
template <typename T> __global__ void gpuKernelRadialSphericalHarmonicsBasisDeriv(T *dx, T *dy, T *dz, T *cx, 
        T *cy, T *cz, int *atomtype, int *neighlist, int *Nnb, int Ntype, int Na, int Nbf, int N, int N2, int N3)
{       
    // Ntype                : number of atom types 
    // Na                   : number of atoms in the simulation domain 
    // Nbf                  : number of basis functions per atom type
    // N = Nnb[Na]-Nnb[0] : total number of neighbors     
    // atomtype [Na]        : a list containing the types of atoms
    // neighborlist [N]     : a list containing the indices of neighbors
    // c [Na*Nbf]           : spectrum components for all atoms
    // d [Nbf*Ntype]        : basis functions based on atom types
    // cx, cy, cz [N*Nbf]   : derivatives of spectrum components  
    // dx, dy, dz [Na*Nbf*Ntype]: derivatives of the basis functions
    
    int tid = threadIdx.x + blockIdx.x * blockDim.x;    
    while (tid < N3) { 
        int s = tid%N2;     // [1, Na*Nbf]           
        int t = (tid-s)/N2; // spherical basis index [1, Ntype] 
        int n = s%Na;       // atom index n [1, Na] 
        int m = (s-n)/Na;   // radial basis index [1, Nbf]       
             
        int nmt = Na*Nbf*t + Na*m + n; //index of the derivatives of the basis function                
        dx[nmt] = (T) 0.0;
        dy[nmt] = (T) 0.0;
        dz[nmt] = (T) 0.0;                
        int tn = atomtype[n]; // type of atom n                
        int numnb = Nnb[n+1]-Nnb[n];  // number of neighbors for atom n                                                   

        for (int q=0; q<numnb; q++) { // loop over each neighbor of atom n                    
            if (tn==t) {// type of atom n is equal to t
                // atom n is self, atom i is neighbor
                int qnm = m*N + (Nnb[n] + q);  // index of the spectrum components of the atom pair (i, n)                    
                dx[nmt] +=  -cx[qnm];
                dy[nmt] +=  -cy[qnm];
                dz[nmt] +=  -cz[qnm]; 
            }                    

            int i = neighlist[Nnb[n] + q]; // atom i 
            int ti = atomtype[i];             // type of atom i                                        
            if (ti==t) { // type of atom i is equal to t                    
                int Nnbi = Nnb[i+1]-Nnb[i];  // number of neighbors for atom i
                int r;
                for (r=0; r<Nnbi; r++) // loop over each neighbor of atom i
                    if (neighlist[Nnb[i] + r] == n) // if a neighbor of atom i matchs atom n
                        break;

                int rim = m*N + (Nnb[i] + r); // index of the spectrum components of the atom pair (n, i)                                        
                // atom n is neighbor, atom i is self
                dx[nmt] += cx[rim];
                dy[nmt] += cy[rim];
                dz[nmt] += cz[rim];
            }
        }                
        tid += blockDim.x * gridDim.x;                
    }
}
template <typename T> void gpuRadialSphericalHarmonicsBasisDeriv(T *dx, T *dy, T *dz, 
        T *cx, T *cy, T *cz, int *atomtype, int *neighlist, int *Nnb, int Ntype, int Na, int Nbf)
{            
    int N = Nnb[Na]-Nnb[0]; // total number of neighbors     
    int N1 = Na;    // number of atoms in the simulation domain
    int N2 = Nbf*Na;
    int N3 = Ntype*Nbf*Na;                            
    int blockDim = 256;
    int gridDim = (N3 + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelRadialSphericalHarmonicsBasisDeriv<<<gridDim, blockDim>>>(dx, dy, dz, cx, cy, cz, atomtype, 
            neighlist, Nnb, Ntype, Na, Nbf, N, N2, N3);
}
template void gpuRadialSphericalHarmonicsBasisDeriv(double*, double*, double*, double*, double*, double*, int*, int*, int*, int, int, int);
template void gpuRadialSphericalHarmonicsBasisDeriv(float*, float*, float*, float*, float*, float*, int*, int*, int*, int, int, int);


#endif


