#ifndef __CPUELEMFACEGEOM
#define __CPUELEMFACEGEOM

template <typename T> void cpuElemGeom(T *Xx, T *jac, T *Jg, int ne, int ng, int nd)
{        
    T *Jg11, *Jg12, *Jg13, *Jg21, *Jg22, *Jg23, *Jg31, *Jg32, *Jg33;
    T *Xx11, *Xx12, *Xx13, *Xx21, *Xx22, *Xx23, *Xx31, *Xx32, *Xx33;
    int ngv = ng*ne;
     
    if (nd==1) {
        #pragma omp parallel for
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
        #pragma omp parallel for
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
        #pragma omp parallel for
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

template <typename T> void cpuElemGeom1D(T *jac, T *Xx, T *Jg, int ngv)
{        
    #pragma omp parallel for
    for (int i=0; i<ngv; i++) {
        jac[i] = Jg[i];
        Xx[i] = 1.0;
    }    
}

template <typename T> void cpuElemGeom2D(T *jac, T *Xx11, T *Xx12, T *Xx21, T *Xx22,
                 T *Jg11, T *Jg12, T *Jg21, T *Jg22, int ngv)
{        
    #pragma omp parallel for 
    for (int i=0; i<ngv; i++) {
        jac[i] = Jg11[i]*Jg22[i] - Jg12[i]*Jg21[i];
        Xx11[i] = Jg22[i]; // dxi/dx
        Xx21[i] = -Jg21[i]; //dxi/dy 
        Xx12[i] = -Jg12[i]; //deta/dx
        Xx22[i] = Jg11[i]; //deta/dy
    }        
}

template <typename T> void cpuElemGeom3D(T *jac, T *Xx11, T *Xx12, T *Xx13, T *Xx21, 
                T *Xx22, T *Xx23, T *Xx31, T *Xx32, T *Xx33,
                T *Jg11, T *Jg12, T *Jg13, T *Jg21, T *Jg22, 
                T *Jg23, T *Jg31, T *Jg32, T *Jg33, int ngv)
{            
    #pragma omp parallel for
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

template <typename T> void cpuFaceGeom(T *nlg, T *jacg, T *Jg, int nf, int ng, int nd)
{        
    T *Jg11, *Jg12, *Jg21, *Jg22, *Jg31, *Jg32;
    int na = ng*nf;        
    
    if (nd==1) {
        #pragma omp parallel for
        for (int i=0; i<na; i++) {
            jacg[i] = 1.0;
            nlg[i] = -1.0;
        }                
    }
    else if (nd==2) {        
        #pragma omp parallel for
        for (int i=0; i<na; i++) {
            int j = i+na;
            jacg[i] = sqrt(Jg[i]*Jg[i] + Jg[j]*Jg[j]);
            nlg[i] = Jg[j]/jacg[i];
            nlg[j] = -Jg[i]/jacg[i];
        }
    }
    else if (nd==3) {
        Jg11 = &Jg[0];
        Jg21 = &Jg[na];
        Jg31 = &Jg[2*na];
        Jg12 = &Jg[3*na];
        Jg22 = &Jg[4*na];
        Jg32 = &Jg[5*na];
        #pragma omp parallel for
        for (int i=0; i<na; i++) {
            int j = i+na;
            int k = i+2*na;
            nlg[i] = Jg21[i]*Jg32[i] - Jg31[i]*Jg22[i];
            nlg[j] = Jg31[i]*Jg12[i] - Jg11[i]*Jg32[i];
            nlg[k] = Jg11[i]*Jg22[i] - Jg21[i]*Jg12[i];
            jacg[i] = sqrt(nlg[i]*nlg[i] + nlg[j]*nlg[j] + nlg[k]*nlg[k]);
            nlg[i] = nlg[i]/jacg[i];
            nlg[j] = nlg[j]/jacg[i];
            nlg[k] = nlg[k]/jacg[i];
        }
    }    
}

template <typename T> void cpuFaceGeom1D(T *jacg, T *nlg, T *Jg, int na)
{       
    #pragma omp parallel for
    for (int i=0; i<na; i++) {
        jacg[i] = 1.0;
        nlg[i] = -1.0;
    }                    
}

template <typename T> void cpuFaceGeom2D(T *jacg, T *nlg, T *Jg, int na)
{       
    #pragma omp parallel for
    for (int i=0; i<na; i++) {
        int j = i+na;
        jacg[i] = sqrt(Jg[i]*Jg[i] + Jg[j]*Jg[j]);
        nlg[i] = Jg[j]/jacg[i];
        nlg[j] = -Jg[i]/jacg[i];
    }
}

template <typename T> void cpuFaceGeom3D(T *jacg, T *nlg, T *Jg, int na)
{       
    T *Jg11, *Jg12, *Jg21, *Jg22, *Jg31, *Jg32;
    
    Jg11 = &Jg[0];
    Jg21 = &Jg[na];
    Jg31 = &Jg[2*na];
    Jg12 = &Jg[3*na];
    Jg22 = &Jg[4*na];
    Jg32 = &Jg[5*na];
    #pragma omp parallel for
    for (int i=0; i<na; i++) {
        int j = i+na;
        int k = i+2*na;
        nlg[i] = Jg21[i]*Jg32[i] - Jg31[i]*Jg22[i];
        nlg[j] = Jg31[i]*Jg12[i] - Jg11[i]*Jg32[i];
        nlg[k] = Jg11[i]*Jg22[i] - Jg21[i]*Jg12[i];
        jacg[i] = sqrt(nlg[i]*nlg[i] + nlg[j]*nlg[j] + nlg[k]*nlg[k]);
        nlg[i] = nlg[i]/jacg[i];
        nlg[j] = nlg[j]/jacg[i];
        nlg[k] = nlg[k]/jacg[i];
    }
}

template void cpuElemGeom(double*, double*, double*, int, int, int);
template void cpuFaceGeom(double*, double*, double*, int, int, int);
template void cpuElemGeom1D(double*, double*, double*, int);
template void cpuElemGeom2D(double*, double*, double*, double*, double*, double*, double*, double*, double*, int);
template void cpuElemGeom3D(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*,
                            double*, double*, double*, double*, double*, double*, double*, double*, double*, int);
template void cpuFaceGeom1D(double*, double*, double*, int);
template void cpuFaceGeom2D(double*, double*, double*, int);
template void cpuFaceGeom3D(double*, double*, double*, int);

template void cpuElemGeom(float*, float*, float*, int, int, int);
template void cpuFaceGeom(float*, float*, float*, int, int, int);
template void cpuElemGeom1D(float*, float*, float*, int);
template void cpuElemGeom2D(float*, float*, float*, float*, float*, float*, float*, float*, float*, int);
template void cpuElemGeom3D(float*, float*, float*, float*, float*, float*, float*, float*, float*, float*,
                            float*, float*, float*, float*, float*, float*, float*, float*, float*, int);
template void cpuFaceGeom1D(float*, float*, float*, int);
template void cpuFaceGeom2D(float*, float*, float*, int);
template void cpuFaceGeom3D(float*, float*, float*, int);

#endif


