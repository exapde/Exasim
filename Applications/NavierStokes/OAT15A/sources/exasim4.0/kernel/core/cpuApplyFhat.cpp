#ifndef __CPUAPPLYFHAT
#define __CPUAPPLYFHAT

template <typename T>  void cpuAverageFlux(T *fg, int N)
{	
    #pragma omp parallel for
    for (int tid=0; tid<N; tid++)         
		fg[tid+N] = 0.5*(fg[tid] + fg[tid+N]);            
}
template void cpuAverageFlux(double *, int);
template void cpuAverageFlux(float *, int);

template <typename T> void cpuAverageFluxDotNormal(T *fg, T *nl, int N, int M, int numPoints, int nd)
{	
    #pragma omp parallel for
    for (int tid=0; tid<M; tid++)                 
	{
        int i = tid%numPoints;                
		fg[tid] = fg[N + 0*M + tid] * nl[i + 0 * numPoints];   
        for (int l = 1; l < nd; l++)
            fg[tid] += fg[N + l*M + tid] * nl[i + l * numPoints];
    }
}
template void cpuAverageFluxDotNormal(double *, double *, int, int, int, int);
template void cpuAverageFluxDotNormal(float *, float *, int, int, int, int);


template <typename T> void cpuAddStabilization1(T *fg, T *ug1, T *ug2, T *tau, int M)
{	
	#pragma omp parallel for
    for (int tid=0; tid<M; tid++)                 
		fg[tid] += tau[0] * (ug1[tid] - ug2[tid]);        
}
template void cpuAddStabilization1(double *, double *, double *, double *, int);
template void cpuAddStabilization1(float *, float *, float *, float *, int);


template <typename T> void cpuAddStabilization2(T *fg, T *ug1, T *ug2, T *tau, int M, int numPoints)
{	
	#pragma omp parallel for
    for (int tid=0; tid<M; tid++)                 
    {
        int i = tid%numPoints;   
        int j = (tid-i)/numPoints;
		fg[tid] += tau[j] * (ug1[tid] - ug2[tid]);     
    }
}
template void cpuAddStabilization2(double *, double *, double *, double *, int, int);
template void cpuAddStabilization2(float *, float *, float *, float *, int, int);


template <typename T>   
void cpuAddStabilization3(T *fg, T *ug1, T *ug2, T *tau, int M, int numPoints, int ncu)
{	
	#pragma omp parallel for
    for (int tid=0; tid<M; tid++)                 
    {
        int i = tid%numPoints;   
        int j = (tid-i)/numPoints;
        for (int k=0; k<ncu; k++) {
            int nm = k * ncu + j;
            int nk = k * numPoints + i;
            fg[tid] += tau[nm] * (ug1[nk] - ug2[nk]);
        }
    }
}
template void cpuAddStabilization3(double *, double *, double *, double *, int, int, int);
template void cpuAddStabilization3(float *, float *, float *, float *, int, int, int);


#endif

