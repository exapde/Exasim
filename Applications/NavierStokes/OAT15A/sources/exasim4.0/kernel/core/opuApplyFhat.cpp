#ifndef __OPUAPPLYFHAT
#define __OPUAPPLYFHAT

template <typename T>  void opuAverageFlux(T *fg, int N)
{	
    for (int tid=0; tid<N; tid++)         
		fg[tid+N] = 0.5*(fg[tid] + fg[tid+N]);            
}
template void opuAverageFlux(double *, int);
template void opuAverageFlux(float *, int);

template <typename T> void opuAverageFluxDotNormal(T *fg, T *nl, int N, int M, int numPoints, int nd)
{	
    for (int tid=0; tid<M; tid++)                 
	{
        int i = tid%numPoints;                
		fg[tid] = fg[N + 0*M + tid] * nl[i + 0 * numPoints];   
        for (int l = 1; l < nd; l++)
            fg[tid] += fg[N + l*M + tid] * nl[i + l * numPoints];
    }
}
template void opuAverageFluxDotNormal(double *, double *, int, int, int, int);
template void opuAverageFluxDotNormal(float *, float *, int, int, int, int);


template <typename T> void opuAddStabilization1(T *fg, T *ug1, T *ug2, T *tau, int M)
{	
    for (int tid=0; tid<M; tid++)                 
		fg[tid] += tau[0] * (ug1[tid] - ug2[tid]);        
}
template void opuAddStabilization1(double *, double *, double *, double *, int);
template void opuAddStabilization1(float *, float *, float *, float *, int);


template <typename T> void opuAddStabilization2(T *fg, T *ug1, T *ug2, T *tau, int M, int numPoints)
{	
    for (int tid=0; tid<M; tid++)                 
    {
        int i = tid%numPoints;   
        int j = (tid-i)/numPoints;
		fg[tid] += tau[j] * (ug1[tid] - ug2[tid]);     
    }
}
template void opuAddStabilization2(double *, double *, double *, double *, int, int);
template void opuAddStabilization2(float *, float *, float *, float *, int, int);


template <typename T>   
void opuAddStabilization3(T *fg, T *ug1, T *ug2, T *tau, int M, int numPoints, int ncu)
{	
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
template void opuAddStabilization3(double *, double *, double *, double *, int, int, int);
template void opuAddStabilization3(float *, float *, float *, float *, int, int, int);


#endif

