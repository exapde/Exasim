template <typename T>  __device__  void devicegpuFbou11(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i<ng) {
		T param1 = param[0];
		T tau1 = tau[0];
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T wdg2 = wdg[1*ng+i];
		T wdg3 = wdg[2*ng+i];
		T nlg1 = nlg[0*ng+i];
		T nlg2 = nlg[1*ng+i];
		f[0*ng+i] = tau1*udg1+nlg1*(param1*udg2+udg1*wdg2)+nlg2*(param1*udg3+udg1*wdg3);
		i += blockDim.x * gridDim.x;
	}
}

template <typename T>  __global__  void kernelgpuFbou11(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	devicegpuFbou11(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

template <typename T> void gpuFbou1(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	int blockDim = 256;
	int gridDim = (ng + blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	if (ib == 1)
		kernelgpuFbou11<<<gridDim, blockDim>>>(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

template void gpuFbou1(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int);
template void gpuFbou1(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int);
#ifdef _ENZYME
template <typename T> __global__ void kernelGradgpuFbou11Enzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *dodg, T *wdg, T *dwdg, T *uhg, T *duhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	__enzyme_fwddiffgpuFbou1((void*)devicegpuFbou11<T>,
			  enzyme_dup, f, df,
			 enzyme_const, xg,
			 enzyme_dup, udg, dudg,
			 enzyme_dup, odg, dodg,
			 enzyme_dup, wdg, dwdg,
			 enzyme_dup, uhg, duhg,
			 enzyme_const, nlg,
			 enzyme_const, tau,
			 enzyme_const, uinf,
			 enzyme_const, param,
			 enzyme_const, time,
			 enzyme_const, modelnumber,
			 enzyme_const, ng,
			 enzyme_const, nc,
			 enzyme_const, ncu,
			 enzyme_const, nd,
			 enzyme_const, ncx,
			 enzyme_const, nco,
			 enzyme_const, ncw);
}

template <typename T> void gpuFbou1Enzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *dodg, T *wdg, T *dwdg, T *uhg, T *duhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	int blockDim = 256;
	int gridDim = (ng + blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	if (ib == 1)
		kernelGradgpuFbou11Enzyme<<<gridDim, blockDim>>>(f, df, xg, udg, dudg, odg, dodg, wdg, dwdg, uhg, duhg, nlg, tau, uinf, param, time, modelnumber, ib, ng, nc, ncu, nd, ncx, nco, ncw);
}

template void gpuFbou1Enzyme(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int);
#endif