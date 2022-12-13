template <typename T>  __device__  void devicegpuUbou11(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i<ng) {
		f[0*ng+i] = 0.0;
		i += blockDim.x * gridDim.x;
	}
}

template <typename T>  __global__  void kernelgpuUbou11(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	devicegpuUbou11(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

template <typename T> void gpuUbou1(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	int blockDim = 256;
	int gridDim = (ng + blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	if (ib == 1)
		kernelgpuUbou11<<<gridDim, blockDim>>>(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

template void gpuUbou1(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int);
template void gpuUbou1(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int);
#ifdef _ENZYME
template <typename T> __global__ void kernelGradgpuUbou11Enzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *dodg, T *wdg, T *dwdg, T *uhg, T *duhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	__enzyme_fwddiffgpuUbou1((void*)devicegpuUbou11<T>,
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

template <typename T> void gpuUbou1Enzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *dodg, T *wdg, T *dwdg, T *uhg, T *duhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	int blockDim = 256;
	int gridDim = (ng + blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	if (ib == 1)
		kernelGradgpuUbou11Enzyme<<<gridDim, blockDim>>>(f, df, xg, udg, dudg, odg, dodg, wdg, dwdg, uhg, duhg, nlg, tau, uinf, param, time, modelnumber, ib, ng, nc, ncu, nd, ncx, nco, ncw);
}

template void gpuUbou1Enzyme(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int);
#endif