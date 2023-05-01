template <typename T>  __device__  void devicegpuFbou1(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i<ng) {
		T tau1 = tau[0];
		T xdg1 = xdg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg4 = udg[3*ng+i];
		T udg6 = udg[5*ng+i];
		T uhg2 = uhg[1*ng+i];
		T odg2 = odg[1*ng+i];
		T odg3 = odg[2*ng+i];
		T nlg1 = nlg[0*ng+i];
		T nlg2 = nlg[1*ng+i];
		f[0*ng+i] = 0.0;
		f[1*ng+i] = tau1*(udg2-uhg2)+nlg1*xdg1*(odg2+udg4)+nlg2*xdg1*(odg3+udg6);
		i += blockDim.x * gridDim.x;
	}
}

template <typename T>  __global__  void kernelgpuFbou1(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	devicegpuFbou1(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

template <typename T>  __device__  void devicegpuFbou2(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i<ng) {
		f[0*ng+i] = 0.0;
		f[1*ng+i] = 0.0;
		i += blockDim.x * gridDim.x;
	}
}

template <typename T>  __global__  void kernelgpuFbou2(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	devicegpuFbou2(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

template <typename T>  __device__  void devicegpuFbou3(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i<ng) {
		T tau1 = tau[0];
		T xdg1 = xdg[0*ng+i];
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg4 = udg[3*ng+i];
		T udg5 = udg[4*ng+i];
		T udg6 = udg[5*ng+i];
		T uhg1 = uhg[0*ng+i];
		T uhg2 = uhg[1*ng+i];
		T odg2 = odg[1*ng+i];
		T odg3 = odg[2*ng+i];
		T nlg1 = nlg[0*ng+i];
		T nlg2 = nlg[1*ng+i];
		T t2 = odg2+udg4;
		T t3 = odg3+udg6;
		T t4 = -uhg1;
		T t5 = t4+udg1;
		T t7 = nlg1*t2*1.0E+3;
		T t8 = nlg2*t3*1.0E+3;
		T t6 = t5*tau1;
		T t9 = t7+t8;
		T t10 = tanh(t9);
		T t11 = t10/2.0;
		f[0*ng+i] = -(t6-udg1*xdg1*(nlg1*t2+nlg2*t3))*(t11-1.0/2.0)+(t11+1.0/2.0)*(t6+nlg1*xdg1*(udg3/1.0E+1-t2*udg1)+nlg2*xdg1*(udg5/1.0E+1-t3*udg1));
		f[1*ng+i] = tau1*(udg2-uhg2)+nlg1*t2*xdg1+nlg2*t3*xdg1;
		i += blockDim.x * gridDim.x;
	}
}

template <typename T>  __global__  void kernelgpuFbou3(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	devicegpuFbou3(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

template <typename T>  __device__  void devicegpuFbou4(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i<ng) {
		T tau1 = tau[0];
		T xdg1 = xdg[0*ng+i];
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg4 = udg[3*ng+i];
		T udg6 = udg[5*ng+i];
		T uhg1 = uhg[0*ng+i];
		T uhg2 = uhg[1*ng+i];
		T odg2 = odg[1*ng+i];
		T odg3 = odg[2*ng+i];
		T nlg1 = nlg[0*ng+i];
		T nlg2 = nlg[1*ng+i];
		T t2 = odg2+udg4;
		T t3 = odg3+udg6;
		f[0*ng+i] = tau1*(udg1-uhg1)-udg1*xdg1*(nlg1*t2+nlg2*t3);
		f[1*ng+i] = tau1*(udg2-uhg2)+nlg1*t2*xdg1+nlg2*t3*xdg1;
		i += blockDim.x * gridDim.x;
	}
}

template <typename T>  __global__  void kernelgpuFbou4(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	devicegpuFbou4(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

template <typename T>  __device__  void devicegpuFbou5(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i<ng) {
		f[0*ng+i] = 0.0;
		f[1*ng+i] = 0.0;
		i += blockDim.x * gridDim.x;
	}
}

template <typename T>  __global__  void kernelgpuFbou5(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	devicegpuFbou5(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

template <typename T> void gpuFbou(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	int blockDim = 256;
	int gridDim = (ng + blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	if (ib == 1)
		kernelgpuFbou1<<<gridDim, blockDim>>>(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 2)
		kernelgpuFbou2<<<gridDim, blockDim>>>(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 3)
		kernelgpuFbou3<<<gridDim, blockDim>>>(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 4)
		kernelgpuFbou4<<<gridDim, blockDim>>>(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 5)
		kernelgpuFbou5<<<gridDim, blockDim>>>(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

template void gpuFbou(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int);
template void gpuFbou(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int);
#ifdef _ENZYME
template <typename T> __global__ void kernelGradgpuFbou1Enzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *dodg, T *wdg, T *dwdg, T *uhg, T *duhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	__enzyme_fwddiffgpuFbou((void*)devicegpuFbou1<T>,
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

template <typename T> __global__ void kernelGradgpuFbou2Enzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *dodg, T *wdg, T *dwdg, T *uhg, T *duhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	__enzyme_fwddiffgpuFbou((void*)devicegpuFbou2<T>,
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

template <typename T> __global__ void kernelGradgpuFbou3Enzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *dodg, T *wdg, T *dwdg, T *uhg, T *duhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	__enzyme_fwddiffgpuFbou((void*)devicegpuFbou3<T>,
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

template <typename T> __global__ void kernelGradgpuFbou4Enzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *dodg, T *wdg, T *dwdg, T *uhg, T *duhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	__enzyme_fwddiffgpuFbou((void*)devicegpuFbou4<T>,
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

template <typename T> __global__ void kernelGradgpuFbou5Enzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *dodg, T *wdg, T *dwdg, T *uhg, T *duhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	__enzyme_fwddiffgpuFbou((void*)devicegpuFbou5<T>,
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

template <typename T> void gpuFbouEnzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *dodg, T *wdg, T *dwdg, T *uhg, T *duhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	int blockDim = 256;
	int gridDim = (ng + blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	if (ib == 1)
		kernelGradgpuFbou1Enzyme<<<gridDim, blockDim>>>(f, df, xg, udg, dudg, odg, dodg, wdg, dwdg, uhg, duhg, nlg, tau, uinf, param, time, modelnumber, ib, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 2)
		kernelGradgpuFbou2Enzyme<<<gridDim, blockDim>>>(f, df, xg, udg, dudg, odg, dodg, wdg, dwdg, uhg, duhg, nlg, tau, uinf, param, time, modelnumber, ib, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 3)
		kernelGradgpuFbou3Enzyme<<<gridDim, blockDim>>>(f, df, xg, udg, dudg, odg, dodg, wdg, dwdg, uhg, duhg, nlg, tau, uinf, param, time, modelnumber, ib, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 4)
		kernelGradgpuFbou4Enzyme<<<gridDim, blockDim>>>(f, df, xg, udg, dudg, odg, dodg, wdg, dwdg, uhg, duhg, nlg, tau, uinf, param, time, modelnumber, ib, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 5)
		kernelGradgpuFbou5Enzyme<<<gridDim, blockDim>>>(f, df, xg, udg, dudg, odg, dodg, wdg, dwdg, uhg, duhg, nlg, tau, uinf, param, time, modelnumber, ib, ng, nc, ncu, nd, ncx, nco, ncw);
}

template void gpuFbouEnzyme(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int);
#endif