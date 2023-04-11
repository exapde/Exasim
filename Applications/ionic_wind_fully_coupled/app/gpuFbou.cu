template <typename T>  __device__  void devicegpuFbou1(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i<ng) {
		T param1 = param[0];
		T param8 = param[7];
		T param9 = param[8];
		T param10 = param[9];
		T tau1 = tau[0];
		T xdg1 = xdg[0*ng+i];
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg4 = udg[3*ng+i];
		T udg7 = udg[6*ng+i];
		T udg8 = udg[7*ng+i];
		T udg11 = udg[10*ng+i];
		T udg12 = udg[11*ng+i];
		T uhg2 = uhg[1*ng+i];
		T uhg3 = uhg[2*ng+i];
		T uhg4 = uhg[3*ng+i];
		T nlg1 = nlg[0*ng+i];
		T nlg2 = nlg[1*ng+i];
		T t2 = nlg1*udg8;
		T t3 = nlg2*udg12;
		T t4 = 1.0/param9;
		T t5 = 1.0/param10;
		T t6 = t2+t3;
		f[0*ng+i] = param9*param10*(t6*udg1-param8*udg2*sqrt(udg8*udg8+udg12*udg12))*(1.89E+2/8.0E+2);
		f[1*ng+i] = t6*udg2*(-6.19047619047619E-3)+tau1*(udg2-uhg2);
		f[2*ng+i] = tau1*(udg3-uhg3)-nlg1*xdg1*((udg3*udg8)/1.4E+2+t4*t5*udg7*(8.0E+2/1.89E+2))-nlg2*xdg1*((udg3*udg12)/1.4E+2+t4*t5*udg11*(8.0E+2/1.89E+2));
		f[3*ng+i] = t2*xdg1+t3*xdg1+tau1*(udg4-uhg4);
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
		f[2*ng+i] = 0.0;
		f[3*ng+i] = 0.0;
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
		T param1 = param[0];
		T param9 = param[8];
		T param10 = param[9];
		T tau1 = tau[0];
		T xdg1 = xdg[0*ng+i];
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg4 = udg[3*ng+i];
		T udg5 = udg[4*ng+i];
		T udg6 = udg[5*ng+i];
		T udg7 = udg[6*ng+i];
		T udg8 = udg[7*ng+i];
		T udg9 = udg[8*ng+i];
		T udg10 = udg[9*ng+i];
		T udg11 = udg[10*ng+i];
		T udg12 = udg[11*ng+i];
		T uhg1 = uhg[0*ng+i];
		T uhg2 = uhg[1*ng+i];
		T uhg3 = uhg[2*ng+i];
		T uhg4 = uhg[3*ng+i];
		T nlg1 = nlg[0*ng+i];
		T nlg2 = nlg[1*ng+i];
		T t2 = nlg1*udg8;
		T t3 = nlg2*udg12;
		T t4 = 1.0/param9;
		T t5 = 1.0/param10;
		T t6 = -uhg1;
		T t7 = -uhg2;
		T t8 = -uhg3;
		T t9 = t6+udg1;
		T t10 = t7+udg2;
		T t11 = t8+udg3;
		T t12 = t2+t3;
		T t13 = tanh(t12);
		T t14 = t9*tau1;
		T t15 = t10*tau1;
		T t16 = t11*tau1;
		T t17 = t13*1.0E+3;
		T t18 = t17-1.0;
		f[0*ng+i] = t13*(-t14+nlg1*xdg1*(udg1*udg8+t4*t5*udg5*(8.0E+2/1.89E+2))+nlg2*xdg1*(udg1*udg12+t4*t5*udg9*(8.0E+2/1.89E+2)))*-1.0E+3-t18*(t14-t12*udg1);
		f[1*ng+i] = t13*(-t15+nlg1*xdg1*(udg2*udg8*6.19047619047619E-3+t4*t5*udg6*(8.0E+2/1.89E+2))+nlg2*xdg1*(udg2*udg12*6.19047619047619E-3+t4*t5*udg10*(8.0E+2/1.89E+2)))*-1.0E+3-t18*(t15-t12*udg2*6.19047619047619E-3);
		f[2*ng+i] = t13*(-t16+nlg1*xdg1*((udg3*udg8)/1.4E+2+t4*t5*udg7*(8.0E+2/1.89E+2))+nlg2*xdg1*((udg3*udg12)/1.4E+2+t4*t5*udg11*(8.0E+2/1.89E+2)))*-1.0E+3-t18*(t16-(t12*udg3)/1.4E+2);
		f[3*ng+i] = t2*xdg1+t3*xdg1+tau1*(udg4-uhg4);
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
		T param1 = param[0];
		T param9 = param[8];
		T param10 = param[9];
		T tau1 = tau[0];
		T xdg1 = xdg[0*ng+i];
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg4 = udg[3*ng+i];
		T udg6 = udg[5*ng+i];
		T udg8 = udg[7*ng+i];
		T udg10 = udg[9*ng+i];
		T udg12 = udg[11*ng+i];
		T uhg1 = uhg[0*ng+i];
		T uhg2 = uhg[1*ng+i];
		T uhg3 = uhg[2*ng+i];
		T uhg4 = uhg[3*ng+i];
		T nlg1 = nlg[0*ng+i];
		T nlg2 = nlg[1*ng+i];
		T t2 = nlg1*udg8;
		T t3 = nlg2*udg12;
		T t4 = 1.0/param9;
		T t5 = 1.0/param10;
		T t6 = t2+t3;
		f[0*ng+i] = -t6*udg1+tau1*(udg1-uhg1);
		f[1*ng+i] = tau1*(udg2-uhg2)-nlg1*xdg1*(udg2*udg8*6.19047619047619E-3+t4*t5*udg6*(8.0E+2/1.89E+2))-nlg2*xdg1*(udg2*udg12*6.19047619047619E-3+t4*t5*udg10*(8.0E+2/1.89E+2));
		f[2*ng+i] = t6*udg3*(-1.0/1.4E+2)+tau1*(udg3-uhg3);
		f[3*ng+i] = t2*xdg1+t3*xdg1+tau1*(udg4-uhg4);
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
		T param1 = param[0];
		T param9 = param[8];
		T param10 = param[9];
		T tau1 = tau[0];
		T xdg1 = xdg[0*ng+i];
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg5 = udg[4*ng+i];
		T udg6 = udg[5*ng+i];
		T udg7 = udg[6*ng+i];
		T udg8 = udg[7*ng+i];
		T udg9 = udg[8*ng+i];
		T udg10 = udg[9*ng+i];
		T udg11 = udg[10*ng+i];
		T udg12 = udg[11*ng+i];
		T uhg1 = uhg[0*ng+i];
		T uhg2 = uhg[1*ng+i];
		T uhg3 = uhg[2*ng+i];
		T nlg1 = nlg[0*ng+i];
		T nlg2 = nlg[1*ng+i];
		T t2 = nlg1*udg8;
		T t3 = nlg2*udg12;
		T t4 = 1.0/param9;
		T t5 = 1.0/param10;
		T t6 = -uhg1;
		T t7 = -uhg2;
		T t8 = -uhg3;
		T t9 = t6+udg1;
		T t10 = t7+udg2;
		T t11 = t8+udg3;
		T t12 = t2+t3;
		T t13 = tanh(t12);
		T t14 = t9*tau1;
		T t15 = t10*tau1;
		T t16 = t11*tau1;
		T t17 = t13*1.0E+3;
		T t18 = t17-1.0;
		f[0*ng+i] = t13*(-t14+nlg1*xdg1*(udg1*udg8+t4*t5*udg5*(8.0E+2/1.89E+2))+nlg2*xdg1*(udg1*udg12+t4*t5*udg9*(8.0E+2/1.89E+2)))*-1.0E+3-t18*(t14-t12*udg1);
		f[1*ng+i] = t13*(-t15+nlg1*xdg1*(udg2*udg8*6.19047619047619E-3+t4*t5*udg6*(8.0E+2/1.89E+2))+nlg2*xdg1*(udg2*udg12*6.19047619047619E-3+t4*t5*udg10*(8.0E+2/1.89E+2)))*-1.0E+3-t18*(t15-t12*udg2*6.19047619047619E-3);
		f[2*ng+i] = t13*(-t16+nlg1*xdg1*((udg3*udg8)/1.4E+2+t4*t5*udg7*(8.0E+2/1.89E+2))+nlg2*xdg1*((udg3*udg12)/1.4E+2+t4*t5*udg11*(8.0E+2/1.89E+2)))*-1.0E+3-t18*(t16-(t12*udg3)/1.4E+2);
		f[3*ng+i] = 0.0;
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