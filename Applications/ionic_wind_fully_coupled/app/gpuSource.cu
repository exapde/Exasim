template <typename T>  __device__  void devicegpuSource(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i<ng) {
		T param1 = param[0];
		T param5 = param[4];
		T param6 = param[5];
		T param10 = param[9];
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg8 = udg[7*ng+i];
		T udg12 = udg[11*ng+i];
		T t2 = udg8*udg8;
		T t3 = udg12*udg12;
		T t4 = 1.0/param5;
		T t5 = t2+t3;
		T t6 = sqrt(t5);
		f[0*ng+i] = param10*t6*udg1*(-2.2681E-19)+param6*t4*udg1*udg2*5.291005291005291E-12;
		f[1*ng+i] = param10*t6*udg1*1.19E-21+param6*t4*udg2*(udg1+udg3)*5.291005291005291E-12;
		f[2*ng+i] = param10*t6*udg1*2.28E-19+param6*t4*udg2*udg3*5.291005291005291E-12;
		f[3*ng+i] = -udg1+udg2-udg3;
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> __global__ void kernelgpuSource(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	devicegpuSource(f, xdg, udg, odg, wdg, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

template <typename T> void gpuSource(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	int blockDim = 256;
	int gridDim = (ng + blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuSource<<<gridDim, blockDim>>>(f, xdg, udg, odg, wdg, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

template void gpuSource(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int);
template void gpuSource(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int);
