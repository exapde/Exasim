template <typename T>  __device__  void devicegpuFlux(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i<ng) {
		T param1 = param[0];
		T param2 = param[1];
		T param3 = param[2];
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
		T t2 = udg2*udg2;
		T t3 = 1.0/(udg1*udg1);
		T t4 = 1.0/udg1;
		T t5 = 1.0/param2;
		T t6 = t2*t3*(1.0/2.0);
		T t7 = udg3*udg3;
		T t8 = t3*t7*(1.0/2.0);
		T t9 = t6+t8;
		T t23 = t9*udg1;
		T t10 = -t23+udg4;
		T t11 = param1-1.0;
		T t22 = t4*udg3*udg5;
		T t12 = -t22+udg7;
		T t13 = t4*t12;
		T t24 = t4*udg2*udg9;
		T t14 = -t24+udg10;
		T t15 = t4*t14;
		T t16 = t13+t15;
		T t21 = t4*udg2*udg5;
		T t17 = -t21+udg6;
		T t18 = t4*t17*2.0;
		T t29 = t4*udg3*udg9;
		T t19 = -t29+udg11;
		T t20 = t18-t4*t19;
		T t25 = t5*t16;
		T t26 = t4*udg2*udg3;
		T t27 = t25+t26;
		T t28 = t10*t11;
		T t30 = t4*udg4;
		T t31 = t4*t10*t11;
		T t32 = t30+t31;
		T t33 = t4*t17;
		T t34 = t33-t4*t19*2.0;
		T t35 = 1.0/param3;
		T t36 = 1.0/t11;
		f[0*ng+i] = udg2;
		f[1*ng+i] = t28+t2*t4+t5*t20*(2.0/3.0);
		f[2*ng+i] = t27;
		f[3*ng+i] = t32*udg2+t4*t5*t16*udg3+t4*t5*t20*udg2*(2.0/3.0)-param1*t3*t5*t35*t36*(t11*udg1*(-udg8+t9*udg5+udg1*(t3*t12*udg3+t3*t17*udg2))+t10*t11*udg5);
		f[4*ng+i] = udg3;
		f[5*ng+i] = t27;
		f[6*ng+i] = t28+t4*t7-t5*t34*(2.0/3.0);
		f[7*ng+i] = t32*udg3+t4*t5*t16*udg2-t4*t5*t34*udg3*(2.0/3.0)-param1*t3*t5*t35*t36*(t11*udg1*(-udg12+t9*udg9+udg1*(t3*t14*udg2+t3*t19*udg3))+t10*t11*udg9);
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> __global__ void kernelgpuFlux(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	devicegpuFlux(f, xdg, udg, odg, wdg, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

template <typename T> void gpuFlux(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	int blockDim = 256;
	int gridDim = (ng + blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuFlux<<<gridDim, blockDim>>>(f, xdg, udg, odg, wdg, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

template void gpuFlux(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int);
template void gpuFlux(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int);
