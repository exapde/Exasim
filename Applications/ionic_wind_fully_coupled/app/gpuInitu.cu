template <typename T>  __global__  void kernelgpuInitu(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i<ng) {
		int j = i%npe;
		int k = (i-j)/npe;
		T param1 = param[0];
		T param2 = param[1];
		T param3 = param[2];
		T param4 = param[3];
		T param10 = param[9];
		T param11 = param[10];
		T xdg1 = xdg[j+npe*0+npe*ncx*k];
		T xdg2 = xdg[j+npe*1+npe*ncx*k];
		T t2 = param10*param10;
		T t3 = 1.0/(param3*param3);
		T t4 = 1.0/param10;
		T t5 = 1.0/param11;
		T t6 = param1*t4;
		T t7 = param2*t4;
		T t12 = pow(t6-xdg1,2.0);
		T t13 = pow(t7-xdg2,2.0);
		T t14 = t12+t13;
		T t15 = (t2*t3*t14)/2.0;
		T t16 = -t15;
		T t17 = exp(t16);
		T t18 = param4*t5*t17;
		f[j+npe*0+npe*nce*k] = t18;
		f[j+npe*1+npe*nce*k] = t18;
		f[j+npe*2+npe*nce*k] = 0.0;
		f[j+npe*3+npe*nce*k] = 0.0;
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> void gpuInitu(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne)
{
	int blockDim = 256;
	int gridDim = (ng + blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuInitu<<<gridDim, blockDim>>>(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
}

template void gpuInitu(double *, double *, double *, double *, int, int, int, int, int, int);
template void gpuInitu(float *, float *, float *, float *, int, int, int, int, int, int);
