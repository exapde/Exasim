template <typename T>  __global__  void kernelgpuInitu(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i<ng) {
		int j = i%npe;
		int k = (i-j)/npe;
		T param5 = param[4];
		T param6 = param[5];
		T param7 = param[6];
		T param8 = param[7];
		T xdg1 = xdg[j+npe*0+npe*ncx*k];
		T xdg2 = xdg[j+npe*1+npe*ncx*k];
		f[j+npe*0+npe*nce*k] = param5;
		f[j+npe*1+npe*nce*k] = param6;
		f[j+npe*2+npe*nce*k] = param7;
		f[j+npe*3+npe*nce*k] = param8;
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
