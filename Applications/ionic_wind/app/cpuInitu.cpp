template <typename T> void cpuInitu(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne)
{
	#pragma omp parallel for
	for (int i = 0; i <ng; i++) {
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
		T t2 = 1.0/param10;
		f[j+npe*0+npe*nce*k] = (param4*exp(1.0/(param3*param3)*1.0/(t2*t2)*(pow(xdg1-param1*t2,2.0)+pow(xdg2-param2*t2,2.0))*(-1.0/2.0)))/param11;
		f[j+npe*1+npe*nce*k] = 0.0;
	}
}

template void cpuInitu(double *, double *, double *, double *, int, int, int, int, int, int);
template void cpuInitu(float *, float *, float *, float *, int, int, int, int, int, int);
