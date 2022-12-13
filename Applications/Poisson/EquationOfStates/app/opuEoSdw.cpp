template <typename T> void opuEoSdw(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce, int npe, int ne)
{
	for (int i = 0; i <ng; i++) {
		int j = i%npe;
		int k = (i-j)/npe;
		T xdg1 = xdg[j+npe*0+npe*ncx*k];
		T xdg2 = xdg[j+npe*1+npe*ncx*k];
		T udg1 = udg[j+npe*0+npe*nc*k];
		T udg2 = udg[j+npe*1+npe*nc*k];
		T udg3 = udg[j+npe*2+npe*nc*k];
		T wdg1 = wdg[j+npe*0+npe*ncw*k];
		f[j+npe*0+npe*nce*k] = exp(wdg1);
	}
}

template void opuEoSdw(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int, int, int);
template void opuEoSdw(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int, int, int);
