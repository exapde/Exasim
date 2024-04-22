template <typename T> void opuInitq(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne)
{
	for (int i = 0; i <ng; i++) {
		int j = i%npe;
		int k = (i-j)/npe;
		T xdg1 = xdg[j+npe*0+npe*ncx*k];
		T xdg2 = xdg[j+npe*1+npe*ncx*k];
		f[j+npe*0+npe*nce*k] = sin(xdg1*3.141592653589793)*sin(xdg2*3.141592653589793);
		f[j+npe*1+npe*nce*k] = 0.0;
		f[j+npe*2+npe*nce*k] = 0.0;
	}
}

template void opuInitq(double *, double *, double *, double *, int, int, int, int, int, int);
template void opuInitq(float *, float *, float *, float *, int, int, int, int, int, int);
