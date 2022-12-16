template <typename T> void opuInitu(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne, Mutation::Mixture *mix)
{
	for (int i = 0; i <ng; i++) {
		int j = i%npe;
		int k = (i-j)/npe;
		T param1 = param[0];
		T param2 = param[1];
		T param3 = param[2];
		T param4 = param[3];
		T param5 = param[4];
		T param6 = param[5];
		T param7 = param[6];
		T param8 = param[7];
		T uinf1 = uinf[0];
		T uinf2 = uinf[1];
		T uinf3 = uinf[2];
		T xdg1 = xdg[j+npe*0+npe*ncx*k];
		T xdg2 = xdg[j+npe*1+npe*ncx*k];
		T t2 = 1.0/uinf1;
		T t3 = 1.0/uinf2;
		f[j+npe*0+npe*nce*k] = param1*t2;
		f[j+npe*1+npe*nce*k] = param2*t2;
		f[j+npe*2+npe*nce*k] = param3*t2;
		f[j+npe*3+npe*nce*k] = param4*t2;
		f[j+npe*4+npe*nce*k] = param5*t2;
		f[j+npe*5+npe*nce*k] = param6*t2*t3;
		f[j+npe*6+npe*nce*k] = param7*t2*t3;
		f[j+npe*7+npe*nce*k] = param8/uinf3;
	}
}

template void opuInitu(double *, double *, double *, double *, int, int, int, int, int, int, Mutation::Mixture *);
template void opuInitu(float *, float *, float *, float *, int, int, int, int, int, int, Mutation::Mixture *);
