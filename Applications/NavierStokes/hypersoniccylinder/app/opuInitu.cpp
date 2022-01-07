template <typename T> void opuInitu(T *f, T *xdg, T *uinf, T *param, int ng, int ncx, int nce, int npe, int ne)
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
		T param9 = param[8];
		T param10 = param[9];
		T param11 = param[10];
		T param12 = param[11];
		T param13 = param[12];
		T param14 = param[13];
		T param15 = param[14];
		T param16 = param[15];
		T param17 = param[16];
		T uinf1 = uinf[0];
		T uinf2 = uinf[1];
		T xdg1 = xdg[j+npe*0+npe*ncx*k];
		T xdg2 = xdg[j+npe*1+npe*ncx*k];
		f[j+npe*0+npe*nce*k] = param5;
		f[j+npe*1+npe*nce*k] = param6;
		f[j+npe*2+npe*nce*k] = param7;
		f[j+npe*3+npe*nce*k] = param8;
	}
}

template void opuInitu(double *, double *, double *, double *, int, int, int, int, int);
template void opuInitu(float *, float *, float *, float *, int, int, int, int, int);
