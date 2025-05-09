void cpuInitu(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)
{
	for (int i = 0; i <ng; i++) {
		int j = i%npe;
		int k = i/npe;
		dstype param5 = param[4];
		dstype param6 = param[5];
		dstype param7 = param[6];
		dstype param8 = param[7];
		f[j+npe*0+npe*nce*k] = param5;
		f[j+npe*1+npe*nce*k] = param6;
		f[j+npe*2+npe*nce*k] = param7;
		f[j+npe*3+npe*nce*k] = param8;
	}
}

