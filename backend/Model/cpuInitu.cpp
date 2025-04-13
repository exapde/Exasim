void cpuInitu(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)
{
	for (int i = 0; i <ng; i++) {
		int j = i%npe;
		int k = i/npe;
		dstype param2 = param[1];
		dstype param3 = param[2];
		dstype xdg1 = xdg[j+npe*0+npe*ncx*k];
		f[j+npe*0+npe*nce*k] = param2/param3;
		f[j+npe*1+npe*nce*k] = exp((xdg1*xdg1)*-1.0E+4);
	}
}

