void cpuInitu(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)
{
	for (int i = 0; i < ng; i++) {
		int j = i % npe;
		int k = i / npe;
		dstype xdg1 = xdg[j+npe*0+npe*ncx*k];
		dstype xdg2 = xdg[j+npe*1+npe*ncx*k];
		dstype x0 = xdg1 - 0.25;
		f[j+npe*0+npe*nce*k] = exp(-100.0*pow(x0, 2) - 100.0*pow(xdg2, 2));
	}
}

