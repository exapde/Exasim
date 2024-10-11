void KokkosInitu(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)
{
	Kokkos::parallel_for("Initu", ng, KOKKOS_LAMBDA(const size_t i) {
		int j = i%npe;
		int k = i/npe;
		dstype param1 = param[0];
		dstype param4 = param[3];
		dstype xdg1 = xdg[j+npe*0+npe*ncx*k];
		dstype xdg2 = xdg[j+npe*1+npe*ncx*k];
		dstype xdg3 = xdg[j+npe*2+npe*ncx*k];
		dstype t2 = cos(xdg1);
		dstype t3 = cos(xdg2);
		dstype t4 = cos(xdg3);
		dstype t5 = sin(xdg1);
		dstype t6 = sin(xdg2);
		dstype t7 = t4*t4;
		f[j+npe*0+npe*nce*k] = 1.0;
		f[j+npe*1+npe*nce*k] = -t3*t4*t5;
		f[j+npe*2+npe*nce*k] = t2*t4*t6;
		f[j+npe*3+npe*nce*k] = 0.0;
		f[j+npe*4+npe*nce*k] = ((cos(xdg1*2.0)/1.6E+1+cos(xdg2*2.0)/1.6E+1)*(cos(xdg3*2.0)+2.0)+1.0/(param4*param4)/param1)/(param1-1.0)+((t2*t2)*(t6*t6)*t7)/2.0+((t3*t3)*(t5*t5)*t7)/2.0;
	});
}

