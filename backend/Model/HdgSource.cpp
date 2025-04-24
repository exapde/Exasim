void HdgSource(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Source", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param2 = param[1];
		dstype param4 = param[3];
		dstype param5 = param[4];
		dstype param7 = param[6];
		dstype param8 = param[7];
		dstype param9 = param[8];
		dstype xdg1 = xdg[0*ng+i];
		dstype xdg2 = xdg[1*ng+i];
		{
		f[0*ng+i] = (param2*param5*param7*xdg1*exp(-(param2*param2)*1.0/(param9*param9)*(pow(xdg2-param8/param2,2.0)+xdg1*xdg1)))/param4;
		}
		{
		f_udg[0*ng+i] = 0.0;
		f_udg[1*ng+i] = 0.0;
		f_udg[2*ng+i] = 0.0;
		}
	});
}

