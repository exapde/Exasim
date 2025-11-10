void HdgSource(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Source", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype xdg2 = xdg[1*ng+i];
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg3 = udg[2*ng+i];
		dstype udg4 = udg[3*ng+i];
		{
		dstype t2 = udg3*udg3;
		dstype t3 = 1.0/udg1;
		dstype t5 = 1.0/xdg2;
		dstype t4 = t3*t3;
		f[0*ng+i] = -t5*udg3;
		f[1*ng+i] = -t3*t5*udg2*udg3;
		f[2*ng+i] = -t2*t3*t5;
		f[3*ng+i] = t5*udg3*((t2*t4)/5.0-t3*udg4*(7.0/5.0)+(t4*(udg2*udg2))/5.0);
		}
		{
		dstype t2 = udg2*udg2;
		dstype t3 = udg3*udg3;
		dstype t4 = 1.0/udg1;
		dstype t7 = 1.0/xdg2;
		dstype t5 = t4*t4;
		dstype t6 = t4*t4*t4;
		f_udg[0*ng+i] = 0.0;
		f_udg[1*ng+i] = t5*t7*udg2*udg3;
		f_udg[2*ng+i] = t3*t5*t7;
		f_udg[3*ng+i] = -t7*udg3*(t2*t6*(2.0/5.0)+t3*t6*(2.0/5.0)-t5*udg4*(7.0/5.0));
		f_udg[4*ng+i] = 0.0;
		f_udg[5*ng+i] = -t4*t7*udg3;
		f_udg[6*ng+i] = 0.0;
		f_udg[7*ng+i] = t5*t7*udg2*udg3*(2.0/5.0);
		f_udg[8*ng+i] = -t7;
		f_udg[9*ng+i] = -t4*t7*udg2;
		f_udg[10*ng+i] = t4*t7*udg3*-2.0;
		f_udg[11*ng+i] = t7*((t2*t5)/5.0+(t3*t5)/5.0-t4*udg4*(7.0/5.0))+t3*t5*t7*(2.0/5.0);
		f_udg[12*ng+i] = 0.0;
		f_udg[13*ng+i] = 0.0;
		f_udg[14*ng+i] = 0.0;
		f_udg[15*ng+i] = t4*t7*udg3*(-7.0/5.0);
		f_udg[16*ng+i] = 0.0;
		f_udg[17*ng+i] = 0.0;
		f_udg[18*ng+i] = 0.0;
		f_udg[19*ng+i] = 0.0;
		f_udg[20*ng+i] = 0.0;
		f_udg[21*ng+i] = 0.0;
		f_udg[22*ng+i] = 0.0;
		f_udg[23*ng+i] = 0.0;
		f_udg[24*ng+i] = 0.0;
		f_udg[25*ng+i] = 0.0;
		f_udg[26*ng+i] = 0.0;
		f_udg[27*ng+i] = 0.0;
		f_udg[28*ng+i] = 0.0;
		f_udg[29*ng+i] = 0.0;
		f_udg[30*ng+i] = 0.0;
		f_udg[31*ng+i] = 0.0;
		f_udg[32*ng+i] = 0.0;
		f_udg[33*ng+i] = 0.0;
		f_udg[34*ng+i] = 0.0;
		f_udg[35*ng+i] = 0.0;
		f_udg[36*ng+i] = 0.0;
		f_udg[37*ng+i] = 0.0;
		f_udg[38*ng+i] = 0.0;
		f_udg[39*ng+i] = 0.0;
		f_udg[40*ng+i] = 0.0;
		f_udg[41*ng+i] = 0.0;
		f_udg[42*ng+i] = 0.0;
		f_udg[43*ng+i] = 0.0;
		f_udg[44*ng+i] = 0.0;
		f_udg[45*ng+i] = 0.0;
		f_udg[46*ng+i] = 0.0;
		f_udg[47*ng+i] = 0.0;
		}
	});
}

