void HdgFlux(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Flux", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype param2 = param[1];
		dstype udg3 = udg[2*ng+i];
		dstype udg4 = udg[3*ng+i];
		dstype udg5 = udg[4*ng+i];
		dstype udg6 = udg[5*ng+i];
		{
		dstype t2 = udg3*udg6;
		dstype t3 = udg4*udg5;
		dstype t4 = -t2;
		dstype t5 = -t3;
		dstype t6 = t2+t5;
		dstype t7 = t3+t4+1.0;
		dstype t8 = 1.0/t6;
		f[0*ng+i] = param1*udg3-param1*t8*udg6-param2*t7*udg6;
		f[1*ng+i] = param1*udg4+param1*t8*udg5+param2*t7*udg5;
		f[2*ng+i] = param1*udg5+param1*t8*udg4+param2*t7*udg4;
		f[3*ng+i] = param1*udg6-param1*t8*udg3-param2*t7*udg3;
		}
		{
		dstype t2 = udg3*udg6;
		dstype t3 = udg4*udg5;
		dstype t4 = udg3*udg3;
		dstype t5 = udg4*udg4;
		dstype t6 = udg5*udg5;
		dstype t7 = udg6*udg6;
		dstype t8 = param2*udg3*udg4;
		dstype t9 = param2*udg3*udg5;
		dstype t12 = param2*udg4*udg6;
		dstype t13 = param2*udg5*udg6;
		dstype t10 = param2*t2;
		dstype t11 = param2*t3;
		dstype t14 = -t2;
		dstype t15 = -t3;
		dstype t16 = -t8;
		dstype t17 = -t9;
		dstype t18 = -t12;
		dstype t19 = -t13;
		dstype t20 = t2+t15;
		dstype t21 = t3+t14+1.0;
		dstype t22 = 1.0/t20;
		dstype t24 = param2*t21;
		dstype t23 = t22*t22;
		dstype t25 = param1*t22;
		dstype t26 = -t24;
		dstype t27 = param1*t23*udg3*udg4;
		dstype t28 = param1*t23*udg3*udg5;
		dstype t29 = param1*t2*t23;
		dstype t30 = param1*t3*t23;
		dstype t31 = param1*t23*udg4*udg6;
		dstype t32 = param1*t23*udg5*udg6;
		dstype t33 = -t25;
		dstype t34 = -t27;
		dstype t35 = -t28;
		dstype t36 = -t31;
		dstype t37 = -t32;
		dstype t42 = t11+t24+t25+t30;
		dstype t43 = t10+t26+t29+t33;
		dstype t38 = t16+t34;
		dstype t39 = t17+t35;
		dstype t40 = t18+t36;
		dstype t41 = t19+t37;
		f_udg[0*ng+i] = 0.0;
		f_udg[1*ng+i] = 0.0;
		f_udg[2*ng+i] = 0.0;
		f_udg[3*ng+i] = 0.0;
		f_udg[4*ng+i] = 0.0;
		f_udg[5*ng+i] = 0.0;
		f_udg[6*ng+i] = 0.0;
		f_udg[7*ng+i] = 0.0;
		f_udg[8*ng+i] = param1+param2*t7+param1*t7*t23;
		f_udg[9*ng+i] = t41;
		f_udg[10*ng+i] = t40;
		f_udg[11*ng+i] = t43;
		f_udg[12*ng+i] = t41;
		f_udg[13*ng+i] = param1+param2*t6+param1*t6*t23;
		f_udg[14*ng+i] = t42;
		f_udg[15*ng+i] = t39;
		f_udg[16*ng+i] = t40;
		f_udg[17*ng+i] = t42;
		f_udg[18*ng+i] = param1+param2*t5+param1*t5*t23;
		f_udg[19*ng+i] = t38;
		f_udg[20*ng+i] = t43;
		f_udg[21*ng+i] = t39;
		f_udg[22*ng+i] = t38;
		f_udg[23*ng+i] = param1+param2*t4+param1*t4*t23;
		}
	});
}

