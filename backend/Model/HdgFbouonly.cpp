void HdgFbouonly1(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbouonly1", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype tau1 = tau[0];
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg3 = udg[2*ng+i];
		dstype udg4 = udg[3*ng+i];
		dstype udg5 = udg[4*ng+i];
		dstype udg6 = udg[5*ng+i];
		dstype udg7 = udg[6*ng+i];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		dstype uhg4 = uhg[3*ng+i];
		dstype uhg5 = uhg[4*ng+i];
		dstype uhg6 = uhg[5*ng+i];
		dstype uhg7 = uhg[6*ng+i];
		dstype nlg1 = nlg[0*ng+i];
		dstype nlg2 = nlg[1*ng+i];
		dstype t2 = udg5*udg6;
		dstype t3 = udg2*udg2;
		dstype t4 = udg3*udg3;
		dstype t5 = udg5*udg5;
		dstype t6 = udg6*udg6;
		dstype t7 = udg7*udg7;
		dstype t8 = param1-1.0;
		dstype t9 = -udg4;
		dstype t10 = 1.0/udg1;
		dstype t11 = t10*t10;
		dstype t12 = t10*udg4;
		dstype t13 = t10*udg2*udg3;
		dstype t14 = t10*udg2*udg5;
		dstype t15 = t10*udg2*udg6;
		dstype t16 = t10*udg3*udg5;
		dstype t17 = t10*udg3*udg6;
		dstype t18 = t5/2.0;
		dstype t19 = t6/2.0;
		dstype t20 = t7/2.0;
		dstype t21 = -t13;
		dstype t22 = -t16;
		dstype t23 = (t3*t11)/2.0;
		dstype t24 = (t4*t11)/2.0;
		dstype t26 = t14+t17;
		dstype t27 = t18+t19;
		dstype t25 = t2+t21;
		dstype t28 = t15+t22;
		dstype t29 = t10*t27;
		dstype t30 = t23+t24;
		dstype t31 = t30*udg1;
		dstype t32 = t9+t20+t27+t31;
		dstype t33 = t8*t32;
		dstype t34 = -t33;
		dstype t36 = t10*t34;
		dstype t37 = t12+t29+t36;
		f[0*ng+i] = nlg1*udg2+nlg2*udg3+tau1*(udg1-uhg1);
		f[1*ng+i] = -nlg1*(t18-t19+t33-t3*t10)-nlg2*t25+tau1*(udg2-uhg2);
		f[2*ng+i] = nlg2*(t18-t19+t34+t4*t10)-nlg1*t25+tau1*(udg3-uhg3);
		f[3*ng+i] = nlg1*(-t26*udg5+t37*udg2+udg5*udg7)+nlg2*(-t26*udg6+t37*udg3+udg6*udg7)+tau1*(udg4-uhg4);
		f[4*ng+i] = -nlg2*t28+nlg1*udg7+tau1*(udg5-uhg5);
		f[5*ng+i] = nlg1*t28+nlg2*udg7+tau1*(udg6-uhg6);
		f[6*ng+i] = nlg1*udg5+nlg2*udg6+tau1*(udg7-uhg7);
	});
}

void HdgFbouonly(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	if (ib == 1)
		HdgFbouonly1(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

