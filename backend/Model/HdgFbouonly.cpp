void HdgFbouonly1(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbouonly1", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype param2 = param[1];
		dstype param3 = param[2];
		dstype tau1 = tau[0];
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg3 = udg[2*ng+i];
		dstype udg4 = udg[3*ng+i];
		dstype udg5 = udg[4*ng+i];
		dstype udg6 = udg[5*ng+i];
		dstype udg7 = udg[6*ng+i];
		dstype udg8 = udg[7*ng+i];
		dstype udg9 = udg[8*ng+i];
		dstype udg10 = udg[9*ng+i];
		dstype udg11 = udg[10*ng+i];
		dstype udg12 = udg[11*ng+i];
		dstype udg13 = udg[12*ng+i];
		dstype udg14 = udg[13*ng+i];
		dstype udg15 = udg[14*ng+i];
		dstype udg16 = udg[15*ng+i];
		dstype udg17 = udg[16*ng+i];
		dstype udg18 = udg[17*ng+i];
		dstype udg19 = udg[18*ng+i];
		dstype udg20 = udg[19*ng+i];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		dstype uhg4 = uhg[3*ng+i];
		dstype uhg5 = uhg[4*ng+i];
		dstype nlg1 = nlg[0*ng+i];
		dstype nlg2 = nlg[1*ng+i];
		dstype nlg3 = nlg[2*ng+i];
		dstype t2 = udg2*udg2;
		dstype t3 = udg3*udg3;
		dstype t4 = udg4*udg4;
		dstype t5 = param1-1.0;
		dstype t6 = 1.0/param2;
		dstype t7 = 1.0/param3;
		dstype t8 = 1.0/udg1;
		dstype t9 = t8*t8;
		dstype t10 = t8*udg5;
		dstype t11 = t8*udg2*udg3;
		dstype t12 = t8*udg2*udg4;
		dstype t13 = t8*udg3*udg4;
		dstype t14 = t8*udg2*udg6;
		dstype t15 = t8*udg3*udg6;
		dstype t16 = t8*udg4*udg6;
		dstype t17 = t8*udg2*udg11;
		dstype t18 = t8*udg3*udg11;
		dstype t19 = t8*udg4*udg11;
		dstype t20 = t8*udg2*udg16;
		dstype t21 = t8*udg3*udg16;
		dstype t22 = t8*udg4*udg16;
		dstype t23 = 1.0/t5;
		dstype t33 = (t2*t9)/2.0;
		dstype t34 = (t3*t9)/2.0;
		dstype t35 = (t4*t9)/2.0;
		dstype t54 = t8*(t14-udg7)*-2.0;
		dstype t55 = t8*(t18-udg13)*-2.0;
		dstype t56 = t8*(t22-udg19)*-2.0;
		dstype t66 = -t6*(t8*(t15-udg8)+t8*(t17-udg12));
		dstype t67 = -t6*(t8*(t16-udg9)+t8*(t20-udg17));
		dstype t68 = -t6*(t8*(t19-udg14)+t8*(t21-udg18));
		dstype t63 = t33+t34+t35;
		dstype t71 = t11+t66;
		dstype t72 = t12+t67;
		dstype t73 = t13+t68;
		dstype t64 = t63*udg1;
		dstype t70 = -t5*(t64-udg5);
		dstype t74 = t8*t70;
		dstype t78 = t10+t74;
		f[0*ng+i] = nlg1*udg2+nlg2*udg3+nlg3*udg4+tau1*(udg1-uhg1);
		f[1*ng+i] = nlg2*t71+nlg3*t72+tau1*(udg2-uhg2)+nlg1*(t70+t2*t8+t6*(t54+t8*(t18-udg13)+t8*(t22-udg19))*(2.0/3.0));
		f[2*ng+i] = nlg1*t71+nlg3*t73+tau1*(udg3-uhg3)+nlg2*(t70+t3*t8+t6*(t55+t8*(t14-udg7)+t8*(t22-udg19))*(2.0/3.0));
		f[3*ng+i] = nlg1*t72+nlg2*t73+tau1*(udg4-uhg4)+nlg3*(t70+t4*t8+t6*(t56+t8*(t14-udg7)+t8*(t18-udg13))*(2.0/3.0));
		f[4*ng+i] = tau1*(udg5-uhg5)+nlg1*(t78*udg2+t8*t66*udg3+t8*t67*udg4+t6*t8*udg2*(t54+t8*(t18-udg13)+t8*(t22-udg19))*(2.0/3.0)+param1*t6*t7*t9*t23*(t5*udg1*(udg10+udg1*(t9*udg2*(t14-udg7)+t9*udg3*(t15-udg8)+t9*udg4*(t16-udg9))-t63*udg6)+t5*udg6*(t64-udg5)))+nlg2*(t78*udg3+t8*t66*udg2+t8*t68*udg4+t6*t8*udg3*(t55+t8*(t14-udg7)+t8*(t22-udg19))*(2.0/3.0)+param1*t6*t7*t9*t23*(t5*udg1*(udg15+udg1*(t9*udg2*(t17-udg12)+t9*udg3*(t18-udg13)+t9*udg4*(t19-udg14))-t63*udg11)+t5*udg11*(t64-udg5)))+nlg3*(t78*udg4+t8*t67*udg2+t8*t68*udg3+t6*t8*udg4*(t56+t8*(t14-udg7)+t8*(t18-udg13))*(2.0/3.0)+param1*t6*t7*t9*t23*(t5*udg1*(udg20+udg1*(t9*udg2*(t20-udg17)+t9*udg3*(t21-udg18)+t9*udg4*(t22-udg19))-t63*udg16)+t5*udg16*(t64-udg5)));
	});
}

void HdgFbouonly(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	if (ib == 1)
		HdgFbouonly1(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

