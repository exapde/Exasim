void HdgFbouonly11(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbouonly11", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param5 = param[4];
		dstype param6 = param[5];
		dstype param7 = param[6];
		dstype param8 = param[7];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		dstype uhg4 = uhg[3*ng+i];
		f[0*ng+i] = param5-uhg1;
		f[1*ng+i] = param6-uhg2;
		f[2*ng+i] = param7-uhg3;
		f[3*ng+i] = param8-uhg4;
	});
}

void HdgFbouonly12(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbouonly12", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg3 = udg[2*ng+i];
		dstype udg4 = udg[3*ng+i];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		dstype uhg4 = uhg[3*ng+i];
		f[0*ng+i] = udg1-uhg1;
		f[1*ng+i] = udg2-uhg2;
		f[2*ng+i] = udg3-uhg3;
		f[3*ng+i] = udg4-uhg4;
	});
}

void HdgFbouonly13(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbouonly13", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype param9 = param[8];
		dstype param10 = param[9];
		dstype param11 = param[10];
		dstype udg1 = udg[0*ng+i];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		dstype uhg4 = uhg[3*ng+i];
		f[0*ng+i] = udg1-uhg1;
		f[1*ng+i] = -uhg2;
		f[2*ng+i] = -uhg3;
		f[3*ng+i] = -uhg4+(param9*param11*uhg1)/param10;
	});
}

void HdgFbouonly14(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbouonly14", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype param2 = param[1];
		dstype param3 = param[2];
		dstype param4 = param[3];
		dstype param10 = param[9];
		dstype tau1 = tau[0];
		dstype udg1 = udg[0*ng+i];
		dstype udg4 = udg[3*ng+i];
		dstype udg5 = udg[4*ng+i];
		dstype udg6 = udg[5*ng+i];
		dstype udg7 = udg[6*ng+i];
		dstype udg8 = udg[7*ng+i];
		dstype udg9 = udg[8*ng+i];
		dstype udg10 = udg[9*ng+i];
		dstype udg11 = udg[10*ng+i];
		dstype udg12 = udg[11*ng+i];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		dstype uhg4 = uhg[3*ng+i];
		dstype odg1 = odg[0*ng+i];
		dstype nlg1 = nlg[0*ng+i];
		dstype nlg2 = nlg[1*ng+i];
		dstype t2 = param4*param4;
		dstype t3 = uhg2*uhg2;
		dstype t4 = uhg3*uhg3;
		dstype t5 = 1.0/3.141592653589793;
		dstype t6 = param1-1.0;
		dstype t7 = 1.0/param2;
		dstype t8 = 1.0/param3;
		dstype t9 = uhg1*1.0E+3;
		dstype t11 = param10+5.52E+2/5.0;
		dstype t12 = uhg1-1.0/1.0E+2;
		dstype t10 = 1.0/t6;
		dstype t13 = t9-1.0E+1;
		dstype t14 = atan(t13);
		dstype t15 = t5*t14;
		dstype t16 = t15+1.0/2.0;
		dstype t17 = t12*t16;
		dstype t18 = t17*1.0E+3;
		dstype t19 = t17+3.183097800805168E-4;
		dstype t21 = t17+1.031830978008052E-2;
		dstype t20 = t19*t19;
		dstype t22 = t18+3.183097800805168E-1;
		dstype t24 = 1.0/t21;
		dstype t23 = atan(t22);
		dstype t25 = t24*t24;
		dstype t26 = t24*uhg4;
		dstype t27 = t20*1.0E+6;
		dstype t28 = t27+1.0;
		dstype t29 = t5*t23;
		dstype t30 = (t3*t25)/2.0;
		dstype t31 = (t4*t25)/2.0;
		dstype t32 = 1.0/t28;
		dstype t34 = t30+t31;
		dstype t33 = t5*t22*t32;
		dstype t35 = t21*t34;
		dstype t38 = t29+t33+1.0/2.0;
		dstype t40 = t6*(t35-uhg4)*-1.0E+3;
		dstype t41 = t40-1.0;
		dstype t46 = t24*t38*udg5*uhg2;
		dstype t47 = t24*t38*udg5*uhg3;
		dstype t48 = t24*t38*udg9*uhg2;
		dstype t49 = t24*t38*udg9*uhg3;
		dstype t42 = atan(t41);
		dstype t44 = t5*t42;
		dstype t45 = t44+1.0/2.0;
		dstype t61 = t45*(t6*(t35-uhg4)+1.0/1.0E+3)*-1.0E+3;
		dstype t63 = pow(t45*(t6*(t35-uhg4)+1.0/1.0E+3)-3.183097800805168E-4,2.0);
		dstype t71 = -t24*(t45*(t6*(t35-uhg4)+1.0/1.0E+3)-1.318309780080517E-3);
		dstype t76 = -1.0/(param1*param10*t2*t24*(t45*(t6*(t35-uhg4)+1.0/1.0E+3)-1.318309780080517E-3)-5.52E+2/5.0);
		dstype t65 = t61+3.183097800805168E-1;
		dstype t67 = t63*1.0E+6;
		dstype t72 = param1*t2*t71;
		dstype t77 = t26+t71;
		dstype t66 = atan(t65);
		dstype t68 = t67+1.0;
		dstype t74 = pow(t72,3.0/2.0);
		dstype t69 = t5*t66;
		dstype t70 = 1.0/t68;
		dstype t79 = -t5*t70*(t45*(t6*(t35-uhg4)+1.0/1.0E+3)*1.0E+3-3.183097800805168E-1);
		dstype t80 = t69+t79+1.0/2.0;
		f[0*ng+i] = udg1-uhg1;
		f[1*ng+i] = -uhg2;
		f[2*ng+i] = -uhg3;
		f[3*ng+i] = nlg1*(odg1*udg8+t77*uhg2+(t7*t11*t24*t74*uhg3*(t24*(t47-udg7)+t24*(t48-udg10)))/(param1*param10*t2*t24*(t45*(t6*(t35-uhg4)+1.0/1.0E+3)-1.318309780080517E-3)-5.52E+2/5.0)+(t7*t11*t24*t74*uhg2*(t24*(t46-udg6)*2.0-t24*(t49-udg11))*(2.0/3.0))/(param1*param10*t2*t24*(t45*(t6*(t35-uhg4)+1.0/1.0E+3)-1.318309780080517E-3)-5.52E+2/5.0)+param1*t7*t8*t10*t11*t25*t74*t76*(t38*udg5*(t45*(t6*(t35-uhg4)+1.0/1.0E+3)-1.318309780080517E-3)+t6*t21*t80*(udg8+t21*(t25*uhg2*(t46-udg6)+t25*uhg3*(t47-udg7))-t34*t38*udg5)))+nlg2*(odg1*udg12+t77*uhg3+(t7*t11*t24*t74*uhg2*(t24*(t47-udg7)+t24*(t48-udg10)))/(param1*param10*t2*t24*(t45*(t6*(t35-uhg4)+1.0/1.0E+3)-1.318309780080517E-3)-5.52E+2/5.0)-(t7*t11*t24*t74*uhg3*(t24*(t46-udg6)-t24*(t49-udg11)*2.0)*(2.0/3.0))/(param1*param10*t2*t24*(t45*(t6*(t35-uhg4)+1.0/1.0E+3)-1.318309780080517E-3)-5.52E+2/5.0)+param1*t7*t8*t10*t11*t25*t74*t76*(t38*udg9*(t45*(t6*(t35-uhg4)+1.0/1.0E+3)-1.318309780080517E-3)+t6*t21*t80*(udg12+t21*(t25*uhg2*(t48-udg10)+t25*uhg3*(t49-udg11))-t34*t38*udg9)))+tau1*(udg4-uhg4);
	});
}

void HdgFbouonly15(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbouonly15", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype udg1 = udg[0*ng+i];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		dstype uhg4 = uhg[3*ng+i];
		f[0*ng+i] = udg1-uhg1;
		f[1*ng+i] = -uhg2;
		f[2*ng+i] = -uhg3;
		f[3*ng+i] = uhg4/uhg1;
	});
}

void HdgFbouonly1(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	if (ib == 1)
		HdgFbouonly11(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 2)
		HdgFbouonly12(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 3)
		HdgFbouonly13(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 4)
		HdgFbouonly14(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 5)
		HdgFbouonly15(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

