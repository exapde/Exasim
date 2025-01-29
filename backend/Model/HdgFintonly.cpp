void HdgFintonly1(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fintonly1", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype param2 = param[1];
		dstype param3 = param[2];
		dstype param4 = param[3];
		dstype param10 = param[9];
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
		dstype uhg4 = uhg[3*ng+i];
		dstype odg1 = odg[0*ng+i];
		dstype nlg1 = nlg[0*ng+i];
		dstype nlg2 = nlg[1*ng+i];
		dstype t2 = 1.0/3.141592653589793;
		dstype t3 = udg1*1.0E3;
		dstype t4 = t3-1.0E1;
		dstype t5 = atan(t4);
		dstype t6 = t2*t5;
		dstype t7 = t6+1.0/2.0;
		dstype t8 = udg1-1.0/1.0E2;
		dstype t9 = t7*t8;
		dstype t10 = t9+1.031830978008052E-2;
		dstype t11 = 1.0/t10;
		dstype t12 = 1.0/(t10*t10);
		dstype t13 = udg2*udg2;
		dstype t14 = t12*t13*(1.0/2.0);
		dstype t15 = udg3*udg3;
		dstype t16 = t12*t15*(1.0/2.0);
		dstype t17 = t14+t16;
		dstype t31 = t10*t17;
		dstype t18 = -t31+udg4;
		dstype t19 = param1-1.0;
		dstype t20 = t7*t8*1.0E3;
		dstype t21 = t20+3.183097800805168E-1;
		dstype t22 = t9+3.183097800805168E-4;
		dstype t23 = atan(t21);
		dstype t24 = t2*t23;
		dstype t25 = t22*t22;
		dstype t26 = t25*1.0E6;
		dstype t27 = t26+1.0;
		dstype t28 = 1.0/t27;
		dstype t29 = t2*t21*t28;
		dstype t30 = t24+t29+1.0/2.0;
		dstype t32 = t18*t19*1.0E3;
		dstype t33 = t32-1.0;
		dstype t34 = atan(t33);
		dstype t35 = t2*t34;
		dstype t36 = t35+1.0/2.0;
		dstype t37 = t18*t19;
		dstype t38 = t37-1.0/1.0E3;
		dstype t39 = t36*t38;
		dstype t40 = t39+1.318309780080517E-3;
		dstype t41 = param4*param4;
		dstype t42 = 1.0/param2;
		dstype t43 = param1*param10*t11*t40*t41;
		dstype t44 = t43+5.52E2/5.0;
		dstype t45 = 1.0/t44;
		dstype t46 = param10+5.52E2/5.0;
		dstype t47 = param1*t11*t40*t41;
		dstype t48 = pow(t47,3.0/2.0);
		dstype t49 = t39+3.183097800805168E-4;
		dstype t50 = t36*t38*1.0E3;
		dstype t51 = t50+3.183097800805168E-1;
		dstype t62 = t11*t30*udg2*udg5;
		dstype t52 = -t62+udg6;
		dstype t57 = t11*t30*udg3*udg5;
		dstype t53 = -t57+udg7;
		dstype t54 = t11*udg4;
		dstype t55 = t11*t40;
		dstype t56 = t54+t55;
		dstype t58 = t11*t53;
		dstype t73 = t11*t30*udg2*udg9;
		dstype t59 = -t73+udg10;
		dstype t60 = t11*t59;
		dstype t61 = t58+t60;
		dstype t74 = t11*t30*udg3*udg9;
		dstype t63 = -t74+udg11;
		dstype t64 = 1.0/param3;
		dstype t65 = atan(t51);
		dstype t66 = t2*t65;
		dstype t67 = t49*t49;
		dstype t68 = t67*1.0E6;
		dstype t69 = t68+1.0;
		dstype t70 = 1.0/t69;
		dstype t71 = t2*t51*t70;
		dstype t72 = t66+t71+1.0/2.0;
		dstype t75 = 1.0/t19;
		f[0*ng+i] = tau1*(udg4-uhg4)+nlg1*(odg1*udg8+t56*udg2+t11*t42*t45*t46*t48*udg2*(t11*t52*2.0-t11*t63)*(2.0/3.0)+t11*t42*t45*t46*t48*t61*udg3-param1*t12*t42*t45*t46*t48*t64*t75*(t30*t40*udg5+t10*t19*t72*(-udg8+t10*(t12*t52*udg2+t12*t53*udg3)+t17*t30*udg5)))+nlg2*(odg1*udg12+t56*udg3-t11*t42*t45*t46*t48*udg3*(t11*t52-t11*t63*2.0)*(2.0/3.0)+t11*t42*t45*t46*t48*t61*udg2-param1*t12*t42*t45*t46*t48*t64*t75*(t30*t40*udg9+t10*t19*t72*(-udg12+t10*(t12*t59*udg2+t12*t63*udg3)+t17*t30*udg9)));
	});
}

void HdgFintonly(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	if (ib == 1)
		HdgFintonly1(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

