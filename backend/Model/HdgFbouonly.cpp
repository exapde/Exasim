void HdgFbouonly1(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbouonly1", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype param2 = param[1];
		dstype param3 = param[2];
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
		dstype nlg1 = nlg[0*ng+i];
		dstype nlg2 = nlg[1*ng+i];
		dstype t2 = 1.0/uhg1;
		dstype t3 = 1.0/(uhg1*uhg1);
		dstype t4 = 1.0/param2;
		dstype t5 = param1-1.0;
		dstype t6 = uhg2*uhg2;
		dstype t7 = t3*t6*(1.0/2.0);
		dstype t8 = uhg3*uhg3;
		dstype t9 = t3*t8*(1.0/2.0);
		dstype t10 = t7+t9;
		dstype t23 = t2*udg5*uhg2;
		dstype t11 = -t23+udg6;
		dstype t18 = t2*udg5*uhg3;
		dstype t12 = -t18+udg7;
		dstype t15 = t10*uhg1;
		dstype t13 = -t15+uhg4;
		dstype t14 = t2*uhg4;
		dstype t16 = t2*t5*t13;
		dstype t17 = t14+t16;
		dstype t19 = t2*t12;
		dstype t26 = t2*udg9*uhg2;
		dstype t20 = -t26+udg10;
		dstype t21 = t2*t20;
		dstype t22 = t19+t21;
		dstype t27 = t2*udg9*uhg3;
		dstype t24 = -t27+udg11;
		dstype t25 = 1.0/param3;
		dstype t28 = 1.0/t5;
		f[0*ng+i] = udg1-uhg1;
		f[1*ng+i] = -uhg2;
		f[2*ng+i] = -uhg3;
		f[3*ng+i] = tau1*(udg4-uhg4)+nlg1*(t17*uhg2+t2*t4*t22*uhg3+t2*t4*uhg2*(t2*t11*2.0-t2*t24)*(2.0/3.0)-param1*t3*t4*t25*t28*(t5*uhg1*(-udg8+t10*udg5+uhg1*(t3*t11*uhg2+t3*t12*uhg3))+t5*t13*udg5))+nlg2*(t17*uhg3+t2*t4*t22*uhg2-t2*t4*uhg3*(t2*t11-t2*t24*2.0)*(2.0/3.0)-param1*t3*t4*t25*t28*(t5*uhg1*(-udg12+t10*udg9+uhg1*(t3*t20*uhg2+t3*t24*uhg3))+t5*t13*udg9));
	});
}

void HdgFbouonly2(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbouonly2", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype param5 = param[4];
		dstype param6 = param[5];
		dstype param7 = param[6];
		dstype param8 = param[7];
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg3 = udg[2*ng+i];
		dstype udg4 = udg[3*ng+i];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		dstype uhg4 = uhg[3*ng+i];
		dstype nlg1 = nlg[0*ng+i];
		dstype nlg2 = nlg[1*ng+i];
		dstype t2 = uhg2*uhg2;
		dstype t3 = uhg3*uhg3;
		dstype t4 = 1.0/param1;
		dstype t5 = 1.0/uhg1;
		dstype t7 = uhg1*uhg4*2.0;
		dstype t6 = t2+t3-t7;
		dstype t8 = 1.0/t6;
		dstype t9 = 1.0/(uhg1*uhg1);
		dstype t10 = param1-1.0;
		dstype t11 = sqrt(2.0);
		dstype t12 = nlg1*uhg2*2.0;
		dstype t13 = nlg2*uhg3*2.0;
		dstype t17 = param1*t6*t9*t10;
		dstype t14 = sqrt(-t17);
		dstype t15 = t11*t14*uhg1;
		dstype t16 = t12+t13+t15;
		dstype t18 = uhg1*uhg1;
		dstype t19 = t6*t9*t10*(1.0/2.0);
		dstype t20 = 1.0/t10;
		dstype t21 = t12+t13-t15;
		dstype t22 = t5*t16*5.0E1;
		dstype t23 = tanh(t22);
		dstype t24 = nlg1*nlg1;
		dstype t25 = nlg2*nlg2;
		dstype t26 = t5*t21*5.0E1;
		dstype t27 = tanh(t26);
		dstype t28 = t24*uhg2;
		dstype t29 = t25*uhg2;
		dstype t30 = nlg1*t11*t14*uhg1*(1.0/2.0);
		dstype t31 = nlg1*uhg2;
		dstype t32 = nlg2*uhg3;
		dstype t33 = t31+t32;
		dstype t34 = t5*t33*1.0E2;
		dstype t35 = tanh(t34);
		dstype t36 = t24*uhg3;
		dstype t37 = t25*uhg3;
		dstype t38 = nlg2*t11*t14*uhg1*(1.0/2.0);
		dstype t39 = t24+t25;
		dstype t40 = param8-udg4;
		dstype t41 = t4*t8*t35*t39*uhg1*uhg2*2.0;
		dstype t42 = param5-udg1;
		dstype t43 = t28+t29+t30;
		dstype t44 = t5*uhg4;
		dstype t45 = t5*t11*t14*t16*t20*(1.0/4.0);
		dstype t46 = t28+t29-t30;
		dstype t47 = t5*t11*t14*t20*t21*(1.0/4.0);
		dstype t48 = -t19+t44+t47;
		dstype t49 = param1*t2;
		dstype t50 = param1*t3;
		dstype t67 = param1*uhg1*uhg4*2.0;
		dstype t51 = t2+t3+t49+t50-t67;
		dstype t52 = param6-udg2;
		dstype t53 = param1*t24*uhg2;
		dstype t54 = param1*t25*uhg2;
		dstype t55 = -t28-t29+t30+t53+t54;
		dstype t56 = param7-udg3;
		dstype t57 = t39*t39;
		dstype t58 = param1*t24*uhg3;
		dstype t59 = param1*t25*uhg3;
		dstype t60 = -t36-t37+t38+t58+t59;
		dstype t61 = t4*t8*t35*t39*uhg1*uhg3*2.0;
		dstype t62 = nlg1*uhg3;
		dstype t73 = nlg2*uhg2;
		dstype t63 = t62-t73;
		dstype t64 = t36+t37+t38;
		dstype t65 = t19-t44+t45;
		dstype t66 = t36+t37-t38;
		dstype t68 = t36+t37+t38-t58-t59;
		dstype t69 = nlg1*nlg2*t35;
		dstype t70 = t4*t8*t35*t57*uhg2*uhg3*2.0;
		dstype t71 = t28+t29+t30-t53-t54;
		dstype t72 = t5*t11*t14*t33*(1.0/2.0);
		dstype t74 = -t19+t44+t72;
		dstype t75 = t19-t44+t72;
		dstype t76 = t2+t3;
		dstype t77 = nlg2*t5*t35*t63;
		dstype t78 = nlg1*t5*t35*t63;
		f[0*ng+i] = param5*(1.0/2.0)+udg1*(1.0/2.0)-uhg1+t52*(t41-t4*t8*t20*t27*t55*uhg1+t4*t8*t20*t23*uhg1*(t28+t29+t30-param1*t24*uhg2-param1*t25*uhg2))*(1.0/2.0)+t56*(t61-t4*t8*t20*t27*t60*uhg1+t4*t8*t20*t23*uhg1*(t36+t37+t38-param1*t24*uhg3-param1*t25*uhg3))*(1.0/2.0)-t42*(t4*t8*t35*t51-t4*t8*t18*t27*t48+t4*t8*t18*t23*(t19+t45-t5*uhg4))*(1.0/2.0)+t4*t8*t18*t40*(t23+t27-t35*2.0)*(1.0/2.0);
		f[1*ng+i] = param6*(1.0/2.0)+udg2*(1.0/2.0)-uhg2+t56*(t69+t70-t4*t8*t20*t27*t46*t60+t4*t8*t20*t23*t43*t68)*(1.0/2.0)-t52*(t25*t35-t2*t4*t8*t35*t57*2.0+t4*t8*t20*t27*t46*t55-t4*t8*t20*t23*t43*t71)*(1.0/2.0)+t40*(-t41+t4*t8*t23*t43*uhg1+t4*t8*t27*t46*uhg1)*(1.0/2.0)-t42*(t77-t4*t8*t27*t46*t48*uhg1+t4*t8*t23*t43*t65*uhg1+t4*t5*t8*t35*t39*t51*uhg2)*(1.0/2.0);
		f[2*ng+i] = param7*(1.0/2.0)+udg3*(1.0/2.0)-uhg3+t52*(t69+t70-t4*t8*t20*t27*t55*t66+t4*t8*t20*t23*t64*t71)*(1.0/2.0)-t56*(t24*t35-t3*t4*t8*t35*t57*2.0+t4*t8*t20*t27*t60*t66-t4*t8*t20*t23*t64*t68)*(1.0/2.0)+t40*(-t61+t4*t8*t23*t64*uhg1+t4*t8*t27*t66*uhg1)*(1.0/2.0)+t42*(t78+t4*t8*t27*t48*t66*uhg1-t4*t8*t23*t64*t65*uhg1-t4*t5*t8*t35*t39*t51*uhg3)*(1.0/2.0);
		f[3*ng+i] = param8*(1.0/2.0)+udg4*(1.0/2.0)-uhg4-t42*(-t9*t35*(t63*t63)+t4*t8*t18*t27*t48*t75+t4*t8*t18*t23*t65*t74+t4*t8*t9*t35*t39*t51*t76*(1.0/2.0))*(1.0/2.0)-t40*(-t4*t8*t18*t23*t74+t4*t8*t18*t27*t75+t4*t8*t35*t39*t76)*(1.0/2.0)+t52*(t77+t4*t5*t8*t35*t57*t76*uhg2+t4*t8*t20*t27*t55*t75*uhg1+t4*t8*t20*t23*t71*t74*uhg1)*(1.0/2.0)+t56*(-t78+t4*t5*t8*t35*t57*t76*uhg3+t4*t8*t20*t27*t60*t75*uhg1+t4*t8*t20*t23*t68*t74*uhg1)*(1.0/2.0);
	});
}

void HdgFbouonly(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	if (ib == 1)
		HdgFbouonly1(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 2)
		HdgFbouonly2(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

