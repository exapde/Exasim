void KokkosUbou1(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Ubou1", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype udg1 = udg[0*ng+i];
		dstype udg4 = udg[3*ng+i];
		f[0*ng+i] = udg1;
		f[1*ng+i] = 0.0;
		f[2*ng+i] = 0.0;
		f[3*ng+i] = udg4;
	});
}

void KokkosUbou2(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Ubou2", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype param5 = param[4];
		dstype param6 = param[5];
		dstype param7 = param[6];
		dstype param8 = param[7];
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg3 = udg[2*ng+i];
		dstype udg4 = udg[3*ng+i];
		dstype nlg1 = nlg[0*ng+i];
		dstype nlg2 = nlg[1*ng+i];
		dstype t2 = 1.0/udg1;
		dstype t3 = 1.0/(udg1*udg1);
		dstype t4 = udg2*udg2;
		dstype t5 = t3*t4;
		dstype t6 = udg3*udg3;
		dstype t7 = t3*t6;
		dstype t8 = t5+t7;
		dstype t10 = t8*udg1*(1.0/2.0);
		dstype t9 = -t10+udg4;
		dstype t11 = param1-1.0;
		dstype t12 = nlg1*udg2;
		dstype t13 = nlg2*udg3;
		dstype t14 = t12+t13;
		dstype t15 = param1*t2*t9*t11;
		dstype t16 = sqrt(t15);
		dstype t17 = 1.0/param1;
		dstype t18 = t2*t14*1.0E3;
		dstype t19 = 1.0/t9;
		dstype t20 = t16*1.0E3;
		dstype t21 = t2*udg4;
		dstype t22 = t2*t9*t11;
		dstype t23 = t2*t14;
		dstype t24 = 1.0/t11;
		dstype t25 = tanh(t18);
		dstype t26 = t18-t20;
		dstype t27 = tanh(t26);
		dstype t28 = nlg1*udg3;
		dstype t32 = nlg2*udg2;
		dstype t29 = t28-t32;
		dstype t30 = t18+t20;
		dstype t31 = tanh(t30);
		dstype t33 = t16*t24;
		dstype t34 = t17*t19*t27*t29*(1.0/2.0);
		dstype t35 = t17*t19*t29*t31*(1.0/2.0);
		dstype t45 = t17*t19*t25*t29;
		dstype t36 = t34+t35-t45;
		dstype t37 = t23+t33;
		dstype t38 = t17*t19*t27*t37*udg1*(1.0/2.0);
		dstype t39 = t23-t33;
		dstype t40 = t17*t19*t31*t39*udg1*(1.0/2.0);
		dstype t42 = t14*t17*t19*t25;
		dstype t41 = t38+t40-t42;
		dstype t43 = t16-t23;
		dstype t44 = t16+t23;
		dstype t46 = nlg2*t36;
		dstype t47 = param8-udg4;
		dstype t48 = t29*t29;
		dstype t49 = t2*t14*t17*t19*t25*t29;
		dstype t50 = param6-udg2;
		dstype t51 = t2*t17*t19*t25*t48;
		dstype t79 = t2*t17*t19*t27*t48*(1.0/2.0);
		dstype t80 = t2*t17*t19*t31*t48*(1.0/2.0);
		dstype t52 = t25+t51-t79-t80;
		dstype t53 = nlg2*t52;
		dstype t54 = t17*t19*t27*t29*t43*(1.0/2.0);
		dstype t81 = t17*t19*t29*t31*t44*(1.0/2.0);
		dstype t55 = t49+t54-t81;
		dstype t56 = t53-nlg1*t55;
		dstype t57 = t14*t14;
		dstype t58 = t2*t17*t19*t25*t57;
		dstype t59 = t17*t19*t27*t37*t43*udg1*(1.0/2.0);
		dstype t82 = t17*t19*t31*t39*t44*udg1*(1.0/2.0);
		dstype t60 = t58+t59-t82;
		dstype t61 = nlg1*t60;
		dstype t62 = t17*t19*t27*t29*t37*(1.0/2.0);
		dstype t63 = t17*t19*t29*t31*t39*(1.0/2.0);
		dstype t64 = -t49+t62+t63;
		dstype t65 = nlg2*t64;
		dstype t66 = t61+t65;
		dstype t67 = param7-udg3;
		dstype t68 = param5-udg1;
		dstype t69 = t2*udg4*2.0;
		dstype t70 = t2*t9*t11*2.0;
		dstype t74 = param1*t2*t9*4.0;
		dstype t71 = t69+t70-t74;
		dstype t72 = -t16+t23;
		dstype t77 = t16*t24*t44;
		dstype t73 = t21+t22-t77;
		dstype t75 = t16*t24*t72;
		dstype t76 = t21+t22+t75;
		dstype t78 = nlg1*t36;
		dstype t83 = t17*t19*t27*t29*t72*(1.0/2.0);
		dstype t84 = -t49+t81+t83;
		dstype t85 = nlg1*t52-nlg2*t84;
		dstype t86 = t17*t19*t27*t37*t72*udg1*(1.0/2.0);
		dstype t87 = -t58+t82+t86;
		dstype t88 = nlg2*t87;
		dstype t89 = nlg1*t64;
		dstype t90 = t88+t89;
		dstype t91 = t17*t19*t27*t72*t76*udg1*(1.0/2.0);
		dstype t92 = t17*t19*t31*t44*t73*udg1*(1.0/2.0);
		dstype t93 = t91+t92-t14*t17*t19*t25*t71*(1.0/2.0);
		dstype t94 = t2*t25*t29;
		dstype t95 = t17*t19*t25*t29*t71*(1.0/2.0);
		dstype t96 = t94+t95-t17*t19*t27*t29*t76*(1.0/2.0)-t17*t19*t29*t31*t73*(1.0/2.0);
		dstype t97 = t2*t14*t16;
		dstype t98 = t3*t57*(1.0/2.0);
		dstype t99 = t3*t48*(1.0/2.0);
		dstype t100 = t98+t99;
		dstype t101 = t21+t22-t97;
		dstype t102 = t21+t22+t97;
		dstype t103 = t17*t19*t25*t29*t100;
		dstype t104 = t94+t103-t17*t19*t27*t29*t101*(1.0/2.0)-t17*t19*t29*t31*t102*(1.0/2.0);
		dstype t105 = t17*t19*t27*t37*t101*udg1*(1.0/2.0);
		dstype t106 = t17*t19*t31*t39*t102*udg1*(1.0/2.0);
		dstype t107 = t105+t106-t14*t17*t19*t25*t100;
		f[0*ng+i] = param5*(1.0/2.0)+udg1*(1.0/2.0)-t68*(t17*t19*t27*udg1*(t21+t22-t16*t24*t43)*(1.0/2.0)-t17*t19*t25*t71*udg1*(1.0/2.0)+t17*t19*t31*t73*udg1*(1.0/2.0))*(1.0/2.0)-t47*(-t17*t19*t25*udg1+t17*t19*t27*udg1*(1.0/2.0)+t17*t19*t31*udg1*(1.0/2.0))*(1.0/2.0)-t50*(t46-nlg1*t41)*(1.0/2.0)+t67*(t78+nlg2*t41)*(1.0/2.0);
		f[1*ng+i] = param6*(1.0/2.0)+udg2*(1.0/2.0)+t47*(t46+nlg1*(t42+t17*t19*t27*t43*udg1*(1.0/2.0)-t17*t19*t31*t44*udg1*(1.0/2.0)))*(1.0/2.0)-t50*(nlg2*t56+nlg1*t66)*(1.0/2.0)+t67*(nlg1*t56-nlg2*t66)*(1.0/2.0)-t68*(nlg1*t93+nlg2*t96)*(1.0/2.0);
		f[2*ng+i] = param7*(1.0/2.0)+udg3*(1.0/2.0)+t50*(nlg2*t85+nlg1*t90)*(1.0/2.0)-t67*(nlg1*t85-nlg2*t90)*(1.0/2.0)-t68*(nlg2*t93-nlg1*t96)*(1.0/2.0)-t47*(t78+nlg2*(-t42+t17*t19*t31*t44*udg1*(1.0/2.0)+t17*t19*t27*t72*udg1*(1.0/2.0)))*(1.0/2.0);
		f[3*ng+i] = param8*(1.0/2.0)+udg4*(1.0/2.0)-t47*(-t17*t19*t25*t100*udg1+t17*t19*t27*t101*udg1*(1.0/2.0)+t17*t19*t31*t102*udg1*(1.0/2.0))*(1.0/2.0)+t68*(t3*t25*t48+t17*t19*t25*t71*t100*udg1*(1.0/2.0)-t17*t19*t27*t76*t101*udg1*(1.0/2.0)-t17*t19*t31*t73*t102*udg1*(1.0/2.0))*(1.0/2.0)+t50*(nlg2*t104+nlg1*t107)*(1.0/2.0)-t67*(nlg1*t104-nlg2*t107)*(1.0/2.0);
	});
}

void KokkosUbou(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	if (ib == 1)
		KokkosUbou1(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 2)
		KokkosUbou2(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

