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
		dstype t2 = nlg1*udg2;
		dstype t3 = nlg1*udg3;
		dstype t4 = nlg2*udg2;
		dstype t5 = nlg2*udg3;
		dstype t6 = udg2*udg2;
		dstype t7 = udg3*udg3;
		dstype t8 = param1-1.0;
		dstype t9 = 1.0/param1;
		dstype t10 = -udg1;
		dstype t11 = -udg2;
		dstype t12 = -udg3;
		dstype t13 = -udg4;
		dstype t14 = 1.0/udg1;
		dstype t15 = t14*t14;
		dstype t16 = -t4;
		dstype t17 = param5+t10;
		dstype t18 = param6+t11;
		dstype t19 = param7+t12;
		dstype t20 = param8+t13;
		dstype t21 = t14*udg4;
		dstype t22 = 1.0/t8;
		dstype t24 = t2+t5;
		dstype t23 = t21*2.0;
		dstype t25 = t6*t15;
		dstype t26 = t7*t15;
		dstype t27 = t24*t24;
		dstype t28 = t3+t16;
		dstype t30 = t14*t24;
		dstype t29 = t28*t28;
		dstype t31 = t30*1.0E+3;
		dstype t32 = t25+t26;
		dstype t33 = (t15*t27)/2.0;
		dstype t34 = tanh(t31);
		dstype t35 = (t15*t29)/2.0;
		dstype t36 = (t32*udg1)/2.0;
		dstype t39 = -1.0/(t36-udg4);
		dstype t40 = t14*t28*t34;
		dstype t42 = param1*t14*(t36-udg4)*4.0;
		dstype t43 = -t8*t14*(t36-udg4);
		dstype t44 = t8*t14*(t36-udg4)*-2.0;
		dstype t46 = t33+t35;
		dstype t59 = (t9*t24*t34)/(t36-udg4);
		dstype t61 = (t9*t28*t34)/(t36-udg4);
		dstype t45 = param1*t43;
		dstype t57 = t9*t24*t34*t39;
		dstype t58 = t9*t28*t34*t39;
		dstype t67 = t9*t14*t27*t34*t39;
		dstype t68 = t9*t14*t29*t34*t39;
		dstype t72 = t23+t42+t44;
		dstype t77 = t46*t59;
		dstype t47 = sqrt(t45);
		dstype t70 = t30*t58;
		dstype t76 = t46*t58;
		dstype t88 = (t59*t72)/2.0;
		dstype t89 = t61*t72*(-1.0/2.0);
		dstype t48 = -t47;
		dstype t49 = t47*1.0E+3;
		dstype t51 = t22*t47;
		dstype t53 = t30*t47;
		dstype t54 = t30+t47;
		dstype t50 = -t49;
		dstype t52 = t22*t48;
		dstype t56 = t30+t48;
		dstype t60 = t30+t51;
		dstype t62 = t31+t49;
		dstype t73 = t21+t43+t53;
		dstype t78 = t51*t54;
		dstype t63 = tanh(t62);
		dstype t64 = t31+t50;
		dstype t65 = t30+t52;
		dstype t80 = t51*t56;
		dstype t66 = tanh(t64);
		dstype t81 = (t9*t28*t63*(-1.0/2.0))/(t36-udg4);
		dstype t85 = (t9*t14*t29*t63)/(t36*2.0-udg4*2.0);
		dstype t91 = t21+t43+t80;
		dstype t96 = (t9*t63*t65*udg1*(-1.0/2.0))/(t36-udg4);
		dstype t101 = (t9*t28*t63*t73)/(t36*2.0-udg4*2.0);
		dstype t115 = (t9*t54*t63*udg1*(-t21+t78+t8*t14*(t36-udg4)))/(t36*2.0-udg4*2.0);
		dstype t82 = (t9*t28*t66*(-1.0/2.0))/(t36-udg4);
		dstype t86 = (t9*t14*t29*t66)/(t36*2.0-udg4*2.0);
		dstype t93 = (t9*t56*t66*udg1*(-1.0/2.0))/(t36-udg4);
		dstype t107 = t81*(-t21+t78+t8*t14*(t36-udg4));
		dstype t109 = (t9*t28*t66*t91)/(t36*2.0-udg4*2.0);
		dstype t110 = t73*t96;
		dstype t111 = (t9*t60*t66*udg1*(-t21+t53+t8*t14*(t36-udg4)))/(t36*2.0-udg4*2.0);
		dstype t114 = -nlg2*(t58+(t9*t28*t63)/(t36*2.0-udg4*2.0)+(t9*t28*t66)/(t36*2.0-udg4*2.0));
		dstype t123 = -nlg1*(t70+(t9*t28*t54*t63)/(t36*2.0-udg4*2.0)+(t9*t28*t56*t66)/(t36*2.0-udg4*2.0));
		dstype t125 = nlg2*(t70+(t9*t28*t54*t63)/(t36*2.0-udg4*2.0)+(t9*t28*t56*t66)/(t36*2.0-udg4*2.0));
		dstype t129 = nlg2*(t70+(t9*t28*t60*t66)/(t36*2.0-udg4*2.0)+(t9*t28*t63*t65)/(t36*2.0-udg4*2.0));
		dstype t131 = -nlg1*(t67+(t9*t54*t63*t65*udg1)/(t36*2.0-udg4*2.0)+(t9*t56*t60*t66*udg1)/(t36*2.0-udg4*2.0));
		dstype t103 = t82*(-t21+t53+t8*t14*(t36-udg4));
		dstype t116 = t91*t93;
		dstype t117 = t34+t68+t85+t86;
		dstype t134 = t77+t110+t111;
		dstype t135 = t40+t89+t107+t109;
		dstype t140 = t129+t131;
		dstype t118 = nlg1*t117;
		dstype t119 = nlg2*t117;
		dstype t133 = t40+t76+t101+t103;
		dstype t138 = t88+t115+t116;
		dstype t136 = t119+t123;
		dstype t137 = t118+t125;
		f[0*ng+i] = param5/2.0+udg1/2.0-(t19*(nlg2*(t57+(t9*t60*t66*udg1)/(t36*2.0-udg4*2.0)+(t9*t63*t65*udg1)/(t36*2.0-udg4*2.0))+nlg1*(t58+(t9*t28*t63)/(t36*2.0-udg4*2.0)+(t9*t28*t66)/(t36*2.0-udg4*2.0))))/2.0+(t20*((t9*t10*t34)/(t36-udg4)+(t9*t63*udg1)/(t36*2.0-udg4*2.0)+(t9*t66*udg1)/(t36*2.0-udg4*2.0)))/2.0-(t18*(t114+nlg1*(t57+(t9*t60*t66*udg1)/(t36*2.0-udg4*2.0)+(t9*t63*t65*udg1)/(t36*2.0-udg4*2.0))))/2.0-(t17*((t9*t34*t72*udg1)/(t36*2.0-udg4*2.0)-(t9*t66*t91*udg1)/(t36*2.0-udg4*2.0)+(t9*t63*udg1*(-t21+t78+t8*t14*(t36-udg4)))/(t36*2.0-udg4*2.0)))/2.0;
		f[1*ng+i] = param6/2.0+udg2/2.0-(t17*(nlg2*t135+nlg1*t138))/2.0-(t18*(nlg2*t136-nlg1*t140))/2.0+(t19*(nlg1*t136+nlg2*t140))/2.0+(t20*(t114+nlg1*(t57+(t9*t54*t63*udg1)/(t36*2.0-udg4*2.0)+(t9*t56*t66*udg1)/(t36*2.0-udg4*2.0))))/2.0;
		f[2*ng+i] = param7/2.0+udg3/2.0+(t20*(nlg2*(t57+(t9*t54*t63*udg1)/(t36*2.0-udg4*2.0)+(t9*t56*t66*udg1)/(t36*2.0-udg4*2.0))+nlg1*(t58+(t9*t28*t63)/(t36*2.0-udg4*2.0)+(t9*t28*t66)/(t36*2.0-udg4*2.0))))/2.0-(t18*(nlg1*(nlg1*(t70+(t9*t28*t60*t66)/(t36*2.0-udg4*2.0)+(t9*t28*t63*t65)/(t36*2.0-udg4*2.0))+nlg2*(t67+(t9*t54*t63*t65*udg1)/(t36*2.0-udg4*2.0)+(t9*t56*t60*t66*udg1)/(t36*2.0-udg4*2.0)))-nlg2*t137))/2.0-(t19*(nlg2*(nlg1*(t70+(t9*t28*t60*t66)/(t36*2.0-udg4*2.0)+(t9*t28*t63*t65)/(t36*2.0-udg4*2.0))+nlg2*(t67+(t9*t54*t63*t65*udg1)/(t36*2.0-udg4*2.0)+(t9*t56*t60*t66*udg1)/(t36*2.0-udg4*2.0)))+nlg1*t137))/2.0+(t17*(nlg1*t135-nlg2*t138))/2.0;
		f[3*ng+i] = param8/2.0+udg4/2.0-(t17*(-t15*t29*t34+(t9*t34*t46*t72*udg1)/(t36*2.0-udg4*2.0)+(t9*t66*t91*udg1*(-t21+t53+t8*t14*(t36-udg4)))/(t36*2.0-udg4*2.0)+(t9*t63*t73*udg1*(-t21+t78+t8*t14*(t36-udg4)))/(t36*2.0-udg4*2.0)))/2.0+(t18*(nlg1*t134+nlg2*t133))/2.0-(t19*(nlg1*t133-nlg2*t134))/2.0+(t20*((t9*t10*t34*t46)/(t36-udg4)+(t9*t63*t73*udg1)/(t36*2.0-udg4*2.0)-(t9*t66*udg1*(-t21+t53+t8*t14*(t36-udg4)))/(t36*2.0-udg4*2.0)))/2.0;
	});
}

void KokkosUbou(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	if (ib == 1)
		KokkosUbou1(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 2)
		KokkosUbou2(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

