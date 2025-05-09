void KokkosFlux1(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Flux1", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype param2 = param[1];
		dstype param3 = param[2];
		dstype param4 = param[3];
		dstype param10 = param[9];
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
		dstype odg1 = odg[0*ng+i];
		dstype t2 = 1.0/3.141592653589793;
		dstype t3 = udg1*1.0E3;
		dstype t4 = t3-1.0E1;
		dstype t5 = atan(t4);
		dstype t6 = t2*t5;
		dstype t7 = t6+1.0/2.0;
		dstype t8 = udg1-1.0/1.0E2;
		dstype t9 = t7*t8*1.0E3;
		dstype t10 = t9+3.183097800805168E-1;
		dstype t12 = t7*t8;
		dstype t11 = t12+3.183097800805168E-4;
		dstype t13 = t12+1.031830978008052E-2;
		dstype t14 = 1.0/(t13*t13);
		dstype t15 = udg2*udg2;
		dstype t16 = t14*t15*(1.0/2.0);
		dstype t17 = udg3*udg3;
		dstype t18 = t14*t17*(1.0/2.0);
		dstype t19 = t16+t18;
		dstype t31 = t13*t19;
		dstype t20 = -t31+udg4;
		dstype t21 = param1-1.0;
		dstype t22 = 1.0/t13;
		dstype t23 = atan(t10);
		dstype t24 = t2*t23;
		dstype t25 = t11*t11;
		dstype t26 = t25*1.0E6;
		dstype t27 = t26+1.0;
		dstype t28 = 1.0/t27;
		dstype t29 = t2*t10*t28;
		dstype t30 = t24+t29+1.0/2.0;
		dstype t32 = t20*t21*1.0E3;
		dstype t33 = t32-1.0;
		dstype t34 = atan(t33);
		dstype t35 = t2*t34;
		dstype t36 = t35+1.0/2.0;
		dstype t37 = t20*t21;
		dstype t38 = t37-1.0/1.0E3;
		dstype t39 = t36*t38;
		dstype t40 = param4*param4;
		dstype t41 = t39+1.318309780080517E-3;
		dstype t42 = 1.0/param2;
		dstype t43 = param1*param10*t22*t40*t41;
		dstype t44 = t43+5.52E2/5.0;
		dstype t45 = 1.0/t44;
		dstype t46 = param10+5.52E2/5.0;
		dstype t47 = param1*t22*t40*t41;
		dstype t48 = pow(t47,3.0/2.0);
		dstype t62 = t22*t30*udg3*udg5;
		dstype t49 = -t62+udg7;
		dstype t50 = t22*t49;
		dstype t64 = t22*t30*udg2*udg9;
		dstype t51 = -t64+udg10;
		dstype t52 = t22*t51;
		dstype t53 = t50+t52;
		dstype t61 = t22*t30*udg2*udg5;
		dstype t54 = -t61+udg6;
		dstype t55 = t22*t54*2.0;
		dstype t66 = t22*t30*udg3*udg9;
		dstype t56 = -t66+udg11;
		dstype t57 = t55-t22*t56;
		dstype t58 = t39+3.183097800805168E-4;
		dstype t59 = t36*t38*1.0E3;
		dstype t60 = t59+3.183097800805168E-1;
		dstype t63 = t22*udg2*udg3;
		dstype t65 = t42*t45*t46*t48*t53;
		dstype t67 = t22*udg4;
		dstype t68 = t22*t41;
		dstype t69 = t67+t68;
		dstype t70 = t22*t54;
		dstype t71 = t70-t22*t56*2.0;
		dstype t72 = 1.0/param3;
		dstype t73 = atan(t60);
		dstype t74 = t2*t73;
		dstype t75 = t58*t58;
		dstype t76 = t75*1.0E6;
		dstype t77 = t76+1.0;
		dstype t78 = 1.0/t77;
		dstype t79 = t2*t60*t78;
		dstype t80 = t74+t79+1.0/2.0;
		dstype t81 = 1.0/t21;
		f[0*ng+i] = udg2+odg1*t30*udg5;
		f[1*ng+i] = t39+odg1*udg6+t15*t22+t42*t45*t46*t48*t57*(2.0/3.0)+1.318309780080517E-3;
		f[2*ng+i] = t63+t65+odg1*udg7;
		f[3*ng+i] = odg1*udg8+t69*udg2+t22*t42*t45*t46*t48*t53*udg3+t22*t42*t45*t46*t48*t57*udg2*(2.0/3.0)-param1*t14*t42*t45*t46*t48*t72*t81*(t30*t41*udg5+t13*t21*t80*(-udg8+t13*(t14*t49*udg3+t14*t54*udg2)+t19*t30*udg5));
		f[4*ng+i] = udg3+odg1*t30*udg9;
		f[5*ng+i] = t63+t65+odg1*udg10;
		f[6*ng+i] = t39+odg1*udg11+t17*t22-t42*t45*t46*t48*t71*(2.0/3.0)+1.318309780080517E-3;
		f[7*ng+i] = odg1*udg12+t69*udg3+t22*t42*t45*t46*t48*t53*udg2-t22*t42*t45*t46*t48*t71*udg3*(2.0/3.0)-param1*t14*t42*t45*t46*t48*t72*t81*(t30*t41*udg9+t13*t21*t80*(-udg12+t13*(t14*t51*udg2+t14*t56*udg3)+t19*t30*udg9));
	});
}

