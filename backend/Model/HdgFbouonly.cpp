void HdgFbouonly1(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbouonly1", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype tau1 = tau[0];
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg4 = udg[3*ng+i];
		dstype udg9 = udg[8*ng+i];
		dstype udg10 = udg[9*ng+i];
		dstype udg12 = udg[11*ng+i];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		dstype uhg4 = uhg[3*ng+i];
		f[0*ng+i] = udg9+tau1*(udg1-uhg1);
		f[1*ng+i] = udg10+tau1*(udg2-uhg2);
		f[2*ng+i] = -tau1*uhg3;
		f[3*ng+i] = udg12+tau1*(udg4-uhg4);
	});
}

void HdgFbouonly2(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbouonly2", ng, KOKKOS_LAMBDA(const size_t i) {
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
		f[0*ng+i] = udg1/1.0E+2-uhg1/1.0E+2;
		f[1*ng+i] = (nlg2*(udg2-uhg2))/1.0E+2-(nlg1*(udg3-uhg3))/1.0E+2;
		f[2*ng+i] = nlg1*uhg2*(-1.0/1.0E+2)-(nlg2*uhg3)/1.0E+2;
		f[3*ng+i] = udg4/1.0E+2-uhg4/1.0E+2;
	});
}

void HdgFbouonly3(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbouonly3", ng, KOKKOS_LAMBDA(const size_t i) {
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
		f[0*ng+i] = udg1/1.0E+2-uhg1/1.0E+2;
		f[1*ng+i] = (nlg2*(udg2-uhg2))/1.0E+2-(nlg1*(udg3-uhg3))/1.0E+2;
		f[2*ng+i] = nlg1*uhg2*(-1.0/1.0E+2)-(nlg2*uhg3)/1.0E+2;
		f[3*ng+i] = udg4/1.0E+2-uhg4/1.0E+2;
	});
}

void HdgFbouonly4(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbouonly4", ng, KOKKOS_LAMBDA(const size_t i) {
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
		dstype t2 = nlg1*uhg2;
		dstype t3 = nlg1*uhg3;
		dstype t4 = nlg2*uhg2;
		dstype t5 = nlg2*uhg3;
		dstype t6 = nlg1*nlg1;
		dstype t7 = nlg2*nlg2;
		dstype t8 = udg2*udg2;
		dstype t9 = udg3*udg3;
		dstype t10 = uhg2*uhg2;
		dstype t11 = uhg3*uhg3;
		dstype t12 = 1.0/3.141592653589793;
		dstype t13 = -udg4;
		dstype t14 = 1.0/udg1;
		dstype t15 = uhg1*2.0E+1;
		dstype t17 = sqrt(3.141592653589793);
		dstype t18 = uhg1*1.0E+2;
		dstype t19 = sqrt(3.5E+1);
		dstype t30 = 3.141592653589793*1.269249157142259E+15;
		dstype t32 = uhg1*3.141592653589793*2.251799813685248E+16;
		dstype t16 = -t4;
		dstype t20 = t15-1.0;
		dstype t21 = t2+t5;
		dstype t22 = t6+t7;
		dstype t23 = t18-5.0;
		dstype t26 = (t8*t14)/2.0;
		dstype t27 = (t9*t14)/2.0;
		dstype t34 = t10*3.141592653589793*2.251799813685248E+16;
		dstype t35 = t11*3.141592653589793*2.251799813685248E+16;
		dstype t24 = atan(t23);
		dstype t25 = t3+t16;
		dstype t44 = t13+t26+t27+3.064634146341463;
		dstype t28 = t24*2.0;
		dstype t38 = t24*2.251799813685248E+15;
		dstype t41 = t24*uhg1*4.503599627370496E+16;
		dstype t29 = t28+3.141592653589793;
		dstype t39 = t38*uhg4;
		dstype t40 = -t38;
		dstype t45 = t12*t20*t29*1.125899906842624E+15;
		dstype t56 = t30+t32+t40+t41;
		dstype t46 = t45+2.395149063984883E+15;
		dstype t57 = t56*t56;
		dstype t47 = 1.0/t46;
		dstype t58 = 1.0/t57;
		dstype t48 = t47*t47;
		dstype t49 = t47*uhg4*6.305039478318694E+16;
		dstype t50 = t21*t47*4.503599627370496E+16;
		dstype t51 = t21*t47*4.503599627370496E+18;
		dstype t52 = nlg1*t25*t47*4.503599627370496E+16;
		dstype t53 = nlg2*t25*t47*4.503599627370496E+16;
		dstype t63 = -t58*(t34+t35+t39-uhg4*3.141592653589793*1.269249157142259E+15-uhg1*uhg4*3.141592653589793*2.251799813685248E+16-t24*uhg1*uhg4*4.503599627370496E+16);
		dstype t54 = tanh(t51);
		dstype t59 = t10*t48*4.056481920730334E+32;
		dstype t60 = t11*t48*4.056481920730334E+32;
		dstype t64 = sqrt(t63);
		dstype t65 = t17*t19*t64*2.68435456E+9;
		dstype t67 = t17*t19*t64*2.68435456E+7;
		dstype t69 = t17*t19*t21*t47*t64*1.208925819614629E+24;
		dstype t66 = -t65;
		dstype t68 = -t67;
		dstype t70 = t50+t67;
		dstype t72 = t51+t65;
		dstype t71 = t50+t68;
		dstype t73 = tanh(t72);
		dstype t74 = t51+t66;
		dstype t75 = tanh(t74);
		f[0*ng+i] = udg1-uhg1+(t12*t44*t57*(t54*-2.0+t73+t75))/(t34*2.522015791327478E+17+t35*2.522015791327478E+17+t39*2.522015791327478E+17-uhg4*3.141592653589793*3.201066417441869E+32-uhg1*uhg4*3.141592653589793*5.679074689022468E+33-t24*uhg1*uhg4*1.135814937804494E+34);
		f[1*ng+i] = udg2-uhg2-(t44*((t12*t57*t73*(t53-nlg1*t70))/(t34*1.261007895663739E+17+t35*1.261007895663739E+17+t39*1.261007895663739E+17-uhg4*3.141592653589793*1.600533208720934E+32-uhg1*uhg4*3.141592653589793*2.839537344511234E+33-t24*uhg1*uhg4*5.679074689022468E+33)+(t12*t57*t75*(t53-nlg1*t71))/(t34*1.261007895663739E+17+t35*1.261007895663739E+17+t39*1.261007895663739E+17-uhg4*3.141592653589793*1.600533208720934E+32-uhg1*uhg4*3.141592653589793*2.839537344511234E+33-t24*uhg1*uhg4*5.679074689022468E+33)+(t22*t54*t56*uhg2*(5.0/7.0))/(t34+t35+t39-uhg4*3.141592653589793*1.269249157142259E+15-uhg1*uhg4*3.141592653589793*2.251799813685248E+16-t24*uhg1*uhg4*4.503599627370496E+16)))/2.0;
		f[2*ng+i] = udg3-uhg3+(t44*((t12*t57*t73*(t52+nlg2*t70))/(t34*1.261007895663739E+17+t35*1.261007895663739E+17+t39*1.261007895663739E+17-uhg4*3.141592653589793*1.600533208720934E+32-uhg1*uhg4*3.141592653589793*2.839537344511234E+33-t24*uhg1*uhg4*5.679074689022468E+33)+(t12*t57*t75*(t52+nlg2*t71))/(t34*1.261007895663739E+17+t35*1.261007895663739E+17+t39*1.261007895663739E+17-uhg4*3.141592653589793*1.600533208720934E+32-uhg1*uhg4*3.141592653589793*2.839537344511234E+33-t24*uhg1*uhg4*5.679074689022468E+33)-(t22*t54*t56*uhg3*(5.0/7.0))/(t34+t35+t39-uhg4*3.141592653589793*1.269249157142259E+15-uhg1*uhg4*3.141592653589793*2.251799813685248E+16-t24*uhg1*uhg4*4.503599627370496E+16)))/2.0;
		f[3*ng+i] = udg4/2.0-uhg4-(t44*((t12*t57*t73*(t49-t59-t60+t69)*(-7.930164461608261E-18))/(t34+t35+t39-uhg4*3.141592653589793*1.269249157142259E+15-uhg1*uhg4*3.141592653589793*2.251799813685248E+16-t24*uhg1*uhg4*4.503599627370496E+16)+(t22*t54*3.141592653589793*(t10+t11)*1.608428438346606E+16)/(t34+t35+t39-uhg4*3.141592653589793*1.269249157142259E+15-uhg1*uhg4*3.141592653589793*2.251799813685248E+16-t24*uhg1*uhg4*4.503599627370496E+16)+(t12*t57*t75*(-t49+t59+t60+t69))/(t34*1.261007895663739E+17+t35*1.261007895663739E+17+t39*1.261007895663739E+17-uhg4*3.141592653589793*1.600533208720934E+32-uhg1*uhg4*3.141592653589793*2.839537344511234E+33-t24*uhg1*uhg4*5.679074689022468E+33)))/2.0+(t8*t14)/4.0+(t9*t14)/4.0+1.532317073170732;
	});
}

void HdgFbouonly(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	if (ib == 1)
		HdgFbouonly1(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 2)
		HdgFbouonly2(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 3)
		HdgFbouonly3(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 4)
		HdgFbouonly4(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

