void HdgFbouonly1(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbouonly1", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg4 = udg[3*ng+i];
		dstype udg9 = udg[8*ng+i];
		dstype udg10 = udg[9*ng+i];
		dstype udg11 = udg[10*ng+i];
		dstype udg12 = udg[11*ng+i];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		dstype uhg4 = uhg[3*ng+i];
		dstype t2 = 1.0/uhg1;
		dstype t3 = t2*t2;
		dstype t4 = t2*udg9*uhg2;
		f[0*ng+i] = udg1+udg9-uhg1;
		f[1*ng+i] = udg2-uhg2-t2*(t4-udg10);
		f[2*ng+i] = -uhg3;
		f[3*ng+i] = udg4-uhg4+t2*(udg12-t2*udg9*uhg4)*(2.0/5.0)-t3*uhg3*(udg11-t2*udg9*uhg3)*(2.0/5.0)+t3*uhg2*(t4-udg10)*(2.0/5.0);
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
		dstype t2 = -uhg2;
		dstype t3 = -uhg3;
		f[0*ng+i] = udg1-uhg1;
		f[1*ng+i] = nlg2*(t2+udg2)-nlg1*(t3+udg3);
		f[2*ng+i] = nlg1*t2+nlg2*t3;
		f[3*ng+i] = udg4-uhg4;
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
		dstype t2 = -uhg2;
		dstype t3 = -uhg3;
		f[0*ng+i] = udg1-uhg1;
		f[1*ng+i] = nlg2*(t2+udg2)-nlg1*(t3+udg3);
		f[2*ng+i] = nlg1*t2+nlg2*t3;
		f[3*ng+i] = udg4-uhg4;
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
		dstype t3 = nlg2*uhg3;
		dstype t4 = nlg1*nlg1;
		dstype t5 = nlg2*nlg2;
		dstype t6 = udg2*udg2;
		dstype t7 = udg3*udg3;
		dstype t8 = uhg1*uhg1;
		dstype t9 = uhg2*uhg2;
		dstype t10 = uhg3*uhg3;
		dstype t11 = uhg1*uhg4*2.0;
		dstype t12 = -udg4;
		dstype t13 = 1.0/udg1;
		dstype t14 = 1.0/uhg1;
		dstype t19 = uhg1*uhg4*1.4E+1;
		dstype t20 = sqrt(7.0);
		dstype t15 = 1.0/t8;
		dstype t16 = t9*7.0;
		dstype t17 = t10*7.0;
		dstype t18 = -t11;
		dstype t21 = nlg1*t2;
		dstype t22 = t4*uhg3;
		dstype t23 = t5*uhg2;
		dstype t24 = nlg2*t3;
		dstype t25 = t2*1.0E+2;
		dstype t26 = t3*1.0E+2;
		dstype t27 = -t19;
		dstype t28 = t2+t3;
		dstype t29 = t4+t5;
		dstype t30 = t14*uhg4*(7.0/5.0);
		dstype t32 = (t6*t13)/2.0;
		dstype t33 = (t7*t13)/2.0;
		dstype t31 = -t30;
		dstype t34 = (t9*t15)/5.0;
		dstype t35 = (t10*t15)/5.0;
		dstype t36 = t9+t10+t18;
		dstype t37 = t25+t26;
		dstype t39 = t16+t17+t27;
		dstype t52 = t12+t32+t33+3.064634146341463;
		dstype t38 = 1.0/t36;
		dstype t40 = t15*t36;
		dstype t41 = t14*t37;
		dstype t44 = t15*t39;
		dstype t42 = -t40;
		dstype t43 = tanh(t41);
		dstype t46 = -t44;
		dstype t45 = sqrt(t42);
		dstype t47 = sqrt(t46);
		dstype t48 = t47*uhg1*2.0E+1;
		dstype t50 = (nlg1*t47*uhg1)/5.0;
		dstype t51 = (nlg2*t47*uhg1)/5.0;
		dstype t53 = (t14*t20*t28*t45)/5.0;
		dstype t49 = -t48;
		dstype t54 = t37+t48;
		dstype t55 = t37+t49;
		dstype t56 = t14*t54;
		dstype t57 = tanh(t56);
		dstype t58 = t14*t55;
		dstype t59 = tanh(t58);
		f[0*ng+i] = udg1-uhg1+t8*t38*t52*(t43*-2.0+t57+t59)*(5.0/1.4E+1);
		f[1*ng+i] = udg2-uhg2+(t52*(t38*t57*uhg1*(t21+t23+t50)*(5.0/7.0)+t38*t59*uhg1*(t21+t23-t50)*(5.0/7.0)-t29*t38*t43*uhg1*uhg2*(1.0E+1/7.0)))/2.0;
		f[2*ng+i] = udg3-uhg3+(t52*(t38*t57*uhg1*(t22+t24+t51)*(5.0/7.0)+t38*t59*uhg1*(t22+t24-t51)*(5.0/7.0)-t29*t38*t43*uhg1*uhg3*(1.0E+1/7.0)))/2.0;
		f[3*ng+i] = udg4/2.0-uhg4+(t6*t13)/4.0+(t7*t13)/4.0-(t52*(t29*t38*t43*(t9+t10)*(5.0/7.0)-t8*t38*t57*(t30-t34-t35+t53)*(5.0/7.0)+t8*t38*t59*(t31+t34+t35+t53)*(5.0/7.0)))/2.0+1.532317073170732;
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

