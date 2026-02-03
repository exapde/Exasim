void HdgFbou1(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbou1", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param5 = param[4];
		dstype param6 = param[5];
		dstype param7 = param[6];
		dstype param8 = param[7];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		dstype uhg4 = uhg[3*ng+i];
		{
		f[0*ng+i] = param5-uhg1;
		f[1*ng+i] = param6-uhg2;
		f[2*ng+i] = param7-uhg3;
		f[3*ng+i] = param8-uhg4;
		}
		{
		f_udg[0*ng+i] = 0.0;
		f_udg[1*ng+i] = 0.0;
		f_udg[2*ng+i] = 0.0;
		f_udg[3*ng+i] = 0.0;
		f_udg[4*ng+i] = 0.0;
		f_udg[5*ng+i] = 0.0;
		f_udg[6*ng+i] = 0.0;
		f_udg[7*ng+i] = 0.0;
		f_udg[8*ng+i] = 0.0;
		f_udg[9*ng+i] = 0.0;
		f_udg[10*ng+i] = 0.0;
		f_udg[11*ng+i] = 0.0;
		f_udg[12*ng+i] = 0.0;
		f_udg[13*ng+i] = 0.0;
		f_udg[14*ng+i] = 0.0;
		f_udg[15*ng+i] = 0.0;
		f_udg[16*ng+i] = 0.0;
		f_udg[17*ng+i] = 0.0;
		f_udg[18*ng+i] = 0.0;
		f_udg[19*ng+i] = 0.0;
		f_udg[20*ng+i] = 0.0;
		f_udg[21*ng+i] = 0.0;
		f_udg[22*ng+i] = 0.0;
		f_udg[23*ng+i] = 0.0;
		f_udg[24*ng+i] = 0.0;
		f_udg[25*ng+i] = 0.0;
		f_udg[26*ng+i] = 0.0;
		f_udg[27*ng+i] = 0.0;
		f_udg[28*ng+i] = 0.0;
		f_udg[29*ng+i] = 0.0;
		f_udg[30*ng+i] = 0.0;
		f_udg[31*ng+i] = 0.0;
		f_udg[32*ng+i] = 0.0;
		f_udg[33*ng+i] = 0.0;
		f_udg[34*ng+i] = 0.0;
		f_udg[35*ng+i] = 0.0;
		f_udg[36*ng+i] = 0.0;
		f_udg[37*ng+i] = 0.0;
		f_udg[38*ng+i] = 0.0;
		f_udg[39*ng+i] = 0.0;
		f_udg[40*ng+i] = 0.0;
		f_udg[41*ng+i] = 0.0;
		f_udg[42*ng+i] = 0.0;
		f_udg[43*ng+i] = 0.0;
		f_udg[44*ng+i] = 0.0;
		f_udg[45*ng+i] = 0.0;
		f_udg[46*ng+i] = 0.0;
		f_udg[47*ng+i] = 0.0;
		}
		{
		f_uhg[0*ng+i] = -1.0;
		f_uhg[1*ng+i] = 0.0;
		f_uhg[2*ng+i] = 0.0;
		f_uhg[3*ng+i] = 0.0;
		f_uhg[4*ng+i] = 0.0;
		f_uhg[5*ng+i] = -1.0;
		f_uhg[6*ng+i] = 0.0;
		f_uhg[7*ng+i] = 0.0;
		f_uhg[8*ng+i] = 0.0;
		f_uhg[9*ng+i] = 0.0;
		f_uhg[10*ng+i] = -1.0;
		f_uhg[11*ng+i] = 0.0;
		f_uhg[12*ng+i] = 0.0;
		f_uhg[13*ng+i] = 0.0;
		f_uhg[14*ng+i] = 0.0;
		f_uhg[15*ng+i] = -1.0;
		}
	});
}

void HdgFbou2(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbou2", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg3 = udg[2*ng+i];
		dstype udg4 = udg[3*ng+i];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		dstype uhg4 = uhg[3*ng+i];
		{
		f[0*ng+i] = udg1-uhg1;
		f[1*ng+i] = udg2-uhg2;
		f[2*ng+i] = udg3-uhg3;
		f[3*ng+i] = udg4-uhg4;
		}
		{
		f_udg[0*ng+i] = 1.0;
		f_udg[1*ng+i] = 0.0;
		f_udg[2*ng+i] = 0.0;
		f_udg[3*ng+i] = 0.0;
		f_udg[4*ng+i] = 0.0;
		f_udg[5*ng+i] = 1.0;
		f_udg[6*ng+i] = 0.0;
		f_udg[7*ng+i] = 0.0;
		f_udg[8*ng+i] = 0.0;
		f_udg[9*ng+i] = 0.0;
		f_udg[10*ng+i] = 1.0;
		f_udg[11*ng+i] = 0.0;
		f_udg[12*ng+i] = 0.0;
		f_udg[13*ng+i] = 0.0;
		f_udg[14*ng+i] = 0.0;
		f_udg[15*ng+i] = 1.0;
		f_udg[16*ng+i] = 0.0;
		f_udg[17*ng+i] = 0.0;
		f_udg[18*ng+i] = 0.0;
		f_udg[19*ng+i] = 0.0;
		f_udg[20*ng+i] = 0.0;
		f_udg[21*ng+i] = 0.0;
		f_udg[22*ng+i] = 0.0;
		f_udg[23*ng+i] = 0.0;
		f_udg[24*ng+i] = 0.0;
		f_udg[25*ng+i] = 0.0;
		f_udg[26*ng+i] = 0.0;
		f_udg[27*ng+i] = 0.0;
		f_udg[28*ng+i] = 0.0;
		f_udg[29*ng+i] = 0.0;
		f_udg[30*ng+i] = 0.0;
		f_udg[31*ng+i] = 0.0;
		f_udg[32*ng+i] = 0.0;
		f_udg[33*ng+i] = 0.0;
		f_udg[34*ng+i] = 0.0;
		f_udg[35*ng+i] = 0.0;
		f_udg[36*ng+i] = 0.0;
		f_udg[37*ng+i] = 0.0;
		f_udg[38*ng+i] = 0.0;
		f_udg[39*ng+i] = 0.0;
		f_udg[40*ng+i] = 0.0;
		f_udg[41*ng+i] = 0.0;
		f_udg[42*ng+i] = 0.0;
		f_udg[43*ng+i] = 0.0;
		f_udg[44*ng+i] = 0.0;
		f_udg[45*ng+i] = 0.0;
		f_udg[46*ng+i] = 0.0;
		f_udg[47*ng+i] = 0.0;
		}
		{
		f_uhg[0*ng+i] = -1.0;
		f_uhg[1*ng+i] = 0.0;
		f_uhg[2*ng+i] = 0.0;
		f_uhg[3*ng+i] = 0.0;
		f_uhg[4*ng+i] = 0.0;
		f_uhg[5*ng+i] = -1.0;
		f_uhg[6*ng+i] = 0.0;
		f_uhg[7*ng+i] = 0.0;
		f_uhg[8*ng+i] = 0.0;
		f_uhg[9*ng+i] = 0.0;
		f_uhg[10*ng+i] = -1.0;
		f_uhg[11*ng+i] = 0.0;
		f_uhg[12*ng+i] = 0.0;
		f_uhg[13*ng+i] = 0.0;
		f_uhg[14*ng+i] = 0.0;
		f_uhg[15*ng+i] = -1.0;
		}
	});
}

void HdgFbou3(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbou3", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype param9 = param[8];
		dstype param10 = param[9];
		dstype param11 = param[10];
		dstype udg1 = udg[0*ng+i];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		dstype uhg4 = uhg[3*ng+i];
		{
		f[0*ng+i] = udg1-uhg1;
		f[1*ng+i] = -uhg2;
		f[2*ng+i] = -uhg3;
		f[3*ng+i] = -uhg4+(param9*param11*uhg1)/param10;
		}
		{
		f_udg[0*ng+i] = 1.0;
		f_udg[1*ng+i] = 0.0;
		f_udg[2*ng+i] = 0.0;
		f_udg[3*ng+i] = 0.0;
		f_udg[4*ng+i] = 0.0;
		f_udg[5*ng+i] = 0.0;
		f_udg[6*ng+i] = 0.0;
		f_udg[7*ng+i] = 0.0;
		f_udg[8*ng+i] = 0.0;
		f_udg[9*ng+i] = 0.0;
		f_udg[10*ng+i] = 0.0;
		f_udg[11*ng+i] = 0.0;
		f_udg[12*ng+i] = 0.0;
		f_udg[13*ng+i] = 0.0;
		f_udg[14*ng+i] = 0.0;
		f_udg[15*ng+i] = 0.0;
		f_udg[16*ng+i] = 0.0;
		f_udg[17*ng+i] = 0.0;
		f_udg[18*ng+i] = 0.0;
		f_udg[19*ng+i] = 0.0;
		f_udg[20*ng+i] = 0.0;
		f_udg[21*ng+i] = 0.0;
		f_udg[22*ng+i] = 0.0;
		f_udg[23*ng+i] = 0.0;
		f_udg[24*ng+i] = 0.0;
		f_udg[25*ng+i] = 0.0;
		f_udg[26*ng+i] = 0.0;
		f_udg[27*ng+i] = 0.0;
		f_udg[28*ng+i] = 0.0;
		f_udg[29*ng+i] = 0.0;
		f_udg[30*ng+i] = 0.0;
		f_udg[31*ng+i] = 0.0;
		f_udg[32*ng+i] = 0.0;
		f_udg[33*ng+i] = 0.0;
		f_udg[34*ng+i] = 0.0;
		f_udg[35*ng+i] = 0.0;
		f_udg[36*ng+i] = 0.0;
		f_udg[37*ng+i] = 0.0;
		f_udg[38*ng+i] = 0.0;
		f_udg[39*ng+i] = 0.0;
		f_udg[40*ng+i] = 0.0;
		f_udg[41*ng+i] = 0.0;
		f_udg[42*ng+i] = 0.0;
		f_udg[43*ng+i] = 0.0;
		f_udg[44*ng+i] = 0.0;
		f_udg[45*ng+i] = 0.0;
		f_udg[46*ng+i] = 0.0;
		f_udg[47*ng+i] = 0.0;
		}
		{
		f_uhg[0*ng+i] = -1.0;
		f_uhg[1*ng+i] = 0.0;
		f_uhg[2*ng+i] = 0.0;
		f_uhg[3*ng+i] = (param9*param11)/param10;
		f_uhg[4*ng+i] = 0.0;
		f_uhg[5*ng+i] = -1.0;
		f_uhg[6*ng+i] = 0.0;
		f_uhg[7*ng+i] = 0.0;
		f_uhg[8*ng+i] = 0.0;
		f_uhg[9*ng+i] = 0.0;
		f_uhg[10*ng+i] = -1.0;
		f_uhg[11*ng+i] = 0.0;
		f_uhg[12*ng+i] = 0.0;
		f_uhg[13*ng+i] = 0.0;
		f_uhg[14*ng+i] = 0.0;
		f_uhg[15*ng+i] = -1.0;
		}
	});
}

void HdgFbou4(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbou4", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype tau1 = tau[0];
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
		{
		dstype t2 = uhg2*uhg2;
		dstype t3 = uhg3*uhg3;
		dstype t4 = param1-1.0;
		dstype t5 = 1.0/uhg1;
		dstype t6 = t5*t5;
		dstype t7 = t5*udg9*uhg2;
		dstype t8 = t5*udg9*uhg3;
		dstype t11 = (t2*t6)/2.0;
		dstype t12 = (t3*t6)/2.0;
		dstype t15 = t11+t12;
		f[0*ng+i] = udg9+tau1*(udg1-uhg1);
		f[1*ng+i] = -t5*(t7-udg10)+tau1*(udg2-uhg2);
		f[2*ng+i] = -tau1*uhg3-t5*(t8-udg11);
		f[3*ng+i] = tau1*(udg4-uhg4)-(t6*(t4*udg9*(uhg4-t15*uhg1)-t4*uhg1*(udg12+uhg1*(t6*uhg2*(t7-udg10)+t6*uhg3*(t8-udg11))-t15*udg9)))/t4;
		}
		{
		dstype t2 = uhg2*uhg2;
		dstype t3 = uhg3*uhg3;
		dstype t4 = param1-1.0;
		dstype t5 = 1.0/uhg1;
		dstype t6 = t5*t5;
		dstype t7 = t5*t5*t5;
		dstype t8 = t6*uhg2;
		dstype t9 = t6*uhg3;
		dstype t12 = (t2*t6)/2.0;
		dstype t13 = (t3*t6)/2.0;
		dstype t10 = -t8;
		dstype t11 = -t9;
		f_udg[0*ng+i] = tau1;
		f_udg[1*ng+i] = 0.0;
		f_udg[2*ng+i] = 0.0;
		f_udg[3*ng+i] = 0.0;
		f_udg[4*ng+i] = 0.0;
		f_udg[5*ng+i] = tau1;
		f_udg[6*ng+i] = 0.0;
		f_udg[7*ng+i] = 0.0;
		f_udg[8*ng+i] = 0.0;
		f_udg[9*ng+i] = 0.0;
		f_udg[10*ng+i] = 0.0;
		f_udg[11*ng+i] = 0.0;
		f_udg[12*ng+i] = 0.0;
		f_udg[13*ng+i] = 0.0;
		f_udg[14*ng+i] = 0.0;
		f_udg[15*ng+i] = tau1;
		f_udg[16*ng+i] = 0.0;
		f_udg[17*ng+i] = 0.0;
		f_udg[18*ng+i] = 0.0;
		f_udg[19*ng+i] = 0.0;
		f_udg[20*ng+i] = 0.0;
		f_udg[21*ng+i] = 0.0;
		f_udg[22*ng+i] = 0.0;
		f_udg[23*ng+i] = 0.0;
		f_udg[24*ng+i] = 0.0;
		f_udg[25*ng+i] = 0.0;
		f_udg[26*ng+i] = 0.0;
		f_udg[27*ng+i] = 0.0;
		f_udg[28*ng+i] = 0.0;
		f_udg[29*ng+i] = 0.0;
		f_udg[30*ng+i] = 0.0;
		f_udg[31*ng+i] = 0.0;
		f_udg[32*ng+i] = 1.0;
		f_udg[33*ng+i] = t10;
		f_udg[34*ng+i] = t11;
		f_udg[35*ng+i] = -(t6*(t4*(uhg4-uhg1*(t12+t13))+t4*uhg1*(t12+t13-uhg1*(t2*t7+t3*t7))))/t4;
		f_udg[36*ng+i] = 0.0;
		f_udg[37*ng+i] = t5;
		f_udg[38*ng+i] = 0.0;
		f_udg[39*ng+i] = t10;
		f_udg[40*ng+i] = 0.0;
		f_udg[41*ng+i] = 0.0;
		f_udg[42*ng+i] = t5;
		f_udg[43*ng+i] = t11;
		f_udg[44*ng+i] = 0.0;
		f_udg[45*ng+i] = 0.0;
		f_udg[46*ng+i] = 0.0;
		f_udg[47*ng+i] = t5;
		}
		{
		dstype t2 = uhg2*uhg2;
		dstype t3 = uhg3*uhg3;
		dstype t4 = param1-1.0;
		dstype t5 = -tau1;
		dstype t6 = -udg12;
		dstype t7 = 1.0/uhg1;
		dstype t8 = t7*t7;
		dstype t9 = t7*t7*t7;
		dstype t12 = t7*udg9*uhg2;
		dstype t13 = t7*udg9*uhg3;
		dstype t16 = 1.0/t4;
		dstype t10 = t8*t8;
		dstype t11 = t8*udg9;
		dstype t14 = t9*udg9*uhg2;
		dstype t15 = t9*udg9*uhg3;
		dstype t17 = t2*t9;
		dstype t18 = t3*t9;
		dstype t22 = (t2*t8)/2.0;
		dstype t23 = (t3*t8)/2.0;
		dstype t35 = -uhg1*(t8*uhg2*(t12-udg10)+t8*uhg3*(t13-udg11));
		dstype t19 = -t11;
		dstype t29 = t17+t18;
		dstype t32 = t22+t23;
		dstype t26 = t5+t19;
		dstype t33 = t32*udg9;
		dstype t36 = t6+t33+t35;
		f_uhg[0*ng+i] = t5;
		f_uhg[1*ng+i] = t14+t8*(t12-udg10);
		f_uhg[2*ng+i] = t15+t8*(t13-udg11);
		f_uhg[3*ng+i] = t8*t16*(-t4*t36+t4*uhg1*(t29*udg9-uhg1*(t2*t10*udg9+t3*t10*udg9+t9*uhg2*(t12-udg10)*2.0+t9*uhg3*(t13-udg11)*2.0)+t8*uhg2*(t12-udg10)+t8*uhg3*(t13-udg11))+t4*udg9*(t32-t29*uhg1))+t9*t16*(t4*udg9*(uhg4-t32*uhg1)+t4*t36*uhg1)*2.0;
		f_uhg[4*ng+i] = 0.0;
		f_uhg[5*ng+i] = t26;
		f_uhg[6*ng+i] = 0.0;
		f_uhg[7*ng+i] = t8*t16*(t4*t12+t4*uhg1*(uhg1*(t14+t8*(t12-udg10))+t19*uhg2));
		f_uhg[8*ng+i] = 0.0;
		f_uhg[9*ng+i] = 0.0;
		f_uhg[10*ng+i] = t26;
		f_uhg[11*ng+i] = t8*t16*(t4*t13+t4*uhg1*(uhg1*(t15+t8*(t13-udg11))+t19*uhg3));
		f_uhg[12*ng+i] = 0.0;
		f_uhg[13*ng+i] = 0.0;
		f_uhg[14*ng+i] = 0.0;
		f_uhg[15*ng+i] = t26;
		}
	});
}

void HdgFbou5(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbou5", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype tau1 = tau[0];
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
		{
		dstype t2 = nlg1*udg2;
		dstype t3 = nlg2*udg3;
		dstype t4 = t2+t3;
		f[0*ng+i] = tau1*(udg1-uhg1);
		f[1*ng+i] = -tau1*(-udg2+uhg2+nlg1*t4);
		f[2*ng+i] = -tau1*(-udg3+uhg3+nlg2*t4);
		f[3*ng+i] = tau1*(udg4-uhg4);
		}
		{
		dstype t2 = nlg1*nlg2*tau1;
		dstype t3 = -t2;
		f_udg[0*ng+i] = tau1;
		f_udg[1*ng+i] = 0.0;
		f_udg[2*ng+i] = 0.0;
		f_udg[3*ng+i] = 0.0;
		f_udg[4*ng+i] = 0.0;
		f_udg[5*ng+i] = -tau1*(nlg1*nlg1-1.0);
		f_udg[6*ng+i] = t3;
		f_udg[7*ng+i] = 0.0;
		f_udg[8*ng+i] = 0.0;
		f_udg[9*ng+i] = t3;
		f_udg[10*ng+i] = -tau1*(nlg2*nlg2-1.0);
		f_udg[11*ng+i] = 0.0;
		f_udg[12*ng+i] = 0.0;
		f_udg[13*ng+i] = 0.0;
		f_udg[14*ng+i] = 0.0;
		f_udg[15*ng+i] = tau1;
		f_udg[16*ng+i] = 0.0;
		f_udg[17*ng+i] = 0.0;
		f_udg[18*ng+i] = 0.0;
		f_udg[19*ng+i] = 0.0;
		f_udg[20*ng+i] = 0.0;
		f_udg[21*ng+i] = 0.0;
		f_udg[22*ng+i] = 0.0;
		f_udg[23*ng+i] = 0.0;
		f_udg[24*ng+i] = 0.0;
		f_udg[25*ng+i] = 0.0;
		f_udg[26*ng+i] = 0.0;
		f_udg[27*ng+i] = 0.0;
		f_udg[28*ng+i] = 0.0;
		f_udg[29*ng+i] = 0.0;
		f_udg[30*ng+i] = 0.0;
		f_udg[31*ng+i] = 0.0;
		f_udg[32*ng+i] = 0.0;
		f_udg[33*ng+i] = 0.0;
		f_udg[34*ng+i] = 0.0;
		f_udg[35*ng+i] = 0.0;
		f_udg[36*ng+i] = 0.0;
		f_udg[37*ng+i] = 0.0;
		f_udg[38*ng+i] = 0.0;
		f_udg[39*ng+i] = 0.0;
		f_udg[40*ng+i] = 0.0;
		f_udg[41*ng+i] = 0.0;
		f_udg[42*ng+i] = 0.0;
		f_udg[43*ng+i] = 0.0;
		f_udg[44*ng+i] = 0.0;
		f_udg[45*ng+i] = 0.0;
		f_udg[46*ng+i] = 0.0;
		f_udg[47*ng+i] = 0.0;
		}
		{
		dstype t2 = -tau1;
		f_uhg[0*ng+i] = t2;
		f_uhg[1*ng+i] = 0.0;
		f_uhg[2*ng+i] = 0.0;
		f_uhg[3*ng+i] = 0.0;
		f_uhg[4*ng+i] = 0.0;
		f_uhg[5*ng+i] = t2;
		f_uhg[6*ng+i] = 0.0;
		f_uhg[7*ng+i] = 0.0;
		f_uhg[8*ng+i] = 0.0;
		f_uhg[9*ng+i] = 0.0;
		f_uhg[10*ng+i] = t2;
		f_uhg[11*ng+i] = 0.0;
		f_uhg[12*ng+i] = 0.0;
		f_uhg[13*ng+i] = 0.0;
		f_uhg[14*ng+i] = 0.0;
		f_uhg[15*ng+i] = t2;
		}
	});
}

void HdgFbou6(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbou6", ng, KOKKOS_LAMBDA(const size_t i) {
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
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		dstype uhg4 = uhg[3*ng+i];
		dstype nlg1 = nlg[0*ng+i];
		dstype nlg2 = nlg[1*ng+i];
		{
		f[0*ng+i] = nlg1*udg5+nlg2*udg9+tau1*(udg1-uhg1);
		f[1*ng+i] = nlg1*udg6+nlg2*udg10+tau1*(udg2-uhg2);
		f[2*ng+i] = nlg1*udg7+nlg2*udg11+tau1*(udg3-uhg3);
		f[3*ng+i] = nlg1*udg8+nlg2*udg12+tau1*(udg4-uhg4);
		}
		{
		f_udg[0*ng+i] = tau1;
		f_udg[1*ng+i] = 0.0;
		f_udg[2*ng+i] = 0.0;
		f_udg[3*ng+i] = 0.0;
		f_udg[4*ng+i] = 0.0;
		f_udg[5*ng+i] = tau1;
		f_udg[6*ng+i] = 0.0;
		f_udg[7*ng+i] = 0.0;
		f_udg[8*ng+i] = 0.0;
		f_udg[9*ng+i] = 0.0;
		f_udg[10*ng+i] = tau1;
		f_udg[11*ng+i] = 0.0;
		f_udg[12*ng+i] = 0.0;
		f_udg[13*ng+i] = 0.0;
		f_udg[14*ng+i] = 0.0;
		f_udg[15*ng+i] = tau1;
		f_udg[16*ng+i] = nlg1;
		f_udg[17*ng+i] = 0.0;
		f_udg[18*ng+i] = 0.0;
		f_udg[19*ng+i] = 0.0;
		f_udg[20*ng+i] = 0.0;
		f_udg[21*ng+i] = nlg1;
		f_udg[22*ng+i] = 0.0;
		f_udg[23*ng+i] = 0.0;
		f_udg[24*ng+i] = 0.0;
		f_udg[25*ng+i] = 0.0;
		f_udg[26*ng+i] = nlg1;
		f_udg[27*ng+i] = 0.0;
		f_udg[28*ng+i] = 0.0;
		f_udg[29*ng+i] = 0.0;
		f_udg[30*ng+i] = 0.0;
		f_udg[31*ng+i] = nlg1;
		f_udg[32*ng+i] = nlg2;
		f_udg[33*ng+i] = 0.0;
		f_udg[34*ng+i] = 0.0;
		f_udg[35*ng+i] = 0.0;
		f_udg[36*ng+i] = 0.0;
		f_udg[37*ng+i] = nlg2;
		f_udg[38*ng+i] = 0.0;
		f_udg[39*ng+i] = 0.0;
		f_udg[40*ng+i] = 0.0;
		f_udg[41*ng+i] = 0.0;
		f_udg[42*ng+i] = nlg2;
		f_udg[43*ng+i] = 0.0;
		f_udg[44*ng+i] = 0.0;
		f_udg[45*ng+i] = 0.0;
		f_udg[46*ng+i] = 0.0;
		f_udg[47*ng+i] = nlg2;
		}
		{
		dstype t2 = -tau1;
		f_uhg[0*ng+i] = t2;
		f_uhg[1*ng+i] = 0.0;
		f_uhg[2*ng+i] = 0.0;
		f_uhg[3*ng+i] = 0.0;
		f_uhg[4*ng+i] = 0.0;
		f_uhg[5*ng+i] = t2;
		f_uhg[6*ng+i] = 0.0;
		f_uhg[7*ng+i] = 0.0;
		f_uhg[8*ng+i] = 0.0;
		f_uhg[9*ng+i] = 0.0;
		f_uhg[10*ng+i] = t2;
		f_uhg[11*ng+i] = 0.0;
		f_uhg[12*ng+i] = 0.0;
		f_uhg[13*ng+i] = 0.0;
		f_uhg[14*ng+i] = 0.0;
		f_uhg[15*ng+i] = t2;
		}
	});
}

void HdgFbou7(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbou7", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype udg1 = udg[0*ng+i];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		{
		f[0*ng+i] = udg1-uhg1;
		f[1*ng+i] = -uhg2;
		f[2*ng+i] = -uhg3;
		f[3*ng+i] = 0.0;
		}
		{
		f_udg[0*ng+i] = 1.0;
		f_udg[1*ng+i] = 0.0;
		f_udg[2*ng+i] = 0.0;
		f_udg[3*ng+i] = 0.0;
		f_udg[4*ng+i] = 0.0;
		f_udg[5*ng+i] = 0.0;
		f_udg[6*ng+i] = 0.0;
		f_udg[7*ng+i] = 0.0;
		f_udg[8*ng+i] = 0.0;
		f_udg[9*ng+i] = 0.0;
		f_udg[10*ng+i] = 0.0;
		f_udg[11*ng+i] = 0.0;
		f_udg[12*ng+i] = 0.0;
		f_udg[13*ng+i] = 0.0;
		f_udg[14*ng+i] = 0.0;
		f_udg[15*ng+i] = 0.0;
		f_udg[16*ng+i] = 0.0;
		f_udg[17*ng+i] = 0.0;
		f_udg[18*ng+i] = 0.0;
		f_udg[19*ng+i] = 0.0;
		f_udg[20*ng+i] = 0.0;
		f_udg[21*ng+i] = 0.0;
		f_udg[22*ng+i] = 0.0;
		f_udg[23*ng+i] = 0.0;
		f_udg[24*ng+i] = 0.0;
		f_udg[25*ng+i] = 0.0;
		f_udg[26*ng+i] = 0.0;
		f_udg[27*ng+i] = 0.0;
		f_udg[28*ng+i] = 0.0;
		f_udg[29*ng+i] = 0.0;
		f_udg[30*ng+i] = 0.0;
		f_udg[31*ng+i] = 0.0;
		f_udg[32*ng+i] = 0.0;
		f_udg[33*ng+i] = 0.0;
		f_udg[34*ng+i] = 0.0;
		f_udg[35*ng+i] = 0.0;
		f_udg[36*ng+i] = 0.0;
		f_udg[37*ng+i] = 0.0;
		f_udg[38*ng+i] = 0.0;
		f_udg[39*ng+i] = 0.0;
		f_udg[40*ng+i] = 0.0;
		f_udg[41*ng+i] = 0.0;
		f_udg[42*ng+i] = 0.0;
		f_udg[43*ng+i] = 0.0;
		f_udg[44*ng+i] = 0.0;
		f_udg[45*ng+i] = 0.0;
		f_udg[46*ng+i] = 0.0;
		f_udg[47*ng+i] = 0.0;
		}
		{
		f_uhg[0*ng+i] = -1.0;
		f_uhg[1*ng+i] = 0.0;
		f_uhg[2*ng+i] = 0.0;
		f_uhg[3*ng+i] = 0.0;
		f_uhg[4*ng+i] = 0.0;
		f_uhg[5*ng+i] = -1.0;
		f_uhg[6*ng+i] = 0.0;
		f_uhg[7*ng+i] = 0.0;
		f_uhg[8*ng+i] = 0.0;
		f_uhg[9*ng+i] = 0.0;
		f_uhg[10*ng+i] = -1.0;
		f_uhg[11*ng+i] = 0.0;
		f_uhg[12*ng+i] = 0.0;
		f_uhg[13*ng+i] = 0.0;
		f_uhg[14*ng+i] = 0.0;
		f_uhg[15*ng+i] = 0.0;
		}
	});
}

void HdgFbou(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	if (ib == 1)
		HdgFbou1(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 2)
		HdgFbou2(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 3)
		HdgFbou3(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 4)
		HdgFbou4(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 5)
		HdgFbou5(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 6)
		HdgFbou6(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 7)
		HdgFbou7(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

