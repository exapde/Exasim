template <typename T> void opuFbou1(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
{
	for (int i = 0; i <ng; i++) {
		T tau1 = tau[0];
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg4 = udg[3*ng+i];
		T udg5 = udg[4*ng+i];
		T udg6 = udg[5*ng+i];
		T udg7 = udg[6*ng+i];
		
		T rho_inf = uinf[0];
		T u_inf = uinf[1];
		T rhoe_inf = rho_inf * u_inf * u_inf;
		double rhovec[5] = {udg1*rho_inf, udg2*rho_inf, udg3*rho_inf, udg4*rho_inf, udg5*rho_inf};
		double rhoe = abs(udg7-(udg6*udg6)*(udg1/2.0+udg2/2.0+udg3/2.0+udg4/2.0+udg5/2.0)*1.0/pow(udg1+udg2+udg3+udg4+udg5,2.0));
		rhoe = rhoe * rhoe_inf;
		// printf("rho_inf: %f\n", rhoe_inf);
		mix->setState(rhovec, &rhoe, 0);
		// printf("opuubou1: wdg6 = %f; p = %f\n", wdg6, mix->P()/rhoe_inf);
		T wdg6 = mix->P()/rhoe_inf;
		// printf("opufbou1: wdg6 = %f\n", wdg6);
		//// Given Tinf, Yinf, uinf and pressure from solution, calculates rhoi, rhou, rhoE
		T Ucons[7] = {udg1, udg2, udg3, udg4, udg5, udg6, udg7}; 
		// for (int ii = 0; ii<7; ii++)
		// {
		// 	printf("fbou1 before uinflow: Ucons = %f\n", Ucons[ii*ng+i]);
		// }
		uinflow((double*) Ucons, (double) wdg6, (double*)param, (double*)uinf, mix);

		udg1 = Ucons[0];
		udg2 = Ucons[1];
		udg3 = Ucons[2];
		udg4 = Ucons[3];
		udg5 = Ucons[4];
		udg6 = Ucons[5];
		udg7 = Ucons[6];
		////
		// for (int ii = 0; ii<7; ii++)
		// {
		// 	printf("fbou1 after uinflow: Ucons = %f\n", Ucons[ii*ng+i]);
		// }

		T uhg1 = uhg[0*ng+i];
		T uhg2 = uhg[1*ng+i];
		T uhg3 = uhg[2*ng+i];
		T uhg4 = uhg[3*ng+i];
		T uhg5 = uhg[4*ng+i];
		T uhg6 = uhg[5*ng+i];
		T uhg7 = uhg[6*ng+i];
		T nlg1 = nlg[0*ng+i];
		T t2 = udg1+udg2+udg3+udg4+udg5;
		T t3 = 1.0/t2;
		f[0*ng+i] = tau1*(udg1-uhg1)+nlg1*t3*udg1*udg6;
		f[1*ng+i] = tau1*(udg2-uhg2)+nlg1*t3*udg2*udg6;
		f[2*ng+i] = tau1*(udg3-uhg3)+nlg1*t3*udg3*udg6;
		f[3*ng+i] = tau1*(udg4-uhg4)+nlg1*t3*udg4*udg6;
		f[4*ng+i] = tau1*(udg5-uhg5)+nlg1*t3*udg5*udg6;
		f[5*ng+i] = nlg1*(wdg6+t3*(udg6*udg6))+tau1*(udg6-uhg6);
		f[6*ng+i] = tau1*(udg7-uhg7)+nlg1*udg6*(t3*udg7+t3*wdg6);
		
		// for (int ii = 0; ii<7; ii++)
		// {
		// 	printf("fbou1: %f\n", f[ii*ng+i]);
		// }
		
	}
}

template <typename T> void opuFbou2(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
{
	for (int i = 0; i <ng; i++) {
		T tau1 = tau[0];
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg4 = udg[3*ng+i];
		T udg5 = udg[4*ng+i];
		T udg6 = udg[5*ng+i];
		T udg7 = udg[6*ng+i];
		// Given rho_i, rhou, modifies rhoE
		T Ucons[7] = {udg1, udg2, udg3, udg4, udg5, udg6, udg7}; 
		uoutflow((double*) Ucons, (double*) param, (double*) uinf, mix);

		udg1 = Ucons[0];
		udg2 = Ucons[1];
		udg3 = Ucons[2];
		udg4 = Ucons[3];
		udg5 = Ucons[4];
		udg6 = Ucons[5];
		udg7 = Ucons[6];

		T uhg1 = uhg[0*ng+i];
		T uhg2 = uhg[1*ng+i];
		T uhg3 = uhg[2*ng+i];
		T uhg4 = uhg[3*ng+i];
		T uhg5 = uhg[4*ng+i];
		T uhg6 = uhg[5*ng+i];
		T uhg7 = uhg[6*ng+i];
		T wdg6 = wdg[5*ng+i];
		T nlg1 = nlg[0*ng+i];
		T t2 = udg1+udg2+udg3+udg4+udg5;
		T t3 = 1.0/t2;
		f[0*ng+i] = tau1*(udg1-uhg1)+nlg1*t3*udg1*udg6;
		f[1*ng+i] = tau1*(udg2-uhg2)+nlg1*t3*udg2*udg6;
		f[2*ng+i] = tau1*(udg3-uhg3)+nlg1*t3*udg3*udg6;
		f[3*ng+i] = tau1*(udg4-uhg4)+nlg1*t3*udg4*udg6;
		f[4*ng+i] = tau1*(udg5-uhg5)+nlg1*t3*udg5*udg6;
		f[5*ng+i] = nlg1*(wdg6+t3*(udg6*udg6))+tau1*(udg6-uhg6);
		f[6*ng+i] = tau1*(udg7-uhg7)+nlg1*udg6*(t3*udg7+t3*wdg6);
		// for (int ii = 0; ii<7; ii++)
		// {
		// 	printf("fbou2: %f\n", f[ii*ng+i]);
		// }
		
	}
}

template <typename T> void opuFbou(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
{
	if (ib == 1)
		opuFbou1(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw, mix);
	else if (ib == 2)
		opuFbou2(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw, mix);
}

template void opuFbou(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int, Mutation::Mixture *);
template void opuFbou(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int, Mutation::Mixture *);
