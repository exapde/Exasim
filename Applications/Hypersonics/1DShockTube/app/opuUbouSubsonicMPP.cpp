template <typename T> void opuUbou1(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
{

	for (int i = 0; i <ng; i++) {
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg4 = udg[3*ng+i];
		T udg5 = udg[4*ng+i];
		T udg6 = udg[5*ng+i];
		T udg7 = udg[6*ng+i];
		
		// T wdg6 = wdg[5*ng+i]; //SOMETHING WRONG HERE!
		// printf("opuubou1: wdg6 = %f\n", wdg6);
				
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
		//// Given Tinf, Yinf, uinf and pressure from solution, calculates rhoi, rhou, rhoE
		T Ucons[7] = {udg1, udg2, udg3, udg4, udg5, udg6, udg7};
		uinflow((double*) Ucons, (double) wdg6, (double*) param, (double*) uinf, mix);
		////
		f[0*ng+i] = Ucons[0];
		f[1*ng+i] = Ucons[1];
		f[2*ng+i] = Ucons[2];
		f[3*ng+i] = Ucons[3];
		f[4*ng+i] = Ucons[4];
		f[5*ng+i] = Ucons[5];
		f[6*ng+i] = Ucons[6];
		// for (int ii = 0; ii<7; ii++)
		// {
		// 	printf("ubou1: %f\n", f[ii*ng+i]);
		// }
	}
}

template <typename T> void opuUbou2(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
{
	for (int i = 0; i <ng; i++) {
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

		f[0*ng+i] = Ucons[0];
		f[1*ng+i] = Ucons[1];
		f[2*ng+i] = Ucons[2];
		f[3*ng+i] = Ucons[3];
		f[4*ng+i] = Ucons[4];
		f[5*ng+i] = Ucons[5];
		f[6*ng+i] = Ucons[6];
		// for (int ii = 0; ii<7; ii++)
		// {
		// 	printf("ubou2: %f\n", f[ii*ng+i]);
		// }
	}
}

template <typename T> void opuUbou(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
{
	if (ib == 1)
		opuUbou1(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw, mix);
	else if (ib == 2)
		opuUbou2(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw, mix);
}

template void opuUbou(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int, Mutation::Mixture *);
template void opuUbou(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int, Mutation::Mixture *);
