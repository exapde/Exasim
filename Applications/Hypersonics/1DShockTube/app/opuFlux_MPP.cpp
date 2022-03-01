template <typename T> void opuFlux(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
{
	for (int i = 0; i <ng; i++) {
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg4 = udg[3*ng+i];
		T udg5 = udg[4*ng+i];
		T udg6 = udg[5*ng+i];
		T udg7 = udg[6*ng+i];

		int nspecies = 5;
		T rho_inf = 0.0085;
		T u_inf = 542.0;
		double rhovec[5] = {udg1*rho_inf, udg2*rho_inf, udg3*rho_inf, udg4*rho_inf, udg5*rho_inf};
		double wdot[5];
		double rhoe = abs(udg7-(udg6*udg6)*(udg1/2.0+udg2/2.0+udg3/2.0+udg4/2.0+udg5/2.0)*1.0/pow(udg1+udg2+udg3+udg4+udg5,2.0));
		rhoe = rhoe * rho_inf * u_inf * u_inf;
		mix->setState(rhovec, &rhoe, 0);

		T wdg6 = mix->P();
		T t2 = udg1+udg2+udg3+udg4+udg5;
		T t3 = 1.0/t2;
		f[0*ng+i] = t3*udg1*udg6;
		f[1*ng+i] = t3*udg2*udg6;
		f[2*ng+i] = t3*udg3*udg6;
		f[3*ng+i] = t3*udg4*udg6;
		f[4*ng+i] = t3*udg5*udg6;
		f[5*ng+i] = wdg6+t3*(udg6*udg6);
		f[6*ng+i] = udg6*(t3*udg7+t3*wdg6);
	}
}

template void opuFlux(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, Mutation::Mixture *);
template void opuFlux(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, Mutation::Mixture *);
