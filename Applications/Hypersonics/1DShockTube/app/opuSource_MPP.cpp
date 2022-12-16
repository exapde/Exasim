template <typename T> void opuSource(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
{
	for (int i = 0; i <ng; i++) {
		int nspecies = 5;
		
		T udg1 = fmax(udg[0*ng+i],0.0);
		T udg2 = fmax(udg[1*ng+i],0.0);
		T udg3 = fmax(udg[2*ng+i],0.0);
		T udg4 = fmax(udg[3*ng+i],0.0);
		T udg5 = fmax(udg[4*ng+i],0.0);
		T udg6 = udg[5*ng+i];
		T udg7 = udg[6*ng+i];

		T rho_scale = uinf[0];
		T u_scale = uinf[1];
		T rhoe_scale = uinf[2];
		T omega_scale = rho_scale*u_scale;
		double Ucons[7] = {udg1, udg2, udg3, udg4, udg5, udg6, udg7};
		double Ustate[6];
		double wdot[5];

		dimensionalizeConsVars(Ucons, (double*)uinf, nspecies, 1);
		conservativeToState(Ucons, Ustate, (double*)uinf, nspecies);
		double rhovec[5] = {Ustate[0],Ustate[1],Ustate[2],Ustate[3],Ustate[4]};
		double rhoe = Ustate[nspecies];

		mix->setState(rhovec, &rhoe, 0);
		mix->netProductionRates(wdot);

		f[0*ng+i] = wdot[0]/(omega_scale);
		f[1*ng+i] = wdot[1]/(omega_scale);
		f[2*ng+i] = wdot[2]/(omega_scale);
		f[3*ng+i] = wdot[3]/(omega_scale);
		f[4*ng+i] = wdot[4]/(omega_scale);
		f[5*ng+i] = 0.0;
		f[6*ng+i] = 0.0;
	}
}

template void opuSource(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, Mutation::Mixture *);
template void opuSource(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, Mutation::Mixture *);
