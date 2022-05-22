template <typename T> void opuUbou1(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
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

		double Ucons[7] = {udg1, udg2, udg3, udg4, udg5, udg6, udg7};
		double Ustate[6];

		dimensionalizeConsVars(Ucons, (double*)uinf, nspecies, 1);
		conservativeToState(Ucons, Ustate, (double*)uinf, nspecies);
		double rhovec[5] = {Ustate[0],Ustate[1],Ustate[2],Ustate[3],Ustate[4]};
		double rhoe = Ustate[nspecies];

		mix->setState(rhovec, &rhoe, 0);

		T wdg6 = mix->P()/rhoe_scale;
		
		uinflow((double*) Ucons, (double) wdg6, (double*) param, (double*) uinf, mix);

		f[0*ng+i] = Ucons[0];
		f[1*ng+i] = Ucons[1];
		f[2*ng+i] = Ucons[2];
		f[3*ng+i] = Ucons[3];
		f[4*ng+i] = Ucons[4];
		f[5*ng+i] = Ucons[5];
		f[6*ng+i] = Ucons[6];

	}
}

template <typename T> void opuUbou2(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
{
	for (int i = 0; i <ng; i++) {
		T udg1 = fmax(udg[0*ng+i],0.0);
		T udg2 = fmax(udg[1*ng+i],0.0);
		T udg3 = fmax(udg[2*ng+i],0.0);
		T udg4 = fmax(udg[3*ng+i],0.0);
		T udg5 = fmax(udg[4*ng+i],0.0);
		T udg6 = udg[5*ng+i];
		T udg7 = udg[6*ng+i];
		// Given rho_i, rhou, modifies rhoE
		T Ucons[7] = {udg1, udg2, udg3, udg4, udg5, udg6, udg7}; 

		f[0*ng+i] = Ucons[0];
		f[1*ng+i] = Ucons[1];
		f[2*ng+i] = Ucons[2];
		f[3*ng+i] = Ucons[3];
		f[4*ng+i] = Ucons[4];
		f[5*ng+i] = Ucons[5];
		uoutflow((double*) Ucons, (double*) param, (double*) uinf, mix);

		f[6*ng+i] = Ucons[6];
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
