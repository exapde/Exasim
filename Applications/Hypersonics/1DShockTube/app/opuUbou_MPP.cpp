template <typename T> void opuUbou1(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
{
	for (int i = 0; i <ng; i++) {
		T param1 = param[0];
		T param2 = param[1];
		T param3 = param[2];
		T param4 = param[3];
		T param5 = param[4];
		T param6 = param[5];
		T param7 = param[6];
		f[0*ng+i] = param1;
		f[1*ng+i] = param2;
		f[2*ng+i] = param3;
		f[3*ng+i] = param4;
		f[4*ng+i] = param5;
		f[5*ng+i] = param6;
		f[6*ng+i] = param7;
	}
}

template <typename T> void opuUbou2(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
{
	for (int i = 0; i <ng; i++) {
		T param8 = param[7];
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg4 = udg[3*ng+i];
		T udg5 = udg[4*ng+i];
		T udg6 = udg[5*ng+i];
		f[0*ng+i] = udg1;
		f[1*ng+i] = udg2;
		f[2*ng+i] = udg3;
		f[3*ng+i] = udg4;
		f[4*ng+i] = udg5;
		f[5*ng+i] = udg6;

        T t2 = udg1+udg2+udg3+udg4+udg5; // rho = sum(rho_i)
		T t3 = 1.0/t2; // rhoinv
        T rhou2 = t3*(udg6*udg6); // rho u^2 = (rhou)^2/rho
        
        // we have density and pressure...I suppose here we can fix density and temperature? 
        T rho_inf = uinf[0];
        T u_inf = uinf[1];
        T rhoe_inf = rho_inf * u_inf * u_inf;
        T rho_dim = t2*rho_inf;

        double rhovec[5] = {udg1*rho_inf, udg2*rho_inf, udg3*rho_inf, udg4*rho_inf, udg5*rho_inf};
        double Tout = 3943.7;

        mix->setState(rhovec, &Tout, 1);
        T rhoeDim = mix->mixtureEnergyMass() * rho_dim;
        
        T rhoe = rhoeDim / rhoe_inf;

        T rhoE = rhoe + rhou2;
        // std::cout << rhoE*rhoe_inf << std::endl;
        // std::cout << mix->T() << std::endl;
		f[6*ng+i] = rhoE;
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
