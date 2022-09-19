template <typename T> void opuSource(T* __restrict__ f, T* __restrict__ xdg, T* __restrict__ udg, T*__restrict__ odg, T*__restrict__ wdg, T*__restrict__ uinf, T*__restrict__ param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
{
		// std::cout << "begin source" << std::endl;
	for (int i = 0; i <ng; i++) {
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg4 = udg[3*ng+i];
		T udg5 = udg[4*ng+i];
		T udg6 = udg[5*ng+i];
		T udg7 = udg[6*ng+i];
		T t2 = 1.0/3.141592653589793;

		udg1 = udg1*(t2*atan(udg1*1.0E+3)+1.0/2.0)+3.183097800805168E-4;
		udg2 = udg2*(t2*atan(udg2*1.0E+3)+1.0/2.0)+3.183097800805168E-4;
		udg3 = udg3*(t2*atan(udg3*1.0E+3)+1.0/2.0)+3.183097800805168E-4;
		udg4 = udg4*(t2*atan(udg4*1.0E+3)+1.0/2.0)+3.183097800805168E-4;
		udg5 = udg5*(t2*atan(udg5*1.0E+3)+1.0/2.0)+3.183097800805168E-4;

		int nspecies = 5;
		double rho_scale = uinf[0];
		double u_scale = uinf[1];
		double rhoe_scale = uinf[2];
		double omega_scale = rho_scale*u_scale;

		double Ucons[7] = {udg1, udg2, udg3, udg4, udg5, udg6, udg7};
		double Ustate[6];
		dimensionalizeConsVars(Ucons, (double*)uinf, nspecies, 1);
		conservativeToState(Ucons, Ustate, (double*)uinf, nspecies);
		double rhovec[5] = {Ustate[0],Ustate[1],Ustate[2],Ustate[3],Ustate[4]};
		double rhoe = Ustate[nspecies];

		mix->setState(rhovec, &rhoe, 0);

		double wdot[5];
		mix->netProductionRates(wdot);

		f[0*ng+i] = wdot[0]/(omega_scale);
		f[1*ng+i] = wdot[1]/(omega_scale);
		f[2*ng+i] = wdot[2]/(omega_scale);
		f[3*ng+i] = wdot[3]/(omega_scale);
		f[4*ng+i] = wdot[4]/(omega_scale);
		f[5*ng+i] = 0.0;
		f[6*ng+i] = 0.0;
	}
		// std::cout << "end source " << std::endl;
}

template void opuSource(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, Mutation::Mixture *);
template void opuSource(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, Mutation::Mixture *);
