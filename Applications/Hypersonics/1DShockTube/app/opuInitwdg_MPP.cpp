template <typename T> void opuInitwdg(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne, Mutation::Mixture *mix)
{
	for (int i = 0; i <ng; i++) {
		int j = i%npe;
		int k = (i-j)/npe;
		
		T param1 = param[0];
		T param2 = param[1];
		T param3 = param[2];
		T param4 = param[3];
		T param5 = param[4];
		T param6 = param[5];
		T param7 = param[6];
		T param8 = param[7];

		T udg1 = param1;
		T udg2 = param2;
		T udg3 = param3;
		T udg4 = param4;
		T udg5 = param5;
		T udg6 = param6;
		T udg7 = param7;

		int nspecies = 5;
		T rho_inf = uinf[0];
		T u_inf = uinf[1];
		T rhoe_inf = rho_inf * u_inf * u_inf;
		double rhovec[5] = {udg1*rho_inf, udg2*rho_inf, udg3*rho_inf, udg4*rho_inf, udg5*rho_inf};

		// double rhovec[5] = {udg1, udg2, udg3, udg4, udg5};
		double wdot[5];
		double rhoe = abs(udg7-(udg6*udg6)*(udg1/2.0+udg2/2.0+udg3/2.0+udg4/2.0+udg5/2.0)*1.0/pow(udg1+udg2+udg3+udg4+udg5,2.0));
		// std::cout << rhoe << std::endl;
		// if (rhoe < 0){
			// std::cout << xdg1 << std::endl;
		// }
		rhoe = rhoe * rhoe_inf;
		// rhoe = fmax(rhoe,200.0);
		// std::cout << rhoe << std::endl;
		mix->setState(rhovec, &rhoe, 0);
		mix->netProductionRates(wdot);
		for (int ispecies = 0; ispecies < nspecies; ispecies++ )
		{
			f[j+npe*ispecies+npe*nce*k] = wdot[ispecies]/(rho_inf*u_inf);
		}
		f[j+npe*nspecies+npe*nce*k] = mix->P()/rhoe_inf;
		
	}
}

template void opuInitwdg(double *, double *, double *, double *, int, int, int, int, int, int, Mutation::Mixture *);
template void opuInitwdg(float *, float *, float *, float *, int, int, int, int, int, int, Mutation::Mixture *);
